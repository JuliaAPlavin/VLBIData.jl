using PyCall
using Dates
using FITSIO: FITSHeader, FITS, read_header
using Unitful
using Parameters: @with_kw
using DataFrames: DataFrame
using AxisKeys

macro S_str(str) :(Symbol($str)) end


@with_kw struct FrequencyWindow
    freq::typeof(1f0u"Hz")
    width::typeof(1f0u"Hz")
    wavelen::typeof(1f0u"m") = u"c" / freq
    sideband::Int32
end

@with_kw struct UVHeader
    fits::FITSHeader
    object::String
    date_obs::Date
    stokes::Vector{Symbol}
    frequency::typeof(1.0u"Hz")
end

function UVHeader(fh::FITSHeader)
    @assert fh["CTYPE6"] == "RA" && fh["CTYPE7"] == "DEC"
    @assert fh["CTYPE4"] == "FREQ"
    @assert fh["CTYPE2"] == "COMPLEX" && fh["NAXIS2"] == 3

    val_to_stokes = Dict(-4=>:LR, -3=>:RL, -2=>:LL, -1=>:RR, 1=>:I, 2=>:Q, 3=>:U, 4=>:V)
    stokes = [val_to_stokes[val] for val in axis_vals(fh, "STOKES")]
    # @show fh
    date_obs = match(r"([\d-]+)(\(\d+\))?", fh["DATE-OBS"]).captures[1]

    return UVHeader(
        fits=fh,
        object=fh["OBJECT"],
        date_obs=Date(date_obs, dateformat"Y-m-d"),
        stokes=stokes,
        frequency=axis_dict(fh, "FREQ")["CRVAL"]*u"Hz",
    )
end

@with_kw struct UVData
    path::String
    header::UVHeader
    freq_windows::Vector{FrequencyWindow}
end

function read_freqs(uvh, fq_table)
    df = DataFrame(fq_table)
    @assert size(df, 1) == 1
    fq = Dict(pairs(df[1, :]))
    @assert !haskey(fq, :SIDEBAND) || all(fq[:SIDEBAND] .== 1)
    @assert all(fq[S"CH WIDTH"] .== fq[S"TOTAL BANDWIDTH"])
    res = ((a, b, c) -> FrequencyWindow(freq=a, width=b, sideband=c)).(uvh.frequency .+ fq[S"IF FREQ"].*u"Hz", fq[S"CH WIDTH"].*u"Hz", fq[:SIDEBAND])
    return isa(res, FrequencyWindow) ? [res] : res  # XXX: 1-sized array turn out as scalars
end

function read_data_raw(uvdata::UVData)
    py"""
    import numpy as np
    import astropy.io.fits

    def to_native_byteorder(arr):
        dt = arr.dtype.newbyteorder('=')
        return arr.astype(dt)

    with astropy.io.fits.open(open($(uvdata.path), 'rb'), memmap=False) as f:
        raw = f[0].data
    """
    raw = Dict(n => PyArray(py"to_native_byteorder(raw[$n])"o) for n in py"raw.dtype.names")
    return raw
end

function read_data_arrays(uvdata::UVData)
    raw = read_data_raw(uvdata)

    uvw_keys = [("UU", "VV", "WW"), ("UU--", "VV--", "WW--"), ("UU---SIN", "VV---SIN", "WW---SIN")]
    matching_keys = [keys for keys in uvw_keys if all(haskey(raw, k) for k in keys)]
    @assert length(matching_keys) == 1
    uvw_keys = first(matching_keys)

    count = size(raw["DATA"], 1)
    n_if = length(uvdata.freq_windows)
    axarr = KeyedArray(raw["DATA"],
        IX=1:count,
        DEC=[1], RA=[1],
        IF=1:n_if, FREQ=[1], STOKES=uvdata.header.stokes,
        COMPLEX=[:re, :im, :wt])
    axarr = axarr[DEC=1, RA=1, FREQ=1]  # drop always-singleton axes

    uvw = [kk => Float32.(raw[k]) .* (u"c" * u"s") .|> u"m"
        for (kk, k) in zip([:u, :v, :w], uvw_keys)]
    data = (
        uvw...,
        array_ix = (raw["BASELINE"] .% 1) .* 100 .|> x->round(Int8, x),
        ant1_ix = round.(Int, raw["BASELINE"]) .÷ 256 .|> Int8,
        ant2_ix = round.(Int, raw["BASELINE"]) .% 256 .|> Int8,
        # int_time = raw["INTTIM"] .* u"s",
        date = KeyedArray(Float64.(raw["DATE"]) .+ raw["_DATE"], IX=1:count),
        visibility = map((r, i) -> r + im * i, axarr(COMPLEX=:re), axarr(COMPLEX=:im)),
        weight = axarr(COMPLEX=:wt),
    )
    data = merge(data, (date_0 = Float32.(data.date .- minimum(data.date)),))
    return data
end

function read_data_dataframe(uvdata::UVData)
    data = read_data_arrays(uvdata)
    @assert ndims(data.visibility) == 3
    df = map(Iterators.product(axiskeys.(Ref(data.visibility), [:IX, :IF, :STOKES])...)) do (ix, iif, stokes)
        if_spec = uvdata.freq_windows[iif]
        uvw = (data.u[ix], data.v[ix], data.w[ix])
        uvw_wl = (uvw ./ if_spec.wavelen)
        (
            array_ix=data.array_ix[ix],
            ant1_ix=data.ant1_ix[ix],
            ant2_ix=data.ant2_ix[ix],
            date=data.date[ix],
            date_0=data.date_0[ix],
            stokes=stokes,
            iif=Int8(iif),
            if_spec=if_spec,
            u=uvw[1], v=uvw[2], w=uvw[3],
            u_wl=uvw_wl[1], v_wl=uvw_wl[2], w_wl=uvw_wl[3],
            visibility=data.visibility(IX=ix, IF=iif, STOKES=stokes),
            weight=data.weight(IX=ix, IF=iif, STOKES=stokes),
        )
    end
    df = filter(x -> x.weight > 0, df) |> DataFrame
    return df
end

function load(::Type{UVData}, path)
    fits = FITS(path)
    fh = read_header(fits[1])
    header = UVHeader(fh)
    freqs = read_freqs(header, fits["AIPS FQ"])
    close(fits)

    UVData(
        path=path,
        header=header,
        freq_windows=freqs,
    )
end
