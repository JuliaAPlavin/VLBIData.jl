using Dates, FITSIO, Unitful, Parameters, AxisArrays, DataFrames, Query, PyCall, SplitApplyCombine



function DataFrames.DataFrame(tbl::TableHDU)
    pairs = [Symbol(n) => read(tbl, n) for n in FITSIO.colnames(tbl)]
    pairs = [
        k => (ndims(v) == 1 ? v : [v[:, i] for i in 1:size(v, 2)])
        for (k, v) in pairs]
    return DataFrame(pairs...)
end

Base.Dict(row::DataFrameRow) = Dict(zip(string.(keys(row)), values(row)))



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
    fq = Dict(df[1, :])
    @assert all(fq["SIDEBAND"] .== 1) && all(fq["CH WIDTH"] .== fq["TOTAL BANDWIDTH"])
    res = ((a, b, c) -> FrequencyWindow(freq=a, width=b, sideband=c)).(uvh.frequency .+ fq["IF FREQ"].*u"Hz", fq["CH WIDTH"].*u"Hz", fq["SIDEBAND"])
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
    axarr = AxisArray(raw["DATA"],
        Axis{:IX}(1:count),
        Axis{:DEC}([1]), Axis{:RA}([1]),
        Axis{:IF}(1:n_if), Axis{:FREQ}([1]), Axis{:STOKES}(uvdata.header.stokes),
        Axis{:COMPLEX}([:re, :im, :wt]))
    axarr = axarr[Axis{:DEC}(1), Axis{:RA}(1), Axis{:FREQ}(1)]  # drop always-singleton axes

    uvw = [kk => Float32.(raw[k]) .* (u"c" * u"s") .|> u"m"
        for (kk, k) in zip([:u, :v, :w], uvw_keys)]
    data = (
        uvw...,
        array_ix = (raw["BASELINE"] .% 1) .* 100 .|> x->round(Int8, x),
        ant1_ix = round.(Int, raw["BASELINE"]) .÷ 256 .|> Int8,
        ant2_ix = round.(Int, raw["BASELINE"]) .% 256 .|> Int8,
        # int_time = raw["INTTIM"] .* u"s",
        date = AxisArray(Float64.(raw["DATE"]) .+ raw["_DATE"], Axis{:IX}(1:count)),
        visibility = map((r, i) -> r + im * i, axarr[Axis{:COMPLEX}(:re)], axarr[Axis{:COMPLEX}(:im)]),
        weight = axarr[Axis{:COMPLEX}(:wt)],
    )
    data = merge(data, (date_0 = Float32.(data.date .- minimum(data.date)),))
    return data
end

function read_data_dataframe(uvdata::UVData)
    data = read_data_arrays(uvdata)
    @assert ndims(data.visibility) == 3
    df = map(Iterators.product(AxisArrays.axes.(Ref(data.visibility), [Axis{:IX}, Axis{:IF}, Axis{:STOKES}])...)) do (ix, iif, stokes)
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
            visibility=data.visibility[Axis{:IX}(ix), Axis{:IF}(iif), Axis{:STOKES}(stokes)],
            weight=data.weight[Axis{:IX}(ix), Axis{:IF}(iif), Axis{:STOKES}(stokes)],
        )
    end
    df = df |> @filter(_.weight > 0) |> DataFrame
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
