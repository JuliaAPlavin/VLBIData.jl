@with_kw struct FrequencyWindow
    freq::typeof(1f0u"Hz")
    width::typeof(1f0u"Hz")
    nchan::Int
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

@with_kw struct Antenna
    name::Symbol
    id::Int
end

function Antenna(hdu_row)
    @assert isempty(hdu_row.ORBPARM)
    Antenna(; name=Symbol(hdu_row.ANNAME), id=hdu_row.NOSTA)
end

@with_kw struct AntArray
    name::String
    freq::typeof(1f0u"Hz")
    antennas::Vector{Antenna}
end

strfloat_to_float(x::AbstractFloat) = x
strfloat_to_float(x::String) = parse(Float64, replace(x, "D" => "E"))

function AntArray(hdu::TableHDU)
    header = read_header(hdu)
    antennas = map(Antenna, hdu |> columntable |> rowtable)
    AntArray(;
        name=header["ARRNAM"],
        freq=strfloat_to_float(header["FREQ"]) * u"Hz",
        antennas,
    )
end

@with_kw struct UVData
    path::String
    header::UVHeader
    freq_windows::Vector{FrequencyWindow}
    ant_arrays::Vector{AntArray}
end

function read_freqs(uvh, fq_table)
    fq_row = fq_table |> rowtable |> only
    fq_row = fq_row[Symbol.(["IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"])]
    fq_row = map(fq_row) do x
        isa(x, Real) ? [x] : x
    end
    res = map(fq_row |> rowtable) do r
        nchan = Int(r.var"TOTAL BANDWIDTH" / r.var"CH WIDTH")
        FrequencyWindow(;
            freq=uvh.frequency + r.var"IF FREQ" * u"Hz",
            width=r.var"TOTAL BANDWIDTH" * u"Hz",
            nchan,
            sideband=r.SIDEBAND,
        )
    end
end

function read_data_raw(uvdata::UVData)
    pyimport_conda("numpy", "numpy")
    pyimport_conda("astropy", "astropy")
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
    n_chan = map(fw -> fw.nchan, uvdata.freq_windows) |> unique |> only
    axarr = KeyedArray(raw["DATA"],
        IX=1:count,
        DEC=[1], RA=[1],
        IF=1:n_if,
        FREQ=1:n_chan,
        STOKES=uvdata.header.stokes,
        COMPLEX=[:re, :im, :wt])
    # drop always-singleton axes
    axarr = if size(axarr, :FREQ) == 1
        axarr[DEC=1, RA=1, FREQ=1]
    else
        axarr[DEC=1, RA=1]
    end

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

function read_data_table(uvdata::UVData)
    data = read_data_arrays(uvdata)
    @assert ndims(data.visibility) ∈ (3, 4)
    # `|> columntable |> rowtable` is faster than `|> rowtable` alone
    df = map(data.visibility |> columntable |> rowtable) do r
        ix = r.IX
        if_spec = uvdata.freq_windows[r.IF]
        uvw = (data.u[ix], data.v[ix], data.w[ix])
        uvw_wl = (uvw ./ if_spec.wavelen)
        (
            array_ix=data.array_ix[ix],
            ant1_ix=data.ant1_ix[ix],
            ant2_ix=data.ant2_ix[ix],
            date=data.date[ix],
            date_0=data.date_0[ix],
            stokes=r.STOKES,
            iif=Int8(r.IF),
            if_spec=if_spec,
            u=uvw[1], v=uvw[2], w=uvw[3],
            u_wl=uvw_wl[1], v_wl=uvw_wl[2], w_wl=uvw_wl[3],
            visibility=r.value,
            weight=data.weight(; Base.structdiff(r, NamedTuple{(:value,)})...),
        )
    end
    filter!(x -> x.weight > 0, df)
    return df
end

function load(::Type{UVData}, path)
    fits = FITS(path)
    fh = read_header(fits[1])
    header = UVHeader(fh)
    freq_windows = read_freqs(header, fits["AIPS FQ"])
    ant_arrays = AntArray[]
    for i in Iterators.countfrom(1)
        hdu = try
            fits["AIPS AN", i]
        catch err
            occursin("illegal HDU number", string(err)) && break
            rethrow()
        end
        push!(ant_arrays, AntArray(hdu))
    end
    close(fits)

    UVData(; path, header, freq_windows, ant_arrays)
end
