Base.@kwdef struct FrequencyWindow
    freq::typeof(1f0u"Hz")
    width::typeof(1f0u"Hz")
    nchan::Int16
    sideband::Int8
end

frequency(fw::FrequencyWindow) = fw.freq

Base.@kwdef struct UVHeader
    fits::FITSHeader
    object::String
    date_obs::Date
    stokes::Vector{Symbol}
    frequency::typeof(1.0u"Hz")
end

frequency(h::UVHeader) = h.frequency
Dates.Date(h::UVHeader) = h.date_obs

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

Base.@kwdef struct Antenna
    name::Symbol
    id::Int
end

function Antenna(hdu_row)
    @assert isempty(hdu_row.ORBPARM)
    Antenna(; name=Symbol(hdu_row.ANNAME), id=hdu_row.NOSTA)
end
Base.@kwdef struct AntArray
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

Base.@kwdef struct UVData
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
            freq=frequency(uvh) + r.var"IF FREQ" * u"Hz",
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

struct UVW{T} <: FieldVector{3, T}
    u::T
    v::T
    w::T
end

struct UV{T} <: FieldVector{2, T}
    u::T
    v::T
end

struct Baseline
    array_ix::Int8
    ants_ix::NTuple{2, Int8}
end

function read_data_arrays(uvdata::UVData)
    raw = read_data_raw(uvdata)

    uvw_keys = @p begin
        [("UU", "VV", "WW"), ("UU--", "VV--", "WW--"), ("UU---SIN", "VV---SIN", "WW---SIN")]
        filter(_ ⊆ keys(raw))
        only()
    end

    count = size(raw["DATA"], 1)
    n_if = length(uvdata.freq_windows)
    n_chan = map(fw -> fw.nchan, uvdata.freq_windows) |> unique |> only
    axarr = KeyedArray(raw["DATA"],
        _=1:count,
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

    uvw_m = UVW.([Float32.(raw[k]) .* (u"c" * u"s") .|> u"m" for k in uvw_keys]...)
    baseline = map(raw["BASELINE"]) do b
        bi = floor(Int, b)
        Baseline((b % 1) * 100 + 1, (bi ÷ 256, bi % 256))
    end
    data = (;
        uvw_m,
        baseline,
        datetime = julian_day.(Float64.(raw["DATE"]) .+ raw["_DATE"]),
        visibility = complex.(axarr(COMPLEX=:re), axarr(COMPLEX=:im)),
        weight = axarr(COMPLEX=:wt),
    )
    if haskey(raw, "INTTIM")
        data = merge(data, (int_time = raw["INTTIM"] .* u"s",))
    end
    return data
end

function Tables.rows(uvdata::UVData)
    data = read_data_arrays(uvdata)
    @assert ndims(data.visibility) ∈ (3, 4)
    # `|> columntable |> rowtable` is faster than `|> rowtable` alone
    df = map(data.visibility |> columntable |> rowtable) do r
        ix = r._
        if_spec = uvdata.freq_windows[r.IF]
        uvw_m = data.uvw_m[ix]
        uvw_wl = ustrip.(Unitful.NoUnits, uvw_m ./ (u"c" / frequency(if_spec)))
        (
            baseline=data.baseline[ix],
            datetime=data.datetime[ix],
            stokes=r.STOKES,
            if_ix=Int8(r.IF),
            if_spec=if_spec,
            uv_m=UV(uvw_m[1:2]), w_m=uvw_m[3],
            uv=UV(uvw_wl[1:2]), w=uvw_wl[3],
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
