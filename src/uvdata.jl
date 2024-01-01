Base.@kwdef struct FrequencyWindow
    freq::typeof(1f0u"Hz")
    width::typeof(1f0u"Hz")
    nchan::Int16
    sideband::Int8
end

frequency(fw::FrequencyWindow, kind::Symbol=:reference) = if kind == :reference
	fw.freq
elseif kind == :average
	@assert fw.sideband == 1
	@assert fw.nchan > 0
	fw.freq + fw.width / 2
else
	error("Unsupported kind: $kind")
end

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


function read_data_raw(uvdata::UVData, ::typeof(identity)=identity)
    fits = FITS(uvdata.path)
    hdu = GroupedHDU(fits.fitsfile, 1)
    read(hdu)
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

UV(uvw::UVW) = UV(uvw.u, uvw.v)

struct Baseline
    array_ix::Int8
    ants_ix::NTuple{2, Int8}
end

function read_data_arrays(uvdata::UVData, impl=identity)
    raw = read_data_raw(uvdata, impl)

    uvw_keys = @p begin
        [(S"UU", S"VV", S"WW"), (S"UU--", S"VV--", S"WW--"), (S"UU---SIN", S"VV---SIN", S"WW---SIN")]
        filteronly(_ ⊆ keys(raw))
    end

    count = size(raw[:DATA])[end]
    n_if = length(uvdata.freq_windows)
    n_chan = map(fw -> fw.nchan, uvdata.freq_windows) |> unique |> only
    axarr = KeyedArray(raw[:DATA],
        COMPLEX=[:re, :im, :wt],
        STOKES=uvdata.header.stokes,
        FREQ=1:n_chan,
        IF=1:n_if,
        RA=[1], DEC=[1],
        _=1:count,
    )
    # drop always-singleton axes
    axarr = axarr[DEC=1, RA=1]

    uvw_m = UVW.([Float32.(raw[k]) .* (u"c" * u"s") .|> u"m" for k in uvw_keys]...)
    baseline = map(raw[:BASELINE]) do b
        bi = floor(Int, b)
        Baseline(round(Int, (b % 1) * 100) + 1, (bi ÷ 256, bi % 256))
    end
    data = (;
        uvw_m,
        baseline,
        datetime = julian_day.(Float64.(raw[:DATE]) .+ raw[:_DATE]),
        visibility = complex.(axarr(COMPLEX=:re), axarr(COMPLEX=:im)),
        weight = axarr(COMPLEX=:wt),
    )
    if haskey(raw, :INTTIM)
        data = merge(data, (int_time = raw[:INTTIM] .* u"s",))
    end
    return data
end

function table(uvdata::UVData, impl=identity)
    data = read_data_arrays(uvdata, impl)
    @assert ndims(data.visibility) == 4
    
    @p begin
        data.visibility
        # `|> columntable |> rowtable` is faster than `|> rowtable` alone
        map(__ |> columntable |> rowtable) do r
            ix = r._
            if_spec = uvdata.freq_windows[r.IF]
            uvw_m = data.uvw_m[ix]
            uvw_wl = UVW(ustrip.(Unitful.NoUnits, uvw_m ./ (u"c" / frequency(if_spec, :average))))
            (;
                baseline=data.baseline[ix],
                datetime=data.datetime[ix],
                stokes=r.STOKES,
                if_ix=Int8(r.IF),
                if_spec=if_spec,
                uv_m=UV(uvw_m), w_m=uvw_m.w,
                uv=UV(uvw_wl), w=uvw_wl.w,
                visibility=r.value,
            )
        end
        StructArray()
        @insert __.weight = collect(vec(data.weight))
        filter!(_.weight > 0)
    end
end

function load(::Type{UVData}, path)
    path = abspath(path)  # for RFC.File
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
