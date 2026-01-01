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

Statistics.mean(xs::AbstractVector{<:VLBI.FrequencyWindow}) = VLBI.FrequencyWindow(
	(@p xs map(_.freq) mean),
	(@p xs map(_.width) sum),
	(@p xs map(_.nchan) sum),
	(@p xs map(_.sideband) uniqueonly),
)

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
    xyz::SVector{3, Float64}
end

function Antenna(hdu_row::NamedTuple)
    if !isempty(hdu_row.ORBPARM) && hdu_row.ORBPARM != 0
        @warn "Antennas with ORBPARM detected, be careful" hdu_row.ORBPARM hdu_row.ANNAME
    end
    Antenna(; name=Symbol(hdu_row.ANNAME), id=hdu_row.NOSTA, xyz=hdu_row.STABXYZ)
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
    antennas = map(Antenna, hdu |> columntable |> StructArray)
    AntArray(;
        name=header["ARRNAM"],
        freq=strfloat_to_float(header["FREQ"]) * u"Hz",
        antennas,
    )
end

Base.length(a::AntArray) = length(a.antennas)
Base.getindex(a::AntArray, i::Int) = a.antennas[i]

Base.@kwdef struct UVData
    path::String
    header::Union{UVHeader,Nothing}
    freq_windows::Vector{FrequencyWindow}
    ant_arrays::Vector{AntArray}
end

function read_freqs(uvh, fq_table)
    fq_row = fq_table |> columntable |> StructArray |> only
    # fq_row = fq_row[Symbol.(["IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"])]
    nrows = @p fq_row values() filter(_ isa AbstractVector) (isempty(__) ? 1 : length(__[1]))
    fq_row = map(fq_row) do x
        isa(x, AbstractArray) ? x : fill(x, nrows)
    end
    ref_freq = @oget frequency(uvh) read_header(fq_table)["REF_FREQ"]*u"Hz"
    res = map(fq_row |> rowtable) do r
        total_bw = @oget r[S"TOTAL BANDWIDTH"] r[S"TOTAL_BANDWIDTH"]
        ch_width = @oget r[S"CH WIDTH"] r[S"CH_WIDTH"]
        curfreq = @oget r[S"IF FREQ"] r[S"BANDFREQ"]
        nchan = Int(total_bw / ch_width)
        FrequencyWindow(;
            freq=ref_freq + curfreq * u"Hz",
            width=total_bw * u"Hz",
            nchan,
            sideband=r.SIDEBAND,
        )
    end
end


function read_data_raw(uvdata::UVData, ::typeof(identity)=identity)
    fits = FITS(uvdata.path)
    if haskey(fits, "UV_DATA")
        return StructArray(fits["UV_DATA"] |> columntable)
    end
    hdu = GroupedHDU(fits.fitsfile, 1)
    read(hdu)
end

struct Baseline
    array_ix::Int8
    ant_ids::NTuple{2, Int8}
    ant_names::NTuple{2, Symbol}
end

function Base.getproperty(b::Baseline, key::Symbol)
    key == :ants_ix && (Base.depwarn("Baseline.ants_ix is deprecated, use ant_ids instead", :ants_ix); return b.ant_ids)
    return getfield(b, key)
end

@deprecate Baseline(array_ix::Integer, ant_ids::NTuple{2, Integer}) Baseline(array_ix, ant_ids, (Symbol(:ANT, ant_ids[1]), Symbol(:ANT, ant_ids[2])))

function Baseline(array_ix::Integer, ant_ids::NTuple{2, Integer}, ant_arrays::Vector{AntArray})
    ants = ant_arrays[array_ix].antennas
    names = map(ant_ids) do id
        ix = findfirst(ant -> ant.id == id, ants)
        if !isnothing(ix)
            ants[ix].name
        else
            @warn "Antenna index out of bounds, assigning generated name" length(ants) ant_ids
            Symbol(:ANT, id)
        end
    end
    Baseline(array_ix, ant_ids, names)
end

antennas(b::Baseline) =
    map(b.ant_ids, b.ant_names) do id, name
        Antenna(name, id, SVector(NaN, NaN, NaN))
    end
Accessors.set(b::Baseline, ::typeof(antennas), ants) = setproperties(b, ant_ids=map(a -> a.id, ants), ant_names=map(a -> a.name, ants))
@accessor antennas(x) = antennas(VLBI.Baseline(x))

antennas(x::Antenna) = (x,)
antennas(x::NTuple{<:Any,Antenna}) = x
antennas(x::AbstractVector{Antenna}) = x

abstract type AbstractSpec end

struct VisSpec{TUV<:UV} <: AbstractSpec
	bl::VLBI.Baseline
	uv::TUV
end
VisSpec(bl::VLBI.Baseline, uv::AbstractVector) = VisSpec(bl, UV(uv))

@accessor VLBI.Baseline(bl::Baseline) = identity(bl)
@accessor VLBI.Baseline(vs::VisSpec) = vs.bl
VLBI.Baseline(x::NamedTuple) = VLBI.Baseline(@oget x.baseline x.spec)

@accessor UV(x::VisSpec) = x.uv
@accessor UV(x::NamedTuple) = UV(x.spec)

AccessorsExtra.hasoptic(::Baseline, ::Type{UV}) = false
AccessorsExtra.hasoptic(::Union{Antenna,NTuple{<:Any,Antenna},AbstractVector{Antenna}}, ::Type{UV}) = false
AccessorsExtra.hasoptic(::UV, ::Type{Baseline}) = false
AccessorsExtra.hasoptic(obj, ::typeof(antennas)) = AccessorsExtra.hasoptic(obj, Baseline)

InterferometricModels.visibility(x::NamedTuple) = x.value

frequency(x::NamedTuple) = frequency(x.if_spec)


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
        Baseline(round(Int, (b % 1) * 100) + 1, (bi ÷ 256, bi % 256), uvdata.ant_arrays)
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
    Base.depwarn("table(::UVData) is deprecated, use uvtable(::UVData) instead", :table, force=true)
    _table(uvdata, impl)
end

function _table(uvdata::UVData, impl=identity)
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


uvtable(uvd::VLBI.UVData; stokes=(:I, :LL, :RR)) = @p uvd _table filter(_.stokes ∈ stokes) map((;
	_.datetime, _.stokes, _.if_ix, _.if_spec,
	spec=VisSpec(_.baseline, UV(_.uv)),
	value=U.Value(_.visibility, 1/√_.weight),
))

function load(::Type{UVData}, path)
    path = abspath(path)  # for RFC.File
    fits = FITS(path)
    fh = read_header(fits[1])
    header = try
        UVHeader(fh)
    catch e
        e isa KeyError || rethrow()
        nothing
    end
    freq_windows = read_freqs(header, @oget fits["AIPS FQ"] fits["FREQUENCY"])
    ant_arrays = AntArray[]
    for i in Iterators.countfrom(1)
        hdu = try
            haskey(fits, "AIPS AN") ? fits["AIPS AN", i] : fits["ARRAY_GEOMETRY", i]
        catch err
            err isa FITSIO.CFITSIO.CFITSIOError && "illegal HDU number" == err.errmsgshort && break
            rethrow()
        end
        push!(ant_arrays, AntArray(hdu))
    end
    close(fits)

    UVData(; path, header, freq_windows, ant_arrays)
end
