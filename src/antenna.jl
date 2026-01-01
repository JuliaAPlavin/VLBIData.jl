@enum AntennaMountType begin
    # codes 0-6 must correspond those from the uvfits standard:
    AltAzimuth = 0
    Equatorial = 1
    Orbiting = 2
    XY = 3
    NaismithR = 4
    NaismithL = 5
    ApertureArray = 6

    # everything else:
    Unknown = 9999
end

Base.@kwdef struct Antenna
    name::Symbol
    xyz::SVector{3, Float64} = SVector(NaN, NaN, NaN)
    mount_type::AntennaMountType = AntennaMountType.Unknown
    poltypes::NTuple{2, Symbol} = (:U, :U)
end
@batteries Antenna eq=false

Antenna(name::Symbol; kwargs...) = Antenna(; name, kwargs...)

antenna_names(a::Antenna) = (a.name,)
antenna_names(a::Union{NTuple{<:Any,Antenna}, AbstractVector{Antenna}}) = map(only∘antenna_names, a)

AccessorsExtra.hasoptic(::Union{Antenna,NTuple{<:Any,Antenna},AbstractVector{Antenna}}, ::Type{UV}) = false

Base.show(io::IO, a::Antenna) = print(io, "Antenna ", a.name)

Base.isless(a::Antenna, b::Antenna) = isless(a.name, b.name)
