Base.@kwdef struct Antenna
    name::Symbol
    id::Int
    xyz::SVector{3, Float64}
end
@batteries Antenna eq=false

Antenna(name::Symbol, id::Int) = Antenna(name, id, SVector(NaN, NaN, NaN))

antennas(x::Antenna) = (x,)
antennas(x::NTuple{<:Any,Antenna}) = x
antennas(x::AbstractVector{Antenna}) = x

AccessorsExtra.hasoptic(::Union{Antenna,NTuple{<:Any,Antenna},AbstractVector{Antenna}}, ::Type{UV}) = false

Base.show(io::IO, a::Antenna) = print(io, "Antenna #", a.id, ": ", a.name)

Base.isless(a::Antenna, b::Antenna) = isless(a.id, b.id)
