Base.@kwdef struct Antenna
    name::Symbol
    xyz::SVector{3, Float64} = SVector(NaN, NaN, NaN)
end
@batteries Antenna eq=false

Antenna(name::Symbol) = Antenna(; name)

antennas(x::Antenna) = (x,)
antennas(x::NTuple{<:Any,Antenna}) = x
antennas(x::AbstractVector{Antenna}) = x

AccessorsExtra.hasoptic(::Union{Antenna,NTuple{<:Any,Antenna},AbstractVector{Antenna}}, ::Type{UV}) = false

Base.show(io::IO, a::Antenna) = print(io, "Antenna ", a.name)

Base.isless(a::Antenna, b::Antenna) = isless(a.name, b.name)
