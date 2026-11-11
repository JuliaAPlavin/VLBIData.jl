@enumx AntennaMountType begin
    # 0-6 are the FITS MNTSTA codes, converted by value
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
    mount_type::AntennaMountType.T = AntennaMountType.Unknown
    poltypes::NTuple{2, Symbol} = (:U, :U)
    feed_offsets::NTuple{2, Float64} = (NaN, NaN)
end
@batteries Antenna eq=false

Antenna(name::Symbol; kwargs...) = Antenna(; name, kwargs...)

antenna_names(a::Antenna) = (a.name,)
antenna_names(a::Union{NTuple{<:Any,Antenna}, AbstractVector{Antenna}}) = map(only∘antenna_names, a)

AccessorsExtra.hasoptic(::Union{Antenna,NTuple{<:Any,Antenna},AbstractVector{Antenna}}, ::Type{UV}) = false

Base.show(io::IO, a::Antenna) = print(io, "Antenna ", a.name)

Base.isless(a::Antenna, b::Antenna) = isless(a.name, b.name)

"""
    lonlatalt(a::Antenna)

Geodetic position from the antenna's ECEF `xyz` (WGS-84), in bare SI: `lon` and `lat` in **radians**,
`alt` (height above the ellipsoid) in **meters**. Single-pass Bowring form; the altitude expression is
division-free, so it holds at the poles too.
"""
@stable function lonlatalt(a::Antenna)
    A = 6378137.0; f = 1 / 298.257223563                       # WGS-84 semi-major axis [m], flattening
    B = A * (1 - f); e2 = f * (2 - f); ep2 = e2 / (1 - e2)     # semi-minor, first/second eccentricity²
    x, y, z = a.xyz
    p = hypot(x, y); θ = atan(z * A, p * B)
    lon = atan(y, x)
    lat = atan(z + ep2 * B * sin(θ)^3, p - e2 * A * cos(θ)^3)
    alt = p * cos(lat) + z * sin(lat) - A * sqrt(1 - e2 * sin(lat)^2)
    (; lon, lat, alt)
end
