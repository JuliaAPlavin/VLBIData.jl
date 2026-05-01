# Earth Rotation Angle (IERS 2003 definition, exact)
# Approximation: uses UTC ≈ UT1 (max error 0.9s → negligible for UV sim)
_earth_rotation_angle(jd_ut1::Real) = 2π * (0.7790572732640 + 1.00273781191135448 * (jd_ut1 - 2451545.0))

# ECEF baseline → UVW in meters
function _ecef_to_uvw(baseline_ecef::SVector{3}, H, δ)
    sin_h, cos_h = sincos(H)
    sin_d, cos_d = sincos(δ)
    Bx, By, Bz = baseline_ecef
    u =  sin_h * Bx + cos_h * By
    v = -sin_d * cos_h * Bx + sin_d * sin_h * By + cos_d * Bz
    w =  cos_d * cos_h * Bx - cos_d * sin_h * By + sin_d * Bz
    UVW(u, v, w)
end

function _source_elevation(xyz_ecef::SVector{3}, H, dec)
    xyz_rot = _rotate_z(xyz_ecef, H)
    source_dir = SVector(cos(dec), zero(dec), sin(dec))
    cos_zenith = dot(xyz_rot, source_dir) / norm(xyz_rot)
    asin(clamp(cos_zenith, -1, 1))
end

function _rotate_z(xyz::SVector{3}, angle)
    s, c = sincos(angle)
    SVector(c * xyz[1] - s * xyz[2], s * xyz[1] + c * xyz[2], xyz[3])
end


function synthesize_uvtable(
    antennas::AbstractVector{<:Antenna},
    times::AbstractVector{<:DateTime},
    coords::AbstractSkyCoords;
    frequency::Unitful.Frequency,
    min_elevation = 0.0,
)
    src = ICRSCoords(coords)
    ra, dec = src.ra, src.dec
    λ = ustrip(u"m", u"c" / frequency)
    ant_pairs = [(antennas[i], antennas[j]) for i in eachindex(antennas) for j in i+1:lastindex(antennas)]

    flatmap(times) do t
        jd = julian_day(t)
        H = _earth_rotation_angle(jd) - ra

        filtermap(ant_pairs) do (ant1, ant2)
            if _source_elevation(ant1.xyz, H, dec) < min_elevation ||
               _source_elevation(ant2.xyz, H, dec) < min_elevation
                return nothing
            end

            uvw_m = _ecef_to_uvw(ant1.xyz - ant2.xyz, H, dec)
            uv_wl = UV(uvw_m[1] / λ, uvw_m[2] / λ)

            (
                datetime = t,
                spec = VisSpec(Baseline((ant1.name, ant2.name)), uv_wl),
                freq_spec = frequency,
            )
        end
    end |> StructArray
end

function synthesize_uvtable(
    model,
    antennas::AbstractVector{<:Antenna},
    times::AbstractVector{<:DateTime},
    coords::AbstractSkyCoords;
    kwargs...
)
    uvtbl = synthesize_uvtable(antennas, times, coords; kwargs...)
    @insert uvtbl.value = visibilities(model, map(UV, uvtbl))
end
