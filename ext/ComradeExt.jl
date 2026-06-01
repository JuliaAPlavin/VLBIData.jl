module ComradeExt

using VLBIData
using VLBIData.StructArrays
using VLBIData.DataManipulation
using VLBIData.Unitful
using VLBIData.UnitfulAngles
using VLBIData.Uncertain
using VLBIData.IntervalSets
import PolarizedTypes
const PolT = PolarizedTypes
using VLBIData.InterferometricModels: visibility
using VLBIData.DateFormats
using VLBIData.StaticArrays
import Comrade

# AIPS Memo 117 §4.1: per-mount feed-rotation coefficients (fr_elevation, fr_parallactic).
# fr_offset is per-station and comes from antenna.feed_offsets (POLAA in radians).
_fr_coefs(m::VLBIData.AntennaMountType.T) =
    m == VLBIData.AntennaMountType.AltAzimuth ? (0.0,  1.0) :
    m == VLBIData.AntennaMountType.Equatorial ? (0.0,  0.0) :
    m == VLBIData.AntennaMountType.XY         ? (0.0,  0.0) :
    m == VLBIData.AntennaMountType.NaismithR  ? (1.0,  1.0) :
    m == VLBIData.AntennaMountType.NaismithL  ? (-1.0, 1.0) :
    (0.0, 0.0)   # Orbiting / ApertureArray / Unknown

# Per-row (elevation, parallactic) tuple for the two antennas of a baseline at given time.
function _pointings(ant_xyz_by_name, baseline_antennas, t, ra_rad, dec_rad)
    H = VLBIData._earth_rotation_angle(julian_day(t)) - ra_rad
    xyz1 = ant_xyz_by_name[baseline_antennas[1]]
    xyz2 = ant_xyz_by_name[baseline_antennas[2]]
    el = (VLBIData._source_elevation(xyz1, H, dec_rad), VLBIData._source_elevation(xyz2, H, dec_rad))
    pa = (VLBIData._parallactic_angle(xyz1, H, dec_rad), VLBIData._parallactic_angle(xyz2, H, dec_rad))
    (el, pa)
end

# Circular-only: this bridge currently assumes circular feeds (R/L) at every station.
const _CIRCULAR_FEEDS = (:RR, :LR, :RL, :LL)

function Comrade.extract_table(uvtbl::AbstractVector{<:NamedTuple}; antennas, obscoords, objname="UNKNOWN")
    @assert allunique(a.name for a in antennas)

    tarr = map(antennas) do a
        fr_el, fr_par = _fr_coefs(a.mount_type)
        (;
            sites = a.name,
            X     = a.xyz[1],
            Y     = a.xyz[2],
            Z     = a.xyz[3],
            SEFD1 = 0.0,
            SEFD2 = 0.0,
            fr_parallactic = fr_par,
            fr_elevation   = fr_el,
            fr_offset      = a.feed_offsets[1],
        )
    end |> StructArray

    scan_ivs = VLBI.scan_intervals(
        hasproperty(StructArray(uvtbl), :scan_id) ? nothing : VLBI.GapBasedScans(),
        uvtbl)
    d0 = DateTime(Date(minimum(r -> r.datetime, uvtbl)))   # midnight of first row's date; Ti = hours since
    scans = @p scan_ivs map((
        start=(leftendpoint(_) - d0) /ₜ Hour - √eps(),
        stop =(rightendpoint(_) - d0) /ₜ Hour + √eps(),
    )) StructArray

    bandwidth = @p uvtbl map(_.freq_spec) unique map(_.width) uniqueonly() ustrip(u"Hz", __)

    ra_rad  = ustrip(u"rad", obscoords.ra)
    dec_rad = ustrip(u"rad", obscoords.dec)
    ant_xyz_by_name = Dict(a.name => a.xyz for a in antennas)

    stokes_set = hasproperty(StructArray(uvtbl), :stokes) ? Set(unique(uvtbl.stokes)) : Set([:I])

    if stokes_set == Set([:I])
        # Stokes-I path
        datatable = map(uvtbl) do r
            el, pa = _pointings(ant_xyz_by_name, Tuple(r.spec.bl.antennas), r.datetime, ra_rad, dec_rad)
            Comrade.EHTArrayBaselineDatum(
                Float64(r.spec.uv[1]),
                Float64(r.spec.uv[2]),
                (r.datetime - d0) /ₜ Hour,
                Float64(ustrip(u"Hz", r.freq_spec.freq)),
                Tuple(r.spec.bl.antennas),
                (PolT.CirBasis(), PolT.CirBasis()),
                el, pa,
            )
        end |> StructArray

        config = Comrade.EHTArrayConfiguration(
            Float64(bandwidth), tarr, scans, mjd(Date(d0)),
            ustrip(u"hourAngle", obscoords.ra), ustrip(u"°", obscoords.dec),
            Symbol(objname), :UTC, datatable
        )

        T = Comrade.EHTVisibilityDatum{:I, Float64, eltype(datatable)}
        return Comrade.EHTObservationTable{T}(
            map(U.value ∘ visibility, uvtbl),
            map(U.uncertainty ∘ visibility, uvtbl),
            config)

    elseif stokes_set == Set(_CIRCULAR_FEEDS)
        # Coherency path — circular feeds (R/L) only.
        cm_tbl = VLBI.uvtable_values_to(CoherencyMatrix, uvtbl)

        datatable = map(cm_tbl) do r
            el, pa = _pointings(ant_xyz_by_name, Tuple(r.spec.bl.antennas), r.datetime, ra_rad, dec_rad)
            Comrade.EHTArrayBaselineDatum(
                Float64(r.spec.uv[1]),
                Float64(r.spec.uv[2]),
                (r.datetime - d0) /ₜ Hour,
                Float64(ustrip(u"Hz", r.freq_spec.freq)),
                Tuple(r.spec.bl.antennas),
                (PolT.CirBasis(), PolT.CirBasis()),
                el, pa,
            )
        end |> StructArray

        meas = map(r -> U.value.(r.value), cm_tbl)
        nse = map(r -> U.uncertainty.(r.value), cm_tbl)

        config = Comrade.EHTArrayConfiguration(
            Float64(bandwidth), tarr, scans, mjd(Date(d0)),
            ustrip(u"hourAngle", obscoords.ra), ustrip(u"°", obscoords.dec),
            Symbol(objname), :UTC, datatable
        )

        T = Comrade.EHTCoherencyDatum{Float64, eltype(datatable),
                                       SMatrix{2, 2, ComplexF64, 4},
                                       SMatrix{2, 2, Float64, 4}}
        return Comrade.EHTObservationTable{T}(meas, nse, config)
    else
        error("ComradeExt: stokes set $stokes_set is not supported — only Stokes-I or the four circular feed products (:RR,:LR,:RL,:LL) are accepted. Linear or mixed feeds aren't handled yet.")
    end
end

function VLBI.uvtable(obs::Comrade.EHTObservationTable{<:Comrade.EHTVisibilityDatum})
    cfg = Comrade.arrayconfig(obs)
    dt = cfg.datatable
    meas = obs.measurement
    nse = obs.noise
    d0 = DateTime(mjd(cfg.mjd))

    StructArray(map(eachindex(meas)) do i
        (;
            datetime = d0 + dt.Ti[i] *ₜ Hour,
            spec = VisSpec(Baseline(dt.sites[i]), UV(dt.U[i], dt.V[i])),
            freq_spec = dt.Fr[i] * u"Hz",
            value = meas[i] ±ᵤ nse[i],
        )
    end)
end

function VLBI.uvtable(obs::Comrade.EHTObservationTable{<:Comrade.EHTCoherencyDatum})
    cfg = Comrade.arrayconfig(obs)
    dt = cfg.datatable
    meas = obs.measurement
    nse = obs.noise
    d0 = DateTime(mjd(cfg.mjd))

    # Always circular per the forward path. Matrix layout per src/polarizedtypes.jl:
    # (1,1)→:RR, (2,1)→:LR, (1,2)→:RL, (2,2)→:LL
    rows = mapreduce(vcat, eachindex(meas)) do i
        M = meas[i]
        N = nse[i]
        time = d0 + dt.Ti[i] *ₜ Hour
        spec = VisSpec(Baseline(dt.sites[i]), UV(dt.U[i], dt.V[i]))
        fs = dt.Fr[i] * u"Hz"
        [(datetime = time, spec = spec, freq_spec = fs, stokes = stk,
          value  = M[ii, jj] ±ᵤ N[ii, jj])
         for (stk, (ii, jj)) in zip(_CIRCULAR_FEEDS, ((1,1), (2,1), (1,2), (2,2)))]
    end
    StructArray(rows)
end

end
