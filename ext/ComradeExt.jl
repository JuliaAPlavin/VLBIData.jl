module ComradeExt

using VLBIData
using VLBIData.StructArrays
using VLBIData.DataManipulation
using VLBIData.Unitful
using VLBIData.UnitfulAngles
using VLBIData.Uncertain
using VLBIData.IntervalSets
using VLBIData.InterferometricModels: visibility
using VLBIData.DateFormats
import Comrade

function Comrade.extract_table(uvtbl::AbstractVector{<:NamedTuple}; antennas, obscoords, objname)
    @assert allunique(a.name for a in antennas)
    tarr = map(antennas) do a
        (;
            sites = a.name,
            X     = a.xyz[1],
            Y     = a.xyz[2],
            Z     = a.xyz[3],
            SEFD1 = 0.0,
            SEFD2 = 0.0,
            fr_parallactic = 0.0,
            fr_elevation   = 0.0,
            fr_offset      = 0.0,
        )
    end |> StructArray

    scan_ivs = VLBI.scan_intervals(VLBI.GapBasedScans(), uvtbl)
    d0 = minimum(r -> r.datetime, uvtbl)
    scans = @p scan_ivs map((
        start=(leftendpoint(_) - d0) /ₜ Hour - √eps(),
        stop =(rightendpoint(_) - d0) /ₜ Hour + √eps(),
    )) StructArray

    datatable = map(uvtbl) do r
        Comrade.EHTArrayBaselineDatum(
            Float64(r.spec.uv[1]),
            Float64(r.spec.uv[2]),
            (r.datetime - d0) /ₜ Hour,
            Float64(ustrip(r.freq_spec.freq)),
            Tuple(r.spec.bl.antennas),
            (Comrade.CirBasis(), Comrade.CirBasis()),
            (0.0, 0.0),  # elevation placeholder
            (0.0, 0.0),  # parallactic angle placeholder
        )
    end |> StructArray

    bandwidth = @p uvtbl map(_.freq_spec) unique sum(_.width) ustrip(u"Hz", __)
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
end

function VLBI.uvtable(obs::Comrade.EHTObservationTable{<:Comrade.EHTVisibilityDatum})
    cfg = Comrade.arrayconfig(obs)
    dt = cfg.datatable
    meas = obs.measurement
    nse = obs.noise
    d0 = DateTime(modified_julian_day(cfg.mjd))

    StructArray(map(eachindex(meas)) do i
        (;
            datetime = d0 + Millisecond(round(Int64, dt.Ti[i] * 3_600_000)),
            spec = VisSpec(Baseline(dt.sites[i]), UV(dt.U[i], dt.V[i])),
            freq_spec = dt.Fr[i] * u"Hz",
            value = meas[i] ±ᵤ nse[i],
        )
    end)
end

end
