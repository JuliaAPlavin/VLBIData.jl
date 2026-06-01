@testitem "vlbiskymodels" begin
    using ComradeBase
    using VLBISkyModels
    using VLBIData
    using VLBIData.StructArrays

    # Create a simple test model
    img = zeros(32, 32)
    img[16, 16] = 1.0  # Point source at center
    gim = RectiGrid((X = range(-1, 1, length=32), Y = range(-1, 1, length=32)))
    model = ContinuousImage(img, gim, Gaussian())

    # Create test UV coordinates
    uvtbl_spec = [UV(0.1, 0.2), UV(-0.3, 0.4), UV(0.5, -0.1)] |> StructArray

    # Test the visibilitymap function
    result = visibilitymap(model, uvtbl_spec)

    # Check that we get a uvtable back with correct structure
    @test result isa StructArray
    @test length(result) == 3
    @test result.spec == [UV(0.1, 0.2), UV(-0.3, 0.4), UV(0.5, -0.1)]

    result_shifted = visibilitymap(shifted(model, 1., 0), uvtbl_spec)
    @test result_shifted isa StructArray
    @test length(result_shifted) == 3
    @test result_shifted.spec == [UV(0.1, 0.2), UV(-0.3, 0.4), UV(0.5, -0.1)]
    @test log.(result_shifted.value ./ result.value) ≈ [10/16, -30/16, 50/16]*im  rtol=1e-2

    # Test with freq_spec propagation
    using Unitful
    uvtbl_spec_freq = [(; spec=UV(0.1, 0.2), freq_spec=230e9u"Hz"), (; spec=UV(-0.3, 0.4), freq_spec=230e9u"Hz"), (; spec=UV(0.5, -0.1), freq_spec=230e9u"Hz")] |> StructArray
    result_freq = visibilitymap(model, uvtbl_spec_freq)
    @test result_freq isa StructArray
    @test length(result_freq) == 3
    @test result_freq.spec == [UV(0.1, 0.2), UV(-0.3, 0.4), UV(0.5, -0.1)]
    @test result_freq.freq_spec ≈ fill(230.0u"GHz", 3)
    @test result_freq.value ≈ result.value
end

@testitem "comrade extract_table: Stokes I round-trip" begin
    import Comrade
    using StaticArrays, Unitful, UnitfulAngles, Dates, SkyCoords
    using StructArrays
    using Uncertain
    using VLBIData.PolT: CirBasis

    # Two-station synthetic observation, hardcoded reference geometry.
    # Mixed mount types so fr_coefs branches are exercised.
    ant_FD = Antenna(
        name = :FD, xyz = SVector(-1324009.337937, -5332181.945962, 3231962.404011),
        mount_type = VLBIData.AntennaMountType.AltAzimuth,
        feed_offsets = (deg2rad(0.0), deg2rad(90.0)),
    )
    ant_LA = Antenna(
        name = :LA, xyz = SVector(-1449752.516, -4975298.376, 3709123.520),
        mount_type = VLBIData.AntennaMountType.NaismithR,
        feed_offsets = (deg2rad(45.0), deg2rad(135.0)),
    )
    antennas = [ant_FD, ant_LA]
    source = ICRSCoords(deg2rad(187.70595), deg2rad(12.39112))   # M87
    t0 = DateTime(2010, 12, 24, 8, 0, 0)
    spec_a = VisSpec(Baseline((:FD, :LA)), UV(1.234e6, 5.678e6))
    spec_b = VisSpec(Baseline((:FD, :LA)), UV(1.3e6, 5.7e6))
    fs = (freq=230e9u"Hz", width=2e9u"Hz")

    uvtbl = StructArray([
        (datetime=t0,             spec=spec_a, freq_spec=fs, stokes=:I, value=(1.5 + 0.0im) ±ᵤ 0.05),
        (datetime=t0 + Minute(1), spec=spec_b, freq_spec=fs, stokes=:I, value=(0.9 + 0.4im) ±ᵤ 0.15),
    ])

    obs = Comrade.extract_table(uvtbl; antennas, obscoords=source, objname="M87")

    @test obs isa Comrade.EHTObservationTable{<:Comrade.EHTVisibilityDatum{:I}}
    @test length(obs) == 2
    @test obs.config.source == :M87
    @test length(obs.config.tarr) == 2

    # tarr: per-station AIPS feed-rotation coefficients from mount_type
    tarr = obs.config.tarr
    fd_idx = findfirst(==(:FD), tarr.sites)
    la_idx = findfirst(==(:LA), tarr.sites)
    # AltAzimuth → (fr_elevation=0, fr_parallactic=1)
    @test tarr.fr_elevation[fd_idx] == 0.0
    @test tarr.fr_parallactic[fd_idx] == 1.0
    # NaismithR → (fr_elevation=1, fr_parallactic=1)
    @test tarr.fr_elevation[la_idx] == 1.0
    @test tarr.fr_parallactic[la_idx] == 1.0
    # fr_offset comes from feed_offsets[1] (POLAA in radians)
    @test tarr.fr_offset[fd_idx] == 0.0
    @test tarr.fr_offset[la_idx] ≈ deg2rad(45.0)

    # Per-row elevation and parallactic — hardcoded values (geocentric, ECEF)
    dt = obs.config.datatable
    @test rad2deg(dt.elevation[1][1]) ≈ 15.588400777523546 atol=1e-10
    @test rad2deg(dt.parallactic[1][1]) ≈ -61.467356702839666 atol=1e-10
    @test rad2deg(dt.elevation[1][2]) ≈ 14.143091412547628 atol=1e-10
    @test rad2deg(dt.parallactic[1][2]) ≈ -56.0016591605084 atol=1e-10

    # Polbasis present (irrelevant to Stokes I fitting but required by schema)
    @test dt.polbasis[1] == (CirBasis(), CirBasis())

    # round-trip: EHTObservationTable → uvtable
    rt = VLBI.uvtable(obs)
    @test length(rt) == 2
    @test issetequal(keys(rt[1]), (:datetime, :spec, :freq_spec, :value))
    @test rt[1].spec == spec_a
    @test rt[2].spec == spec_b
    @test U.value(rt[1].value) == 1.5 + 0.0im
    @test U.uncertainty(rt[1].value) == 0.05
    @test U.value(rt[2].value) == 0.9 + 0.4im
    @test U.uncertainty(rt[2].value) == 0.15
    @test rt[1].freq_spec ≈ 230e9u"Hz"
    @test rt[1].datetime == t0
    @test rt[2].datetime == t0 + Minute(1)
end

@testitem "comrade extract_table: coherency round-trip" begin
    import Comrade
    using StaticArrays, Unitful, UnitfulAngles, Dates, SkyCoords
    using StructArrays
    using Uncertain
    using VLBIData.PolT: CirBasis

    antennas = [
        Antenna(name=:FD, xyz=SVector(-1324009.337937, -5332181.945962, 3231962.404011),
                mount_type=VLBIData.AntennaMountType.AltAzimuth, feed_offsets=(0.0, deg2rad(90.0))),
        Antenna(name=:LA, xyz=SVector(-1449752.516, -4975298.376, 3709123.520),
                mount_type=VLBIData.AntennaMountType.AltAzimuth, feed_offsets=(0.0, deg2rad(90.0))),
    ]
    source = ICRSCoords(deg2rad(187.70595), deg2rad(12.39112))
    t = DateTime(2010, 12, 24, 8, 0, 0)
    spec = VisSpec(Baseline((:FD, :LA)), UV(1.234e6, 5.678e6))
    fs = (freq=230e9u"Hz", width=2e9u"Hz")

    uvtbl = StructArray([
        (datetime=t, spec, freq_spec=fs, stokes=:RR, value=(1.5 + 0.1im) ±ᵤ 0.05),
        (datetime=t, spec, freq_spec=fs, stokes=:LR, value=(0.2 - 0.05im) ±ᵤ 0.05),
        (datetime=t, spec, freq_spec=fs, stokes=:RL, value=(0.2 + 0.05im) ±ᵤ 0.05),
        (datetime=t, spec, freq_spec=fs, stokes=:LL, value=(1.4 - 0.05im) ±ᵤ 0.05),
    ])
    obs = Comrade.extract_table(uvtbl; antennas, obscoords=source, objname="M87")

    @test obs isa Comrade.EHTObservationTable{<:Comrade.EHTCoherencyDatum}
    @test length(obs) == 1
    @test obs.config.datatable.polbasis[1] == (CirBasis(), CirBasis())

    # Measurement is 2x2 SMatrix [RR RL; LR LL]
    M = obs.measurement[1]
    @test M isa SMatrix{2, 2, ComplexF64}
    @test M[1, 1] == 1.5 + 0.1im
    @test M[2, 1] == 0.2 - 0.05im
    @test M[1, 2] == 0.2 + 0.05im
    @test M[2, 2] == 1.4 - 0.05im
    N = obs.noise[1]
    @test N isa SMatrix{2, 2, Float64}
    @test all(==(0.05), N)

    # Reverse path: 2x2 → 4 rows with stokes labels
    rt = VLBI.uvtable(obs)
    @test length(rt) == 4
    rows_by_stokes = Dict(r.stokes => r for r in rt)
    @test issetequal(keys(rows_by_stokes), [:RR, :LR, :RL, :LL])
    @test rows_by_stokes[:RR].value == (1.5 + 0.1im) ±ᵤ 0.05
    @test rows_by_stokes[:LR].value == (0.2 - 0.05im) ±ᵤ 0.05
    @test rows_by_stokes[:RL].value == (0.2 + 0.05im) ±ᵤ 0.05
    @test rows_by_stokes[:LL].value == (1.4 - 0.05im) ±ᵤ 0.05
    @test rows_by_stokes[:RR].datetime == t
    @test rows_by_stokes[:RR].spec.bl == Baseline((:FD, :LA))
end

@testitem "comrade extract_table: coherency with partial polarization" begin
    import Comrade
    using StaticArrays, Unitful, UnitfulAngles, Dates, SkyCoords
    using StructArrays
    using Uncertain
    using VLBIData.PolT: CirBasis

    antennas = [
        Antenna(name=:FD, xyz=SVector(-1324009.337937, -5332181.945962, 3231962.404011),
                mount_type=VLBIData.AntennaMountType.AltAzimuth, feed_offsets=(0.0, deg2rad(90.0))),
        Antenna(name=:LA, xyz=SVector(-1449752.516, -4975298.376, 3709123.520),
                mount_type=VLBIData.AntennaMountType.AltAzimuth, feed_offsets=(0.0, deg2rad(90.0))),
    ]
    source = ICRSCoords(deg2rad(187.70595), deg2rad(12.39112))
    t1 = DateTime(2010, 12, 24, 8, 0, 0)
    t2 = t1 + Minute(5)
    spec = VisSpec(Baseline((:FD, :LA)), UV(1.234e6, 5.678e6))
    fs = (freq=230e9u"Hz", width=2e9u"Hz")

    # t1 group: only parallel hands. t2 group: complete 4. Globally stokes_set == _CIRCULAR_FEEDS,
    # so the dispatcher routes to the coherency path; the partial t1 group flows through with NaNs.
    uvtbl = StructArray([
        (datetime=t1, spec, freq_spec=fs, stokes=:RR, value=(1.5 + 0.1im) ±ᵤ 0.05),
        (datetime=t1, spec, freq_spec=fs, stokes=:LL, value=(1.4 - 0.05im) ±ᵤ 0.05),
        (datetime=t2, spec, freq_spec=fs, stokes=:RR, value=(2.5 + 0.1im) ±ᵤ 0.05),
        (datetime=t2, spec, freq_spec=fs, stokes=:LR, value=(0.3 - 0.05im) ±ᵤ 0.05),
        (datetime=t2, spec, freq_spec=fs, stokes=:RL, value=(0.3 + 0.05im) ±ᵤ 0.05),
        (datetime=t2, spec, freq_spec=fs, stokes=:LL, value=(2.4 - 0.05im) ±ᵤ 0.05),
    ])
    obs = Comrade.extract_table(uvtbl; antennas, obscoords=source, objname="M87")

    @test obs isa Comrade.EHTObservationTable{<:Comrade.EHTCoherencyDatum}
    @test length(obs) == 2
    @test obs.config.datatable.polbasis[1] == (CirBasis(), CirBasis())

    # Partial group (t1): finite parallel hands, NaN cross-hands in both measurement and noise.
    M1 = obs.measurement[1]
    @test M1 isa SMatrix{2, 2, ComplexF64}
    @test M1[1, 1] == 1.5 + 0.1im   # :RR
    @test M1[2, 2] == 1.4 - 0.05im  # :LL
    @test isnan(real(M1[2, 1]))     # :LR slot — NaN
    @test isnan(real(M1[1, 2]))     # :RL slot — NaN
    N1 = obs.noise[1]
    @test N1 isa SMatrix{2, 2, Float64}
    @test N1[1, 1] == 0.05
    @test N1[2, 2] == 0.05
    @test N1[2, 1] == 0
    @test N1[1, 2] == 0

    # Complete group (t2): all four finite, no NaNs.
    M2 = obs.measurement[2]
    @test M2[1, 1] == 2.5 + 0.1im
    @test M2[2, 1] == 0.3 - 0.05im
    @test M2[1, 2] == 0.3 + 0.05im
    @test M2[2, 2] == 2.4 - 0.05im
    @test all(==(0.05), obs.noise[2])
end

@testitem "comrade extract_table: error on unsupported stokes" begin
    import Comrade
    using StaticArrays, Unitful, Dates, SkyCoords, StructArrays
    using Uncertain

    ant1 = Antenna(name=:A1, xyz=SVector(6.378e6, 0.0, 0.0),
                   mount_type=VLBIData.AntennaMountType.AltAzimuth, feed_offsets=(0.0, deg2rad(90.0)))
    ant2 = Antenna(name=:A2, xyz=SVector(0.0, 6.378e6, 0.0),
                   mount_type=VLBIData.AntennaMountType.AltAzimuth, feed_offsets=(0.0, deg2rad(90.0)))
    source = ICRSCoords(0.0, 0.0)
    t = DateTime(2020, 1, 1)
    spec = VisSpec(Baseline((:A1, :A2)), UV(1e6, 1e6))
    fs = (freq=230e9u"Hz", width=2e9u"Hz")

    # Only parallel hands → error (partial coherency)
    uvtbl_partial = StructArray([
        (datetime=t, spec, freq_spec=fs, stokes=:RR, value=(1.0 + 0.0im) ±ᵤ 0.1),
        (datetime=t, spec, freq_spec=fs, stokes=:LL, value=(1.0 + 0.0im) ±ᵤ 0.1),
    ])
    @test_throws Exception Comrade.extract_table(uvtbl_partial; antennas=[ant1, ant2], obscoords=source)

    # Linear feed products → error (this bridge is circular-only for now)
    uvtbl_linear = StructArray([
        (datetime=t, spec, freq_spec=fs, stokes=:XX, value=(1.0 + 0.0im) ±ᵤ 0.1),
        (datetime=t, spec, freq_spec=fs, stokes=:YX, value=(0.1 + 0.0im) ±ᵤ 0.1),
        (datetime=t, spec, freq_spec=fs, stokes=:XY, value=(0.1 + 0.0im) ±ᵤ 0.1),
        (datetime=t, spec, freq_spec=fs, stokes=:YY, value=(1.0 + 0.0im) ±ᵤ 0.1),
    ])
    @test_throws Exception Comrade.extract_table(uvtbl_linear; antennas=[ant1, ant2], obscoords=source)
end

@testitem "comrade extract_table: scan_id" begin
    import Comrade
    using Dates
    using Unitful
    using Uncertain

    antennas = [
        Antenna(:AA, xyz=(1e6, 2e6, 3e6)),
        Antenna(:BB, xyz=(4e6, 5e6, 6e6)),
    ]

    d0 = DateTime(2017, 4, 11, 0, 0, 0)
    # Two rows close in time (no gap) but with different scan_ids.
    # GapBasedScans would put them in the same scan; pre-set scan_id splits them.
    uvtbl = [
        (
            datetime = d0,
            spec = VisSpec(Baseline((:AA, :BB)), UV(1e9, 2e9)),
            freq_spec = (freq=230e9u"Hz", width=2e9u"Hz"),
            value = (1.0 + 0.5im) ±ᵤ0.1,
            scan_id = 1,
        ),
        (
            datetime = d0 + Second(10),
            spec = VisSpec(Baseline((:AA, :BB)), UV(1.1e9, 2.1e9)),
            freq_spec = (freq=230e9u"Hz", width=2e9u"Hz"),
            value = (0.9 + 0.4im) ±ᵤ0.15,
            scan_id = 2,
        ),
    ]

    obs = Comrade.extract_table(uvtbl; antennas, obscoords=(ra=180.0u"°", dec=30.0u"°"))
    # Should have 2 scans (from scan_id), not 1 (from GapBasedScans)
    @test length(obs.config.scans) == 2
end
