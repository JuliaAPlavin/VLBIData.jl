@testitem "selfcal" begin
    using Uncertain
    using Dates
    using StaticArrays
    using InterferometricModels
    using VLBIData.StructArrays

    stations = [:A, :B, :C, :D]
    bls = [(stations[i], stations[j]) for i in 1:4 for j in i+1:4]
    dt = DateTime(2020, 1, 1, 0, 0, 0)
    model_vis = Dict(
        (:A,:B) => 1.0+0.5im, (:A,:C) => 0.8-0.3im, (:A,:D) => 0.6+0.2im,
        (:B,:C) => 0.9+0.1im, (:B,:D) => 0.7-0.4im, (:C,:D) => 0.5+0.6im,
    )
    true_gains = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    make_uvtbl(gains) = StructArray([
        (datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
         freq_spec=230, stokes=:RR, source="cal",
         value=U.Value(gains[bl[1]] * gains[bl[2]] * model_vis[bl], 0.01))
        for (idx, bl) in enumerate(bls)
    ])
    mvis = [complex(model_vis[bl]) for bl in bls]

    # --- Exact gain recovery ---
    result = VLBI.selfcal(VLBI.AmplitudeAnalytic(), mvis, make_uvtbl(true_gains))
    @test length(result.gains) == 4
    for g in result.gains
        @test g.gain ≈ true_gains[g.station] rtol=1e-10
    end

    # --- Calibrated data matches model ---
    for (i, bl) in enumerate(bls)
        @test U.value(result.data.value[i]) ≈ model_vis[bl] rtol=1e-10
    end

    # --- Uncertainty scales correctly ---
    for (i, bl) in enumerate(bls)
        g1, g2 = true_gains[bl[1]], true_gains[bl[2]]
        @test U.uncertainty(result.data.value[i]) ≈ 0.01 / (g1 * g2) rtol=1e-10
    end

    # --- Point source model (unresolved calibrator) ---
    point_gains = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    point_uvtbl = StructArray([
        (datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
         freq_spec=230, stokes=:RR, source="cal",
         value=U.Value(complex(point_gains[bl[1]] * point_gains[bl[2]]), 0.01))
        for (idx, bl) in enumerate(bls)
    ])
    gains_scalar = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), Point(1.0, SVector(0.0, 0.0)), point_uvtbl)
    @test length(gains_scalar) == 4
    for g in gains_scalar
        @test g.gain ≈ point_gains[g.station] rtol=1e-10
    end

    # --- solve_gains + apply_gains roundtrip ---
    uvtbl_test = make_uvtbl(true_gains)
    gains_rt = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mvis, uvtbl_test)
    cal_data = VLBI.apply_gains(VLBI.AmplitudeAnalytic(),gains_rt, uvtbl_test)
    for (i, bl) in enumerate(bls)
        @test U.value(cal_data.value[i]) ≈ model_vis[bl] rtol=1e-10
    end

    # --- Multiple time groups ---
    dt2 = DateTime(2020, 1, 1, 1, 0, 0)
    gains_t1 = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    gains_t2 = Dict(:A => 0.7, :B => 1.3, :C => 0.95, :D => 1.15)
    multi_dt_uvtbl = StructArray(vcat(
        [(datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:RR, source="cal",
          value=U.Value(gains_t1[bl[1]] * gains_t1[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
        [(datetime=dt2, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:RR, source="cal",
          value=U.Value(gains_t2[bl[1]] * gains_t2[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
    ))
    multi_dt_mvis = vcat(mvis, mvis)
    multi_dt_result = VLBI.selfcal(VLBI.AmplitudeAnalytic(), multi_dt_mvis, multi_dt_uvtbl)
    @test length(multi_dt_result.gains) == 8  # 4 stations × 2 times
    for g in multi_dt_result.gains
        expected = g.datetime == dt ? gains_t1[g.station] : gains_t2[g.station]
        @test g.gain ≈ expected rtol=1e-10
    end
    for (row, mv) in zip(multi_dt_result.data, multi_dt_mvis)
        @test U.value(row.value) ≈ mv rtol=1e-10
    end

    # --- Multiple stokes groups ---
    gains_rr = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    gains_ll = Dict(:A => 1.05, :B => 1.15, :C => 0.85, :D => 0.95)
    multi_stokes_uvtbl = StructArray(vcat(
        [(datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:RR, source="cal",
          value=U.Value(gains_rr[bl[1]] * gains_rr[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
        [(datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:LL, source="cal",
          value=U.Value(gains_ll[bl[1]] * gains_ll[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
    ))
    multi_stokes_mvis = vcat(mvis, mvis)
    multi_stokes_gains = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), multi_stokes_mvis, multi_stokes_uvtbl)
    @test length(multi_stokes_gains) == 8  # 4 stations × 2 stokes
    for g in multi_stokes_gains
        expected = g.stokes == :RR ? gains_rr[g.station] : gains_ll[g.station]
        @test g.gain ≈ expected rtol=1e-10
    end
    multi_stokes_result = VLBI.selfcal(VLBI.AmplitudeAnalytic(), multi_stokes_mvis, multi_stokes_uvtbl)
    for (row, mv) in zip(multi_stokes_result.data, multi_stokes_mvis)
        @test U.value(row.value) ≈ mv rtol=1e-10
    end

    # --- Non-uniform weights ---
    # One baseline has wildly corrupted value and huge uncertainty.
    # WLS should still recover gains accurately due to low weight on that baseline.
    corrupted_uvtbl = StructArray(map(enumerate(bls)) do (idx, bl)
        σ = bl == (:A, :B) ? 1e4 : 0.01
        v = bl == (:A, :B) ? 100.0+0im : true_gains[bl[1]] * true_gains[bl[2]] * model_vis[bl]
        (datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
         freq_spec=230, stokes=:RR, source="cal",
         value=U.Value(v, σ))
    end)
    corrupted_gains = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mvis, corrupted_uvtbl)
    # Gains from non-corrupted baselines should still be approximately correct
    for g in corrupted_gains
        @test g.gain ≈ true_gains[g.station] rtol=0.1
    end

    # --- Vector-of-NamedTuple input ---
    vec_gains = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mvis, collect(make_uvtbl(true_gains)))
    @test length(vec_gains) == 4
    for g in vec_gains
        @test g.gain ≈ true_gains[g.station] rtol=1e-10
    end

    # --- Idempotence: second pass gives gains ≈ 1.0 ---
    pass1 = VLBI.selfcal(VLBI.AmplitudeAnalytic(), mvis, make_uvtbl(true_gains))
    pass2_gains = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mvis, pass1.data)
    for g in pass2_gains
        @test g.gain ≈ 1.0 rtol=1e-10
    end

    # --- apply_gains preserves fields ---
    uvtbl_orig = make_uvtbl(true_gains)
    gains_for_fields = VLBI.solve_gains(VLBI.AmplitudeAnalytic(), mvis, uvtbl_orig)
    cal_for_fields = VLBI.apply_gains(VLBI.AmplitudeAnalytic(),gains_for_fields, uvtbl_orig)
    for (orig, cal) in zip(uvtbl_orig, cal_for_fields)
        @test cal.datetime == orig.datetime
        @test cal.spec == orig.spec
        @test cal.freq_spec == orig.freq_spec
        @test cal.stokes == orig.stokes
        @test cal.source == orig.source
    end
end
