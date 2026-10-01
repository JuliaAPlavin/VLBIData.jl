@testitem "selfcal" begin
    using Uncertain
    using Dates
    using StaticArrays
    using InterferometricModels
    using VLBIData.StructArrays

    antennas = [:A, :B, :C, :D]
    bls = [(antennas[i], antennas[j]) for i in 1:4 for j in i+1:4]
    dt = DateTime(2020, 1, 1, 0, 0, 0)
    model_vis = Dict(
        (:A,:B) => 1.0+0.5im, (:A,:C) => 0.8-0.3im, (:A,:D) => 0.6+0.2im,
        (:B,:C) => 0.9+0.1im, (:B,:D) => 0.7-0.4im, (:C,:D) => 0.5+0.6im,
    )
    true_gains = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    make_uvtbl(gains) = StructArray([
        (datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
         freq_spec=230, stokes=:RR, source=:cal,
         value=U.Value(gains[bl[1]] * gains[bl[2]] * model_vis[bl], 0.01))
        for (idx, bl) in enumerate(bls)
    ])
    mvis = [complex(model_vis[bl]) for bl in bls]

    # --- Exact gain recovery ---
    result = VLBIData.selfcal(VLBIData.AmplitudeAnalytic(), mvis, make_uvtbl(true_gains))
    @test length(result.gains) == 4
    for g in result.gains
        @test g.gain ≈ true_gains[g.antenna] rtol=1e-10
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
         freq_spec=230, stokes=:RR, source=:cal,
         value=U.Value(complex(point_gains[bl[1]] * point_gains[bl[2]]), 0.01))
        for (idx, bl) in enumerate(bls)
    ])
    gains_scalar = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), Point(1.0, SVector(0.0, 0.0)), point_uvtbl)
    @test length(gains_scalar) == 4
    for g in gains_scalar
        @test g.gain ≈ point_gains[g.antenna] rtol=1e-10
    end

    # --- solve_gains + apply_gains roundtrip ---
    uvtbl_test = make_uvtbl(true_gains)
    gains_rt = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), mvis, uvtbl_test)
    cal_data = VLBIData.apply_gains(VLBIData.AmplitudeAnalytic(),gains_rt, uvtbl_test)
    for (i, bl) in enumerate(bls)
        @test U.value(cal_data.value[i]) ≈ model_vis[bl] rtol=1e-10
    end

    # --- Multiple time groups ---
    dt2 = DateTime(2020, 1, 1, 1, 0, 0)
    gains_t1 = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    gains_t2 = Dict(:A => 0.7, :B => 1.3, :C => 0.95, :D => 1.15)
    multi_dt_uvtbl = StructArray(vcat(
        [(datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:RR, source=:cal,
          value=U.Value(gains_t1[bl[1]] * gains_t1[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
        [(datetime=dt2, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:RR, source=:cal,
          value=U.Value(gains_t2[bl[1]] * gains_t2[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
    ))
    multi_dt_mvis = vcat(mvis, mvis)
    multi_dt_result = VLBIData.selfcal(VLBIData.AmplitudeAnalytic(), multi_dt_mvis, multi_dt_uvtbl)
    @test length(multi_dt_result.gains) == 8  # 4 antennas × 2 times
    for g in multi_dt_result.gains
        expected = g.datetime == dt ? gains_t1[g.antenna] : gains_t2[g.antenna]
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
          freq_spec=230, stokes=:RR, source=:cal,
          value=U.Value(gains_rr[bl[1]] * gains_rr[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
        [(datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
          freq_spec=230, stokes=:LL, source=:cal,
          value=U.Value(gains_ll[bl[1]] * gains_ll[bl[2]] * model_vis[bl], 0.01))
         for (idx, bl) in enumerate(bls)],
    ))
    multi_stokes_mvis = vcat(mvis, mvis)
    multi_stokes_gains = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), multi_stokes_mvis, multi_stokes_uvtbl)
    @test length(multi_stokes_gains) == 8  # 4 antennas × 2 stokes
    for g in multi_stokes_gains
        expected = g.stokes == :RR ? gains_rr[g.antenna] : gains_ll[g.antenna]
        @test g.gain ≈ expected rtol=1e-10
    end
    multi_stokes_result = VLBIData.selfcal(VLBIData.AmplitudeAnalytic(), multi_stokes_mvis, multi_stokes_uvtbl)
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
         freq_spec=230, stokes=:RR, source=:cal,
         value=U.Value(v, σ))
    end)
    corrupted_gains = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), mvis, corrupted_uvtbl)
    # Gains from non-corrupted baselines should still be approximately correct
    for g in corrupted_gains
        @test g.gain ≈ true_gains[g.antenna] rtol=0.1
    end

    # --- Vector-of-NamedTuple input ---
    vec_gains = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), mvis, collect(make_uvtbl(true_gains)))
    @test length(vec_gains) == 4
    for g in vec_gains
        @test g.gain ≈ true_gains[g.antenna] rtol=1e-10
    end

    # --- Idempotence: second pass gives gains ≈ 1.0 ---
    pass1 = VLBIData.selfcal(VLBIData.AmplitudeAnalytic(), mvis, make_uvtbl(true_gains))
    pass2_gains = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), mvis, pass1.data)
    for g in pass2_gains
        @test g.gain ≈ 1.0 rtol=1e-10
    end

    # --- Bipartite baseline graph (chain: A-B, B-C, no A-C) ---
    chain_bls = [(:A, :B), (:B, :C)]
    chain_gains = Dict(:A => 1.3, :B => 0.8, :C => 1.1)
    chain_model = Dict((:A,:B) => 1.0+0.5im, (:B,:C) => 0.9-0.2im)
    chain_uvtbl = StructArray([
        (datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
         freq_spec=230, stokes=:RR, source=:cal,
         value=U.Value(chain_gains[bl[1]] * chain_gains[bl[2]] * chain_model[bl], 0.01))
        for (idx, bl) in enumerate(chain_bls)
    ])
    chain_mvis = [complex(chain_model[bl]) for bl in chain_bls]
    chain_result = VLBIData.selfcal(VLBIData.AmplitudeAnalytic(), chain_mvis, chain_uvtbl)
    # Gain products on each baseline must match true products
    for (i, bl) in enumerate(chain_bls)
        g1 = only(g for g in chain_result.gains if g.antenna == bl[1]).gain
        g2 = only(g for g in chain_result.gains if g.antenna == bl[2]).gain
        @test g1 * g2 ≈ chain_gains[bl[1]] * chain_gains[bl[2]] rtol=1e-10
    end

    # --- Square bipartite graph (4-antenna cycle: A-B, B-C, C-D, D-A) ---
    cycle_bls = [(:A, :B), (:B, :C), (:C, :D), (:A, :D)]
    cycle_gains = Dict(:A => 1.2, :B => 0.9, :C => 1.1, :D => 0.8)
    cycle_model = Dict((:A,:B) => 1.0+0.5im, (:B,:C) => 0.9-0.2im, (:C,:D) => 0.7+0.3im, (:A,:D) => 0.6-0.1im)
    cycle_uvtbl = StructArray([
        (datetime=dt, spec=VisSpec(Baseline(bl), UV(Float64(idx), Float64(idx+1))),
         freq_spec=230, stokes=:RR, source=:cal,
         value=U.Value(cycle_gains[bl[1]] * cycle_gains[bl[2]] * cycle_model[bl], 0.01))
        for (idx, bl) in enumerate(cycle_bls)
    ])
    cycle_mvis = [complex(cycle_model[bl]) for bl in cycle_bls]
    cycle_result = VLBIData.selfcal(VLBIData.AmplitudeAnalytic(), cycle_mvis, cycle_uvtbl)
    for (i, bl) in enumerate(cycle_bls)
        g1 = only(g for g in cycle_result.gains if g.antenna == bl[1]).gain
        g2 = only(g for g in cycle_result.gains if g.antenna == bl[2]).gain
        @test g1 * g2 ≈ cycle_gains[bl[1]] * cycle_gains[bl[2]] rtol=1e-10
    end

    # --- apply_gains preserves fields ---
    uvtbl_orig = make_uvtbl(true_gains)
    gains_for_fields = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), mvis, uvtbl_orig)
    cal_for_fields = VLBIData.apply_gains(VLBIData.AmplitudeAnalytic(),gains_for_fields, uvtbl_orig)
    for (orig, cal) in zip(uvtbl_orig, cal_for_fields)
        @test cal.datetime == orig.datetime
        @test cal.spec == orig.spec
        @test cal.freq_spec == orig.freq_spec
        @test cal.stokes == orig.stokes
        @test cal.source == orig.source
    end
end

@testitem "solve_gains on a table without grouping keys" begin
    using Uncertain
    using VLBIData.StructArrays

    # no datetime/stokes/freq_spec: one solution over the whole table
    # unit model, |V| = 6, 15, 10 on A-B, B-C, A-C ⇒ g_A² = 6·10/15 = 4 ⇒ g = (2, 3, 5)
    tbl = StructArray([
        (spec=VisSpec(Baseline((:A, :B)), UV(1, 2)), value=U.Value(6.0im, 0.01)),
        (spec=VisSpec(Baseline((:B, :C)), UV(3, 4)), value=U.Value(-15.0 + 0im, 0.01)),
        (spec=VisSpec(Baseline((:A, :C)), UV(5, 6)), value=U.Value(10.0 + 0im, 0.01)),
    ])
    gains = VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), ones(ComplexF64, 3), tbl)
    @test gains isa StructArray
    @test propertynames(gains) == (:antenna, :gain)
    @test gains.antenna == [:A, :B, :C]
    @test gains.gain ≈ [2, 3, 5] rtol=1e-12
    @test VLBIData.solve_gains(VLBIData.AmplitudeAnalytic(), ones(ComplexF64, 3), collect(tbl)) == gains
end
