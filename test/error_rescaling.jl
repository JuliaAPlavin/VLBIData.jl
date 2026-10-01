@testitem "exp_scale" begin
    using VLBIData: exp_scale, exp_scale!

    # m = 3: k = 2, and the 2nd of 3 iid Exp(1) has expectation 1/2 + 1/3
    @test exp_scale([1.0, 4.0, 16.0]) ≈ 4.8
    @test (@inferred exp_scale([1.0, 4.0, 16.0])) ≈ 4.8
    @test exp_scale([9.0, 9.0, 9.0]) ≈ 10.8
    # q = 1/3 → k = 1, expectation 1/3; q = 1 → k = 3, expectation 1 + 1/2 + 1/3
    @test exp_scale([1.0, 4.0, 16.0]; q=1/3) ≈ 3.0
    @test exp_scale([1.0, 4.0, 16.0]; q=1.0) ≈ 8.727272727272727
    # k rounds to 0
    @test_throws AssertionError exp_scale([1.0, 2.0, 3.0]; q=0.1)

    x = [16.0, 1.0, 4.0]
    @test exp_scale!(x) ≈ 4.8
    @test x[2] == 4.0
    @test sort(x) == [1.0, 4.0, 16.0]
end

@testitem "find_errmul ConsecutiveDifferencesStandard" begin
    using Dates, Uncertain
    using VLBIData.StructArrays
    using VLBIData: ConsecutiveDifferencesStandard, find_errmul

    # One baseline, one channel, equally spaced in time, σ = 1/√2 per point so every σ_ΔV = 1.
    # Successive |ΔV| are 1, 2, 4 ⇒ (|ΔV|/σ_ΔV)² = 1, 4, 16 ⇒ scale = 4/(1/2 + 1/3) = 4.8,
    # and the multiplier is √(4.8/2) = √2.4.
    uvtbl = StructArray([
        (datetime=DateTime(2020, 1, 1, 0, 0, 10i), freq_spec=123, stokes=:RR,
         spec=VisSpec(Baseline((:A, :B)), UV(0, 1)), value=Uncertain.Value(complex(v), 1/√2))
        for (i, v) in enumerate([0.0, 1.0, 3.0, 7.0])
    ])
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=3), uvtbl) ≈ 1.5491933384829668
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=4), uvtbl) |> isnothing
end

@testitem "compute_diffs and typical_Δt group by source and scan" begin
    using Dates, Unitful, Uncertain
    using VLBIData.StructArrays
    using VLBIData: compute_diffs, typical_Δt, intersect_nt_type, _compute_diffs_inner, _typical_Δt_inner

    # One baseline/channel, σ = 1/√2 per point so every σ_ΔV = 1. Values 0,1,3 | 10,12 | 0,4 in three
    # (source, scan) blocks 10 s apart throughout: within-block |ΔV| are 1,2 | 2 | 4, while the two
    # block boundaries would contribute 7 and 12 if the blocks were one group.
    row(t, v, src, sc) = (datetime=DateTime(2020, 1, 1) + Second(t), source=src, scan_id=sc,
                          freq_spec=123, stokes=:RR, spec=VisSpec(Baseline((:A, :B)), UV(0, 1)),
                          value=Uncertain.Value(complex(v), 1/√2))
    ts = 0:10:60
    vs = [0.0, 1.0, 3.0, 10.0, 12.0, 0.0, 4.0]
    srcs = [:A, :A, :A, :B, :B, :B, :B]
    scans = [1, 1, 1, 2, 2, 3, 3]
    tbl = StructArray(map(row, ts, vs, srcs, scans))

    diffs = compute_diffs(tbl; maxΔt=1u"minute")
    @test map(r -> r.abs, diffs) ≈ [1.0, 2.0, 2.0, 4.0]
    @test map(r -> r.err, diffs) ≈ [1.0, 1.0, 1.0, 1.0]

    # the same rows without the two key columns: every consecutive pair is taken
    tbl_plain = StructArray(map(r -> (; r.datetime, r.freq_spec, r.stokes, r.spec, r.value), tbl))
    @test map(r -> r.abs, compute_diffs(tbl_plain; maxΔt=1u"minute")) ≈ [1.0, 2.0, 7.0, 2.0, 12.0, 4.0]

    # scan 1 is sampled every 10 s, scan 2 every 60 s: the median of the per-scan gaps is 10 s,
    # while the 580 s gap between the scans would pull an ungrouped median to 35 s
    tsp = StructArray(map(row, [0, 10, 20, 600, 660], fill(1.0, 5), fill(:A, 5), [1, 1, 1, 2, 2]))
    @test typical_Δt(tsp) ≈ 10u"s"

    NT = intersect_nt_type(eltype(tbl), NamedTuple{(:source, :freq_spec, :stokes, :scan_id)})
    @test @inferred(_compute_diffs_inner(NT, tbl; maxΔt=1u"minute")) == diffs
    @test @inferred(_typical_Δt_inner(NT, tsp)) ≈ 10u"s"
end
