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

    # values 0, 1, 3, 7, 8 ⇒ pair powers 1, 4, 16, 1; the gaps' median is 10 s ⇒ the cap is 30 s
    capped(ts) = StructArray([
        (datetime=DateTime(2020, 1, 1) + Second(t), freq_spec=123, stokes=:RR,
         spec=VisSpec(Baseline((:A, :B)), UV(0, 1)), value=Uncertain.Value(complex(v), 1/√2))
        for (t, v) in zip(ts, [0.0, 1.0, 3.0, 7.0, 8.0])
    ])
    # 70 s gap dropped: powers 1, 4, 16 ⇒ √2.4 as above
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=3), capped([0, 10, 20, 30, 100])) ≈ 1.5491933384829668
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=4), capped([0, 10, 20, 30, 100])) |> isnothing
    # 30 s gap kept: powers 1, 1, 4, 16, k = 2 ⇒ scale = 1/(1/3 + 1/4) = 12/7, multiplier √(6/7)
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=4), capped([0, 10, 20, 30, 60])) ≈ 0.9258200997725514
end

@testitem "consecutive_pairs and typical_Δt group by source and scan" begin
    using Dates, Unitful, Uncertain
    using VLBIData.StructArrays
    using VLBIData: ConsecutiveDifferencesStandard, consecutive_pairs, typical_Δt, find_errmul

    # One baseline/channel, σ = 1/√2 per point so every σ_ΔV = 1. Values 0,1,3 | 10,12 | 0,4 in three
    # (source, scan) blocks 10 s apart throughout.
    row(t, v, src, sc) = (datetime=DateTime(2020, 1, 1) + Second(t), source=src, scan_id=sc,
                          freq_spec=123, stokes=:RR, spec=VisSpec(Baseline((:A, :B)), UV(0, 1)),
                          value=Uncertain.Value(complex(v), 1/√2))
    ts = 0:10:60
    vs = [0.0, 1.0, 3.0, 10.0, 12.0, 0.0, 4.0]
    srcs = [:A, :A, :A, :B, :B, :B, :B]
    scans = [1, 1, 1, 2, 2, 3, 3]
    tbl = StructArray(map(row, ts, vs, srcs, scans))

    pairs = @inferred consecutive_pairs(tbl)
    @test pairs.i == [1, 2, 4, 6]
    @test pairs.j == [2, 3, 5, 7]
    # powers 1, 4, 4, 16, k = 2 ⇒ scale = 4/(1/3 + 1/4) = 48/7, multiplier √(24/7)
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=4), tbl) ≈ 1.8516401995451028

    # the same rows without the two key columns: every consecutive pair is taken
    tbl_plain = StructArray(map(r -> (; r.datetime, r.freq_spec, r.stokes, r.spec, r.value), tbl))
    pairs_plain = @inferred consecutive_pairs(tbl_plain)
    @test pairs_plain.i == 1:6
    @test pairs_plain.j == 2:7
    # powers 1, 4, 49, 4, 144, 16, k = 3 ⇒ scale = 4/(1/4 + 1/5 + 1/6) = 240/37, multiplier √(120/37)
    @test find_errmul(ConsecutiveDifferencesStandard(; min_cnt=6), tbl_plain) ≈ 1.8009006755629926

    # scan 1 is sampled every 10 s, scan 2 every 60 s: the median of the per-scan gaps is 10 s,
    # while the 580 s gap between the scans would pull an ungrouped median to 35 s
    tsp = StructArray(map(row, [0, 10, 20, 600, 660], fill(1.0, 5), fill(:A, 5), [1, 1, 1, 2, 2]))
    @test @inferred(typical_Δt(tsp)) ≈ 10u"s"
    @test typical_Δt(tsp, consecutive_pairs(tsp)) ≈ 10u"s"
end

