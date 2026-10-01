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
