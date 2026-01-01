@testitem "visibilties" begin
    using Unitful
    using Dates
    using Uncertain

    avgd = VLBI.average_bytime([
        (datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 1)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 0, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 2)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 0, 20), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 3)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 20, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 4)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5)), freq_spec=123, stokes=:LL, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5)), freq_spec=124, stokes=:LL, value=1±ᵤ0.1),
    ], 10u"minute")

    @test avgd == [
        (freq_spec = 123, stokes = :RR, count = 2, value = 1±ᵤ0.07071067811865475, datetime = DateTime("2020-01-01T00:05:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 1.5))),
        (freq_spec = 123, stokes = :RR, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:05:00"), spec = VisSpec(Baseline(1, (1, 3)), UV(0, 3.))),
        (freq_spec = 123, stokes = :RR, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 4.))),
        (freq_spec = 123, stokes = :LL, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 5.))),
        (freq_spec = 124, stokes = :LL, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 5.))),
    ]
end
