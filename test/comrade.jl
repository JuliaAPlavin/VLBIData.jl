@testitem "comrade extract_table" begin
    import Comrade
    using Dates
    using Unitful
    using Uncertain

    antennas = [
        Antenna(:AA, xyz=(1e6, 2e6, 3e6)),
        Antenna(:BB, xyz=(4e6, 5e6, 6e6)),
    ]

    d0 = DateTime(2017, 4, 11, 0, 0, 0)
    uvtbl = [
        (
            datetime = d0,
            spec = VisSpec(Baseline((:AA, :BB)), UV(1e9, 2e9)),
            freq_spec = (freq=230e9u"Hz", width=2e9u"Hz"),
            value = (1.0 + 0.5im) ±ᵤ0.1,
        ),
        (
            datetime = d0 + Minute(1),
            spec = VisSpec(Baseline((:AA, :BB)), UV(1.1e9, 2.1e9)),
            freq_spec = (freq=230e9u"Hz", width=2e9u"Hz"),
            value = (0.9 + 0.4im) ±ᵤ0.15,
        ),
    ]

    obscoords = (ra=180.0u"°", dec=30.0u"°")
    obs = Comrade.extract_table(uvtbl; antennas, obscoords, objname="TestSrc")

    @test obs isa Comrade.EHTObservationTable
    @test length(obs) == 2
    config = obs.config
    @test config.source == :TestSrc
    @test length(config.tarr) == 2
end
