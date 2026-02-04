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
end

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

    # round-trip: EHTObservationTable → uvtable
    result = VLBI.uvtable(obs)
    @test length(result) == 2
    @test issetequal(keys(result[1]), (:datetime, :spec, :freq_spec, :value))

    @test result[1].spec == VisSpec(Baseline((:AA, :BB)), UV(1e9, 2e9))
    @test result[2].spec == VisSpec(Baseline((:AA, :BB)), UV(1.1e9, 2.1e9))

    @test U.value(result[1].value) == 1.0 + 0.5im
    @test U.uncertainty(result[1].value) == 0.1
    @test U.value(result[2].value) == 0.9 + 0.4im
    @test U.uncertainty(result[2].value) == 0.15

    @test result[1].freq_spec ≈ 230e9u"Hz"

    @test result[1].datetime == d0
    @test result[2].datetime == d0 + Minute(1)
end
