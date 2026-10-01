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
