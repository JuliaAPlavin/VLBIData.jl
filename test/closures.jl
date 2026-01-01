@testitem "types" begin
    cas = ClosureAmpSpec((
        VisSpec(Baseline((2, 9)), UV([4.282665f6, -2.8318948f7])),
        VisSpec(Baseline((9, 3)), UV([-1.425433f6, 2.8318948f7])),
        VisSpec(Baseline((3, 5)), UV([1.425433f6, 2.8318948f7])),
        VisSpec(Baseline((5, 2)), UV([-4.282665f6, -2.8318948f7]))
    ))
    @test UVs(cas) == (UV(4.282665f6, -2.8318948f7), UV(-1.425433f6, 2.8318948f7), UV(1.425433f6, 2.8318948f7), UV(-4.282665f6, -2.8318948f7))
    @test antenna_names(cas) == (:ANT2, :ANT9, :ANT3, :ANT5)
    @test string(cas) == "ClosureAmp ANT2 - ANT9 - ANT3 - ANT5"
    @test UVarea(cas) ≈ 3.2329467f14

    cps = ClosurePhaseSpec((
        VisSpec(Baseline((2, 9)), UV([4.282665f6, -2.8318948f7])),
        VisSpec(Baseline((9, 3)), UV([-1.425433f6, 2.8318948f7])),
        VisSpec(Baseline((3, 2)), UV([1.425433f6, 2.8318948f7]))
    ))
    @test UVs(cps) == (UV(4.282665f6, -2.8318948f7), UV(-1.425433f6, 2.8318948f7), UV(1.425433f6, 2.8318948f7))
    @test antenna_names(cps) == (:ANT2, :ANT9, :ANT3)
    @test string(cps) == "ClosurePhase ANT2 - ANT9 - ANT3"
    @test UVarea(cps) ≈ 8.073353f13
end

@testitem "closure amplitude calculations" begin
    using Uncertain
    using DataManipulation
    using DataManipulation.FlexiGroups: GroupArray

    # nbl = bin(n, 2) ~ n^2/2
    # n_clamp_dof = nbl - n ~ n^2/2
    # n_clamp_max = 3*bin(n, 4) ~ n^4/8
    # redundancy_factor = n_clamp_max / n_clamp_dof ~ n^2/4
    # n_bl_in_clamp_max = (n-2)(n-3) ~ n^2
    # n_bl_in_clamp_max = 4*n_clamp_max / n_bl ~ n^4/2 / (n^2/2) ~ n^2
    # n_bl_clamp_min ~ 4*n_clamp_dof / n_bl ~ 4

    uvgr = GroupArray((;), [
        (spec=VisSpec(Baseline((1, 2)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((3, 4)), UV(0, 1)), value=1±ᵤ0.1,),
    ])
    clas = VLBI.closures_scan(ClosureAmpSpec, uvgr)
    nbls = @p clas flatmap(_.spec.vses, antenna_names(_2.bl)) map(extrema) groupmap(length)
    @test length(clas) == binomial(4, 4) * Int(factorial(3) / 2) == 3
    @test all(==(2), nbls)
    @test all(c -> c.value === 1±ᵤ0.2, clas)

    uvgr = GroupArray((;), [
        (spec=VisSpec(Baseline((1, 2)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 5)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 5)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((3, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((3, 5)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((4, 5)), UV(0, 1)), value=1±ᵤ0.1,),
    ])
    clas = VLBI.closures_scan(ClosureAmpSpec, uvgr)
    nbls = @p clas flatmap(_.spec.vses, antenna_names(_2.bl)) map(extrema) groupmap(length)
    @test length(clas) == binomial(5, 4) * Int(factorial(3) / 2) == 15
    @test all(==(6), nbls)
    @test all(c -> c.value === 1±ᵤ0.2, clas)
end

@testitem "closure phase calculations" begin
    using Uncertain
    using DataManipulation
    using DataManipulation.FlexiGroups: GroupArray

    uvgr = GroupArray((;), [
        (spec=VisSpec(Baseline((1, 2)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((3, 4)), UV(0, 1)), value=1±ᵤ0.1,),
    ])
    clps = VLBI.closures_scan(ClosurePhaseSpec, uvgr)
    nbls = @p clps flatmap(_.spec.vses, antenna_names(_2.bl)) map(extrema) groupmap(length)
    @test length(clps) == binomial(4, 3) * Int(factorial(2) / 2) == 4
    @test all(==(2), nbls)
    @test all(c -> c.value ≈ 1 ±ᵤ 0.17320508075688773, clps)


    uvgr = GroupArray((;), [
        (spec=VisSpec(Baseline((1, 2)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((1, 5)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 3)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((2, 5)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((3, 4)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((3, 5)), UV(0, 1)), value=1±ᵤ0.1,),
        (spec=VisSpec(Baseline((4, 5)), UV(0, 1)), value=1±ᵤ0.1,),
    ])
    clps = VLBI.closures_scan(ClosurePhaseSpec, uvgr)
    nbls = @p clps flatmap(_.spec.vses, antenna_names(_2.bl)) map(extrema) groupmap(length)
    @test length(clps) == binomial(5, 3) * Int(factorial(2) / 2) == 10
    @test all(==(3), nbls)
    @test all(c -> c.value ≈ 1 ±ᵤ 0.17320508075688773, clps)
end

@testitem "model evaluation" begin
    using InterferometricModels
    using StaticArrays
    using AccessorsExtra

    model = construct(CircularGaussian, flux=>1.0, fwhm_average=>0.1, coords=>SVector(0., 0.1))
    visspecs = (
        VisSpec(Baseline((2, 9)), UV([1, 2])),
        VisSpec(Baseline((9, 3)), UV([3, 4])),
        VisSpec(Baseline((3, 5)), UV([5, 6])),
        VisSpec(Baseline((5, 2)), UV([7, 8])),
        VisSpec(Baseline((3, 2)), UV([9, 10])),
    )
    cas = ClosureAmpSpec(visspecs[1:4])
    cps = ClosurePhaseSpec(visspecs[[1,2,5]])
    vises = (0.2586332168367465 + 0.7959911936251323im, -0.3322517918545724 + 0.24139505678833204im, -0.09224042204143246 - 0.06701658941426379im, 0.005534089140668198 - 0.01703217504149487im, 0.0015915127983935817 - 3.8980821087094936e-19im)
    @test visibility.(model, visspecs) == vises
    avis = vises[1] * vises[3] / (vises[2] * vises[4])
    pvis = vises[1] * vises[2] * vises[5]
    @test visibility(model, cas) ≈ avis
    @test visibility(visibility(model), cas) ≈ avis
    @test visibility(model, cps) ≈ pvis
end
