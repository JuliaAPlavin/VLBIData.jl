using TestItems
using TestItemRunner
@run_package_tests

@testitem "Antenna" begin
    using AccessorsExtra

    a = Antenna(:MyAnt, 2, (1.0, 2.0, 3.0))
    b = Antenna(:XXAnt, 4, (2, 3, 5))
    @test antennas(a) == (a,)
    @test antennas((a, b)) == (a, b)
    @test antennas([a, b]) == [a, b]
    @test (@oget UV(a) UV((a, b)) 123) == 123
    
    @test Antenna(:A, 1) < Antenna(:B, 2)
    @test Antenna(:B, 1) < Antenna(:A, 2)
    @test !(Antenna(:B, 2) < Antenna(:A, 1))
    @test Antenna(:A, 1) ≥ Antenna(:A, 1)
    @test Antenna(:A, 1) ≥ Antenna(:A, 0)

    @test string(a) == "Antenna #2: MyAnt"
end

@testitem "Baseline" begin
    using AccessorsExtra
    using StaticArrays

    bl = Baseline(1, (2, 3))
    @test bl.ant_ids == (2, 3)
    @test bl.ant_names == (:ANT2, :ANT3)
    @test bl == Baseline(1, (2, 3))
    @test bl != Baseline(1, (2, 4))
    @test bl != Baseline(1, (2, 3), (:ANT2, :ANT4))

    bl = Baseline(1, (2, 3), (:A, :B))
    @test bl.ant_ids == (2, 3)
    @test bl.ant_names == (:A, :B)
    @test bl == Baseline(1, (2, 3), (:A, :B))
    @test bl != Baseline(1, (2, 4), (:A, :B))
    @test bl != Baseline(1, (2, 3), (:A, :C))

    @test antennas(bl) == (Antenna(:A, 2, SVector(NaN, NaN, NaN)), Antenna(:B, 3, SVector(NaN, NaN, NaN)))
    Accessors.test_getset_laws(antennas, bl, (Antenna(:XX, 4, SVector(NaN, NaN, NaN)), Antenna(:YY, 5, SVector(NaN, NaN, NaN))), (Antenna(:ZZ, 6, SVector(NaN, NaN, NaN)), Antenna(:WW, 7, SVector(NaN, NaN, NaN))))

    @test (@oget UV(bl) 123) == 123

    @test Baseline((;spec=bl)) == bl
    @test antennas((;spec=bl)) == antennas(bl)

    @test string(bl) == "Baseline 2-3: A - B"
end

@testitem "UV" begin
    using AccessorsExtra

    @test (@oget UV([1,2]) 123) == UV(1, 2)

    uv = UV(1, 2)
    @test (@oget UV(uv) 123) == uv
    @test (@oget antennas(uv) 123) == 123
    @test (@oget Baseline(uv) 123) == 123

    @test UV((;spec=uv)) == uv
end

@testitem "visspec + visampspec" begin
    using AccessorsExtra
    using StaticArrays

    spec = VLBI.VisSpec0(UV(10, -20))
    @test UV(spec) == UV(10, -20)
    @test !hasoptic(spec, Baseline)
    @test hasoptic(spec, UV)
    @test !hasoptic(spec, antennas)
    @test_broken conj(spec) == VLBI.VisSpec0(UV(-10, 20))
    @test string(spec) == "VisSpec0: [10, -20]"
    @test spec == VLBI.VisSpec0(UV(10, -20.))
    @test spec != VLBI.VisSpec0(UV(10, -21))

    spec = VisSpec(Baseline(1, (2, 3)), UV(10, -20))
    @test Baseline(spec) == Baseline(1, (2, 3))
    @test UV(spec) == UV(10, -20)
    @test antennas(spec) == (Antenna(:ANT2, 2, SVector(NaN, NaN, NaN)), Antenna(:ANT3, 3, SVector(NaN, NaN, NaN)))
    @test hasoptic(spec, Baseline)
    @test hasoptic(spec, UV)
    @test hasoptic(spec, antennas)
    @test conj(spec) == VisSpec(Baseline(1, (3, 2)), UV(-10, 20))
    @test string(spec) == "Vis 2-3: ANT2 - ANT3"
    @test spec == VisSpec(Baseline(1, (2, 3)), UV(10, -20.))
    @test spec != VisSpec(Baseline(1, (2, 4)), UV(10, -20))
    @test spec != VisSpec(Baseline(1, (2, 3)), UV(10, -21))

    aspec = VisAmpSpec(spec)
    @test aspec == VisAmpSpec(Baseline(1, (2, 3)), UV(10, -20))
    @test Baseline(aspec) == Baseline(1, (2, 3))
    @test UV(aspec) == UV(10, -20)
    @test antennas(aspec) == (Antenna(:ANT2, 2, SVector(NaN, NaN, NaN)), Antenna(:ANT3, 3, SVector(NaN, NaN, NaN)))
    @test VisSpec(aspec) == spec
    @test conj(aspec) == VisAmpSpec(VisSpec(Baseline(1, (3, 2)), UV(-10, 20)))
    @test string(aspec) == "VisAmp 2-3: ANT2 - ANT3"
    @test aspec == VisAmpSpec(Baseline(1, (2, 3)), UV(10, -20.))
    @test aspec != VisAmpSpec(Baseline(1, (2, 4)), UV(10, -20))
    @test aspec != VisAmpSpec(Baseline(1, (2, 3)), UV(10, -21))

    @test Baseline((;spec)) == Baseline(spec)
    @test antennas((;spec)) == antennas(spec)
    @test UV((;spec)) == UV(spec)
    @test VLBI.conjvis((;spec, value=1+2im)) == (spec=conj(spec), value=1-2im)
    @test VLBI.conjvis((;spec, value=1+2im, stokes=:RR)) == (spec=conj(spec), value=1-2im, stokes=:RR)
    @test VLBI.conjvis((;spec, value=1+2im, stokes=:RL)) == (spec=conj(spec), value=1-2im, stokes=:LR)
end

@testitem "frequencies" begin
    using Unitful
    struct MyFreq
        freq::Float64
    end
    VLBI.frequency(x::MyFreq) = x.freq
    @test frequency((;freq_spec=MyFreq(10))) == 10
    @test frequency(10u"GHz") == 10u"GHz"
    @test frequency((;freq_spec=10u"GHz")) == 10u"GHz"
end

@testitem "renumbering" begin
    uvtbl_orig = [
        (spec=VLBI.VisSpec0(UV(0, 1)), freq_spec=123, stokes=:RR, value=1),
        (spec=VLBI.VisSpec0(UV(0, 2)), freq_spec=123, stokes=:RR, value=2),
        (spec=VLBI.VisSpec0(UV(0, 5)), freq_spec=124, stokes=:LL, value=3),
    ]
    @test VLBI.uv_reindex(uvtbl_orig) == [
        (spec=VLBI.VisSpec0(UV(0, 1)), freq_spec=123, stokes=:RR, value=1),
        (spec=VLBI.VisSpec0(UV(0, 2)), freq_spec=123, stokes=:RR, value=2),
        (spec=VLBI.VisSpec0(UV(0, 5)), freq_spec=124, stokes=:LL, value=3),
    ]

    uvtbl_orig = [
        (spec=VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 1)), stokes=:RR, value=1),
        (spec=VisSpec(Baseline(1, (2, 2), (:A, :B)), UV(0, 2)), stokes=:RR, value=2),
        (spec=VisSpec(Baseline(1, (1, 2), (:C, :B)), UV(0, 5)), stokes=:LL, value=3),
    ]
    @test VLBI.uv_reindex(uvtbl_orig) == [
        (spec=VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 1)), stokes=:RR, value=1),
        (spec=VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 2)), stokes=:RR, value=2),
        (spec=VisSpec(Baseline(1, (3, 2), (:C, :B)), UV(0, 5)), stokes=:LL, value=3),
    ]

    uvtbl_orig = [
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 1))), freq_spec=123, stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (2, 2), (:A, :B)), UV(0, 2))), freq_spec=123, stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:C, :B)), UV(0, 5))), freq_spec=124, stokes=:LL, value=3),
    ]
    @test VLBI.uv_reindex(uvtbl_orig) == [
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 1))), freq_spec=123, stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 2))), freq_spec=123, stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (3, 2), (:C, :B)), UV(0, 5))), freq_spec=124, stokes=:LL, value=3),
    ]

    uvtbl_orig = [
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 1))), freq_spec=(;ix=123), stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (2, 2), (:A, :B)), UV(0, 2))), freq_spec=(;ix=123), stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:C, :B)), UV(0, 5))), freq_spec=(;ix=124), stokes=:LL, value=3),
    ]
    @test VLBI.uv_reindex(uvtbl_orig) == [
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 1))), freq_spec=(;ix=1), stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (1, 2), (:A, :B)), UV(0, 2))), freq_spec=(;ix=1), stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline(1, (3, 2), (:C, :B)), UV(0, 5))), freq_spec=(;ix=2), stokes=:LL, value=3),
    ]
end

@testitem "visibility error rescaling" begin
    using VLBIFiles
    using VLBIFiles: VLBI
    using Unitful
    
    uvtbl = VLBI.load(joinpath(pkgdir(VLBIFiles), "test/data/SR1_3C279_2017_101_hi_hops_netcal_StokesI.uvfits")) |> uvtable
    @test VLBI.find_errmul(uvtbl, VLBI.ConsecutiveDifferencesStandard()) ≈ 0.6789 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.ConsecutiveDifferencesStandard(2u"minute")) ≈ 0.6789 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.ConsecutiveDifferencesStandard(2u"minute", rayleigh_q=0.2)) ≈ 0.605 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.CoherentAverageScatter()) ≈ 0.7293 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.CoherentAverageScatter(2u"minute")) ≈ 0.7260 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.CoherentAverageScatter(2u"minute", min_cnt_avg=5)) ≈ 0.7275 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.CoherentAverageScatter(0.7u"minute", min_cnt_avg=5)) ≈ 0.6855 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.ErrMulSame(VLBI.ConsecutiveDifferencesStandard(), VLBI.ConsecutiveDifferencesStandard(2u"minute"), rtol=0.2)) ≈ 0.6789 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.ErrMulSame(VLBI.ConsecutiveDifferencesStandard(), VLBI.CoherentAverageScatter(), rtol=0.2)) ≈ 0.7042 rtol=1e-3
    @test_throws "don't agree" VLBI.find_errmul(uvtbl, VLBI.ErrMulSame(VLBI.ConsecutiveDifferencesStandard(), VLBI.CoherentAverageScatter(), rtol=0.05))
    @test length(VLBI.rescale_visibility_errors(uvtbl, VLBI.ConsecutiveDifferencesStandard())) == length(uvtbl)

    uvtbl = VLBI.load(joinpath(pkgdir(VLBIFiles), "test/data/vis_multichan.vis")) |> uvtable
    @test VLBI.find_errmul(uvtbl, VLBI.ConsecutiveDifferencesStandard()) ≈ 0.4298 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.CoherentAverageScatter()) ≈ 0.4306 rtol=1e-3

    uvtbl = VLBI.load(joinpath(pkgdir(VLBI), "test/data/datafile_01-01_230GHz.uvfits")) |> uvtable
    @test VLBI.find_errmul(uvtbl, VLBI.ConsecutiveDifferencesStandard()) ≈ 27.79 rtol=1e-3
    @test VLBI.find_errmul(uvtbl, VLBI.ConsecutiveDifferencesStandard(2u"minute")) |> isnothing
    @test VLBI.find_errmul(uvtbl, VLBI.CoherentAverageScatter()) ≈ 26.30 rtol=1e-3
    @test_throws "couldn't be estimated" VLBI.find_errmul(uvtbl, VLBI.ErrMulSame(VLBI.ConsecutiveDifferencesStandard(), VLBI.ConsecutiveDifferencesStandard(2u"minute"), rtol=0.2))
end

@testitem "comradebase" begin
    import ComradeBase as CB
    using Unitful

    p = (; U = 0.1 * randn(60), V = 0.1 * randn(60), Ti = collect(Float64, 1:60), Fr = fill(230.0e9, 60))
    d = rand(60)
    m = CB.UnstructuredMap(d, CB.UnstructuredDomain(p))
    uvtbl = VLBI.uvtable(m)
    @test length(uvtbl) == 60
    @test issetequal(keys(uvtbl[10]), (:spec, :freq_spec, :value, :datetime))
    @test uvtbl[10].spec::UV == UV(p.U[10], p.V[10])
    @test uvtbl[10].freq_spec ≈ 230u"GHz"
    @test uvtbl[10].value == d[10]
    @test uvtbl[10].datetime == p.Ti[10]

    m = CB.UnstructuredMap(d, CB.UnstructuredDomain(p[(:U, :V)]))
    uvtbl = VLBI.uvtable(m)
    @test length(uvtbl) == 60
    @test issetequal(keys(uvtbl[10]), (:spec, :value))
    @test uvtbl[10].spec::UV == UV(p.U[10], p.V[10])
    @test uvtbl[10].value == d[10]

    m = CB.UnstructuredMap(d, CB.UnstructuredDomain(p[(:U, :V, :Fr)]))
    uvtbl = VLBI.uvtable(m)
    @test length(uvtbl) == 60
    @test issetequal(keys(uvtbl[10]), (:spec, :value, :freq_spec))
    @test uvtbl[10].spec::UV == UV(p.U[10], p.V[10])
    @test uvtbl[10].value == d[10]
    @test uvtbl[10].freq_spec ≈ 230u"GHz"
end

@testitem "polarizedtypes extension (manual table)" begin
    using PolarizedTypes
    using VLBIData

    # Create a minimal fake uvtable with all 4 stokes
    uvtbl = [
        (datetime=1, freq_spec=100.0, spec=1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=1, stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=1, stokes=:RL, value=13.0),
        (datetime=1, freq_spec=100.0, spec=1, stokes=:LL, value=14.0),
        (datetime=2, freq_spec=200.0, spec=2, stokes=:RR, value=21.0),
        (datetime=2, freq_spec=200.0, spec=2, stokes=:LR, value=22.0),
        (datetime=2, freq_spec=200.0, spec=2, stokes=:RL, value=23.0),
        (datetime=2, freq_spec=200.0, spec=2, stokes=:LL, value=24.0),
    ]
    cm_tbl = VLBI.uvtable_values_to(CoherencyMatrix, uvtbl)

    @test length(cm_tbl) == 2
    @test cm_tbl[1].value isa CoherencyMatrix
    @test cm_tbl[2].value isa CoherencyMatrix
    @test cm_tbl[1].value == [11 13; 12 14]
    @test cm_tbl[2].value == [21 23; 22 24]
    @test !haskey(cm_tbl[1], :stokes)

    # Error case: group with not 4 elements (missing :LL)
    uvtbl_err = [
        (datetime=1, freq_spec=100.0, spec=1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=1, stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=1, stokes=:RL, value=13.0),
    ]
    @test_throws "expected 4" VLBI.uvtable_values_to(CoherencyMatrix, uvtbl_err)
end

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
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(VLBIData; piracies=(;broken=true))

    import CompatHelperLocal as CHL
    CHL.@check()
end
