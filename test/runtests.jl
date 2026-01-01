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

    spec = VisSpec(Baseline(1, (2, 3)), UV(10, -20))
    @test Baseline(spec) == Baseline(1, (2, 3))
    @test UV(spec) == UV(10, -20)
    @test antennas(spec) == (Antenna(:ANT2, 2, SVector(NaN, NaN, NaN)), Antenna(:ANT3, 3, SVector(NaN, NaN, NaN)))
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
end

@testitem "frequencies" begin
    struct MyFreq
        freq::Float64
    end
    VLBI.frequency(x::MyFreq) = x.freq
    @test frequency((;freq_spec=MyFreq(10))) == 10
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(VLBIData; piracies=(;broken=true))

    import CompatHelperLocal as CHL
    CHL.@check()
end
