using TestItems
using TestItemRunner
@run_package_tests


@testitem "generic loading" begin
    cd(dirname(@__FILE__))
    @test VLBI.guess_type("./data/map.fits") == VLBI.FitsImage
    @test VLBI.guess_type("./data/vis.fits") == VLBI.UVData
    @test VLBI.guess_type("./data/vis_multichan.vis") == VLBI.UVData
    @test VLBI.guess_type("./data/difmap_model.mod") == MultiComponentModel
    @test VLBI.guess_type("./data/difmap_model_empty.mod") == MultiComponentModel
end

@testitem "img don't read data" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    using VLBIData: frequency
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/map.fits", read_data=false)
    @test_throws "not loaded" KeyedArray(img)

    @test VLBI.pixel_size(img) ≈ 0.2u"mas"  rtol=1e-5
    @test VLBI.pixel_steps(img) ≈ [-0.2u"mas", 0.2u"mas"]  rtol=1e-5
    @test VLBI.pixel_area(img) ≈ 0.04u"mas^2"  rtol=1e-5
    @test frequency(img) ≈ 4.344u"GHz"

    bm = beam(img)
    @test VLBI.load(Beam, "./data/map.fits") === bm
    @test effective_area(bm) ≈ 6.5519451u"mas^2"  rtol=1e-5
    @test fwhm_max(bm) ≈ 4.02887356u"mas"  rtol=1e-5
    @test fwhm_min(bm) ≈ 1.43523228u"mas"  rtol=1e-5
    @test position_angle(bm) ≈ -0.046499863  rtol=1e-5
    @test intensity(ModelComponent(bm))(SVector(0.5u"mas", 0u"mas")) ≈ 0.714721  rtol=1e-5
end
    
@testitem "img read data" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    using VLBIData: frequency
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/map.fits", read_data=true)
    @test img.data == VLBI.load(VLBI.FitsImage, "./data/map.fits").data
    KA = KeyedArray(img)
    @test KA == VLBI.load(KeyedArray, "./data/map.fits")
    @test size(KA) == (512, 512)
    @test mean(KA) ≈ 2.4069064f-5
    @test maximum(KA) ≈ 0.021357307f0
    @test KA[123, 456] ≈ -7.372097f-5
    @test axiskeys(KA, :ra)  .|> ustrip ≈ 51:-0.2:-51.2  atol=1e-3
    @test axiskeys(KA, :dec) .|> ustrip ≈ -51.2:0.2:51   atol=1e-3
    @test axiskeys(KA, :ra) isa AbstractRange
    @test axiskeys(KA, :dec) isa AbstractRange
end

@testitem "img save" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using AxisKeys

    tmpf = tempname()
    
    img = KeyedArray(rand(512, 256), ra=range(-10, 10, length=512)u"mas", dec=range(-20, 20, length=256)u"mas")
    @test_throws Exception VLBI.save(tmpf, img)

    img = KeyedArray(rand(512, 256), ra=range(-10, 10, length=512)u"mas", dec=range(-20, 20, length=256)u"mas")*u"Jy"
    VLBI.save(tmpf, img; freq=123)

    img_r = VLBI.load(tmpf).data
    @test axiskeys(img_r, 1) ≈ axiskeys(img, 1)
    @test axiskeys(img_r, 2) ≈ axiskeys(img, 2)
    @test AxisKeys.keyless_unname(img_r) ≈ AxisKeys.keyless_unname(ustrip.(u"Jy", img))
    @test_broken VLBI.frequency(VLBI.load(tmpf)) == 123

    img = KeyedArray(rand(512, 256), ra=range(-5, 10, length=512)u"mas", dec=range(-7, 1, length=256)u"mas")*u"Jy"
    VLBI.save(tmpf, img)

    img_r = VLBI.load(tmpf).data
    @test axiskeys(img_r, 1) ≈ axiskeys(img, 1)
    @test axiskeys(img_r, 2) ≈ axiskeys(img, 2)
    @test AxisKeys.keyless_unname(img_r) ≈ AxisKeys.keyless_unname(ustrip.(u"Jy", img))
end

@testitem "img read clean" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    using VLBIData: frequency
    cd(dirname(@__FILE__))

    mod = VLBI.load(MultiComponentModel, "./data/map.fits")
    @test length(components(mod)) == 361
    @test sum(flux, components(mod)) ≈ 0.038659f0u"Jy"
    @test mean(first ∘ coords, components(mod)) ≈ -0.913573508654587u"mas"
    @test mean(last ∘ coords, components(mod)) ≈ 8.145706523588233u"mas"
end

@testitem "img stacked" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    using VLBIData: frequency
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/map_stacked.fits", read_data=true)
    KA = KeyedArray(img)
    @test size(KA) == (512, 512)
    @test maximum(KA) ≈ 0.5073718945185344
    @test KA[123, 456] ≈ -1.0848011697817128e-5
    @test axiskeys(KA, :ra)  .|> ustrip ≈ 25.5:-0.1:-25.6  atol=1e-3
    bm = beam(img)
    @test fwhm_min(bm) == fwhm_max(bm)

    @test_broken mod = VLBI.load(MultiComponentModel, "./data/map_stacked.fits")
end

@testitem "img nonstandard header names" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using StaticArrays
    using Statistics
    using AxisKeys
    using VLBIData: frequency
    cd(dirname(@__FILE__))

    img = VLBI.load("./data/sampling_mean.fits", read_data=true)
    KA = KeyedArray(img)
    @test size(KA) == (256, 256)
    @test maximum(KA) ≈ 0.0022970616631588308
    @test KA[123, 222] ≈ 1.0227093348779382e-6
    @test axiskeys(KA, :ra)  .|> ustrip ≈ -0.3984375:0.003125:0.3984375  atol=1e-3
    @test_throws "key \"BMAJ\" not found" beam(img)
end

@testitem "uvf simple" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIData: frequency
    using VLBIData.Uncertain
    using Statistics
    using StaticArrays
    using DataManipulation
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/vis.fits")
    @test uv.header.object == "J1033+6051"
    @test Date(uv.header) === Date(2010, 12, 24)
    @test uv.header.stokes == [:RR, :LL, :RL, :LR]
    @test frequency(uv.header) ≈ 15.33u"GHz" rtol=1e-2
    @test length(uv.freq_windows) == 8
    @test length(uv.ant_arrays) == 1
    antarr = only(uv.ant_arrays)
    @test antarr.name == "VLBA"
    @test map(a -> a.id, antarr.antennas) == 1:10
    @test map(a -> a.name, antarr.antennas) == [:BR, :FD, :HN, :KP, :LA, :MK, :NL, :OV, :PT, :SC]

    # old table extraction:
    df = table(uv)
    @test Tables.rowaccess(df)
    @test length(df) == 160560
    @test all(∈(Tables.schema(df).names), [:uv, :visibility, :if_ix, :datetime])
    @test all(isconcretetype, Tables.schema(df).types)
    @test df[1].uv === VLBI.UV(-6.4672084f7, 8.967541f7)
    target = (
        baseline = VLBIData.Baseline(1, (8, 9), (:OV, :PT)),
        datetime = DateTime("2010-12-24T04:57:34.999"),
        stokes = :RR,
        if_ix = Int8(1),
        if_spec = VLBIData.FrequencyWindow(1.5329522f10u"Hz", 8.0f6u"Hz", 1, 1),
        uv_m = VLBIData.UV(417953.38f0u"m", 860236.9f0u"m"),
        w_m = -179245.12f0u"m",
        uv = VLBIData.UV(2.1377114f7, 4.3998644f7),
        w = -9.167873f6,
        visibility = 0.4755191f0 + 0.0043343306f0im,
        weight = 946.5216f0
    )
    res = filter(r -> r.baseline == target.baseline && r.datetime == target.datetime && r.if_ix == target.if_ix && r.stokes == target.stokes, df) |> only
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    df_cols = Tables.columntable(df)
    @test mean(df_cols.uv_m) ≈ SVector(-174221.2u"m", 314413.6u"m")
    @test mean(v->abs.(v), df_cols.uv) ≈ VLBI.UV(1.0778852f8, 7.967933f7)
    @test mean(df_cols.visibility) ≈ 0.2495917 + 0.0010398296im
    @test first(df) == first(table(uv, pyimport))
    @test df == table(uv, pyimport)

    # new table extraction:
    tbl = uvtable(uv; stokes=UniversalSet)
    @test length(tbl) == 160560
    tbl = uvtable(uv)
    @test length(tbl) == 80280
    @test isconcretetype(eltype(tbl))
    @test tbl[12345] == (
        datetime = DateTime("2010-12-24T08:06:25"),
        stokes = :RR,
        if_ix = 6,
        if_spec = VLBI.FrequencyWindow(1.5369459f10u"Hz", 8.0f6u"Hz", 1, 1),
        spec = VLBI.VisSpec(VLBI.Baseline(1, (2, 9), (:FD, :PT)), [4.282665f6, -2.8318948f7]),
        value = (0.48758677f0 - 0.09014242f0im) ±ᵤ 0.040410895f0
    )
    @test VLBI.antennas(tbl[12345]) == (VLBI.Antenna(:FD, 2, VLBI.SVector(NaN, NaN, NaN)), VLBI.Antenna(:PT, 9, VLBI.SVector(NaN, NaN, NaN)))
    @test VLBI.Baseline(tbl[12345]) == VLBI.Baseline(1, (2, 9), (:FD, :PT))
    @test VLBI.UV(tbl[12345]) == VLBI.UV(4.282665f6, -2.8318948f7)
    @test VLBI.visibility(tbl[12345]) == (0.48758677f0 - 0.09014242f0im) ±ᵤ 0.040410895f0
end

@testitem "uvf multichannel" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIData: frequency
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/vis_multichan.vis")
    @test length(uv.freq_windows) == 8
    df = table(uv)
    target = (
        baseline = VLBIData.Baseline(1, (3, 5), (:HN, :LA)),
        datetime = DateTime("1996-06-05T19:16:45.001"),
        stokes = :LL,
        if_ix = Int8(3),
        if_spec = VLBIData.FrequencyWindow(4.84099f9u"Hz", 8.0f6u"Hz", 16, 1),
        uv_m = VLBIData.UV(2.2424022f6u"m", 794705.06f0u"m"),
        w_m = -1.838948f6u"m",
        uv = VLBIData.UV(3.6239796f7, 1.2843347f7),
        w = -2.9719512f7,
        visibility = -0.21484283f0 - 0.35979474f0im,
        weight = 3.0233376f0
    )
    res = filter(r -> r.baseline.ant_ids == target.baseline.ant_ids && r.datetime == target.datetime && r.if_ix == target.if_ix && r.stokes == target.stokes, df)[1]
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    @test df == table(uv, pyimport)
end

@testitem "uvf EHT 1" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIData: frequency
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/hops_3600_OJ287_LO+HI.medcal_dcal_full.uvfits")
    @test length(uv.freq_windows) == 2
    df = table(uv)
    target = (baseline = VLBIData.Baseline(1, (2, 3), (:AP, :AZ)), datetime = Dates.DateTime("2017-04-10T00:58:45.005"), stokes = :RR, if_ix = Int8(1), if_spec = VLBIData.FrequencyWindow(2.270707f11u"Hz", 1.856f9u"Hz", 1, 1), uv_m = VLBI.UV(3.9143352f6u"m", -5.933179f6u"m"), w_m = 0.0f0u"m", uv = VLBI.UV(2.9769375f9, -4.5123123f9), w = 0.0f0, visibility = 0.5662228f0 + 0.014944182f0im, weight = 117.818825f0)
    res = filter(r -> r.baseline.ant_ids == target.baseline.ant_ids && r.datetime == target.datetime && r.if_ix == target.if_ix && r.stokes == target.stokes, df)[1]
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    @test df == table(uv, pyimport)
end

@testitem "uvf EHT 2" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIData: frequency
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/SR1_3C279_2017_101_hi_hops_netcal_StokesI.uvfits")
    @test length(uv.freq_windows) == 1
    df = table(uv)
    target = (baseline = VLBIData.Baseline(1, (1, 2), (:AA, :AP)), datetime = Dates.DateTime("2017-04-11T02:14:55"), stokes = :RR, if_ix = Int8(1), if_spec = VLBIData.FrequencyWindow(2.290707f11u"Hz", 1.856f9u"Hz", 1, 1), uv_m = VLBI.UV(895.6223f0u"m", -2436.5188f0u"m"), w_m = 0.0f0u"m", uv = VLBI.UV(687115.25f0, -1.8692805f6), w = 0.0f0, visibility = -1.1796279f0 - 7.8919725f0im, weight = 40441.19f0)
    res = filter(r -> r.baseline.ant_ids == target.baseline.ant_ids && r.datetime == target.datetime && r.if_ix == target.if_ix && r.stokes == target.stokes, df)[1]
    @test res == target
    @test map(typeof, res) == map(typeof, target)
    @test df == table(uv, pyimport)
end

@testitem "uvf EHT 3" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    using Dates
    using PyCall
    using VLBIData: frequency
    using Statistics
    using StaticArrays
    using Tables
    cd(dirname(@__FILE__))

    uv = VLBI.load(VLBI.UVData, "./data/datafile_01-01_230GHz.uvfits")
    @test length(uv.freq_windows) == 1
    df = table(uv)
    @test df[1].visibility ≈ -0.0061733737f0 + 0.052109245f0im
    @test df == table(uv, pyimport)
end

@testitem "difmap model" begin
    using Unitful, UnitfulAstro, UnitfulAngles
    cd(dirname(@__FILE__))

    mod = VLBI.load("./data/difmap_model.mod")
    @test length(components(mod)) == 4
    @test map(flux, components(mod)) |> collect ≈ [0.522, 0.0217, 0.0176, 0.2145]u"Jy"  rtol=0.01
    @test map(fwhm_max, components(mod)) |> collect ≈ [0.135, 0.463, 1.99, 0]u"mas"  rtol=0.01
    @test coords(components(mod)[1]) ≈ [-0.000662, -0.00123]u"mas" rtol=0.01

    mod_empty = VLBI.load("./data/difmap_model_empty.mod")
    @test isempty(components(mod_empty))

    mod_clean = VLBI.load("./data/difmap_model_clean.mod")
    @test length(components(mod_clean)) == 631
    @test isconcretetype(eltype(components(mod_clean)))
    @test eltype(components(mod_clean)) <: Point

    mod_map = VLBI.load(MultiComponentModel, "./data/map.fits")

    tmpf = tempname()
    VLBI.save(tmpf, mod)
    @test VLBI.load(tmpf) == mod
    VLBI.save(tmpf, mod_empty)
    @test VLBI.load(tmpf) == mod_empty
    VLBI.save(tmpf, mod_clean)
    @test VLBI.load(tmpf) == mod_clean
    VLBI.save(tmpf, mod_map)
    @test VLBI.load(tmpf) ≈ mod_map  rtol=1e-4  # approx because saving to *.mod involves rounding; also Float32 vs 64
end

# @testitem "RFC" begin
#     using AstroRFC: RFC
#     using ProgressMeter
#     using PyCall

#     @testset "vis" begin
#         rfci = RFC.Files()
#         @showprogress for f in RFC.files(rfci, suffix="vis")
#             try
#                 uv = VLBI.load(f)
#                 df = table(uv)
#                 @test df == table(uv, pyimport)
#             catch e
#                 @show abspath(f) e
#                 rethrow()
#             end
#         end
#     end

#     @testset "maps" begin
#         rfci = RFC.Files()
#         @showprogress for f in RFC.files(rfci, suffix="map", extension="fits")
#             try
#                 VLBI.load(f)
#             catch e
#                 @show f e
#                 rethrow()
#             end
#             try
#                 VLBI.load(MultiComponentModel, f)
#             catch e
#                 @show f e
#                 rethrow()
#             end
#             try
#                 VLBI.load(Beam, abspath(f))
#             catch e
#                 if e isa KeyError && e.key == "BMAJ"
#                     @warn "Missing BMAJ" f
#                 else
#                     @show f e
#                     rethrow()
#                 end
#             end
#         end
#     end
# end

@testitem "Baseline type" begin
    using Accessors

    bl = Baseline(1, (2, 3), (:A, :B))
    @test bl.ant_ids == (2, 3)
    @test bl.ant_names == (:A, :B)
    @test VLBI.antennas(bl) == (VLBI.Antenna(:A, 2, VLBI.SVector(NaN, NaN, NaN)), VLBI.Antenna(:B, 3, VLBI.SVector(NaN, NaN, NaN)))

    Accessors.test_getset_laws(VLBI.antennas, bl, (VLBI.Antenna(:XX, 4, VLBI.SVector(NaN, NaN, NaN)), VLBI.Antenna(:YY, 5, VLBI.SVector(NaN, NaN, NaN))), (VLBI.Antenna(:ZZ, 6, VLBI.SVector(NaN, NaN, NaN)), VLBI.Antenna(:WW, 7, VLBI.SVector(NaN, NaN, NaN))))

    # backwards compatibility
    bl = Baseline(1, (2, 3))
    @test bl == Baseline(1, (2, 3))
    @test bl != Baseline(1, (2, 4))
    @test bl.array_ix == 1
    @test bl.ants_ix == (2, 3)
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(VLBIData; ambiguities=false, piracies=false)
    Aqua.test_ambiguities(VLBIData)

    import CompatHelperLocal as CHL
    CHL.@check()
end
