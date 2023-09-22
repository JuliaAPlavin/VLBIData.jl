using Test
using Statistics
using AxisKeys
using Unitful, UnitfulAstro, UnitfulAngles
using StaticArrays
using Dates
using Tables
using AstroRFC: RFC
using VLBIData
using VLBIData: frequency
using ProgressMeter
using PyCall


# @testset "generic loading" begin
#     @test VLBI.guess_type("./data/map.fits") == VLBI.FitsImage
#     @test VLBI.guess_type("./data/vis.fits") == VLBI.UVData
#     @test VLBI.guess_type("./data/vis_multichan.vis") == VLBI.UVData
#     @test VLBI.guess_type("./data/difmap_model.mod") == MultiComponentModel
#     @test VLBI.guess_type("./data/difmap_model_empty.mod") == MultiComponentModel
# end


# @testset "fits image" begin
#     @testset "read nothing" begin
#         img = VLBI.load("./data/map.fits", read_data=false)
#         @test img.data === nothing

#         @test VLBI.pixel_size(img) ≈ 0.2u"mas"  rtol=1e-5
#         @test VLBI.pixel_steps(img) ≈ [-0.2u"mas", 0.2u"mas"]  rtol=1e-5
#         @test VLBI.pixel_area(img) ≈ 0.04u"mas^2"  rtol=1e-5
#         @test frequency(img) ≈ 4.344u"GHz"

#         bm = beam(img)
#         @test beam("./data/map.fits") === bm
#         @test effective_area(bm) ≈ 6.5519451u"mas^2"  rtol=1e-5
#         @test fwhm_max(bm) ≈ 4.02887356u"mas"  rtol=1e-5
#         @test fwhm_min(bm) ≈ 1.43523228u"mas"  rtol=1e-5
#         @test position_angle(bm) ≈ -0.046499863  rtol=1e-5
#         @test intensity(bm)(SVector(0.5u"mas", 0u"mas")) ≈ 0.714721  rtol=1e-5
#     end
    
#     @testset "read data" begin
#         img = VLBI.load("./data/map.fits", read_data=true)
#         @test size(img.data) == (512, 512)
#         @test mean(img.data) ≈ 2.4069064f-5
#         @test maximum(img.data) ≈ 0.021357307f0
#         @test img.data[123, 456] ≈ -7.372097f-5
#         @test axiskeys(img.data, :ra)  .|> ustrip ≈ 51:-0.2:-51.2  atol=1e-3
#         @test axiskeys(img.data, :dec) .|> ustrip ≈ -51.2:0.2:51   atol=1e-3
#         @test axiskeys(img.data, :ra) isa AbstractRange
#         @test axiskeys(img.data, :dec) isa AbstractRange
#     end
    
#     @testset "read clean" begin
#         mod = VLBI.load(MultiComponentModel, "./data/map.fits")
#         @test length(components(mod)) == 361
#         @test sum(flux, components(mod)) ≈ 0.038659f0u"Jy"
#         @test mean(first ∘ coords, components(mod)) ≈ -0.913573508654587u"mas"
#         @test mean(last ∘ coords, components(mod)) ≈ 8.145706523588233u"mas"
#     end

#     @testset "stacked image" begin
#         img = VLBI.load("./data/map_stacked.fits", read_data=true)
#         @test size(img.data) == (512, 512)
#         @test maximum(img.data) ≈ 0.5073718945185344
#         @test img.data[123, 456] ≈ -1.0848011697817128e-5
#         @test axiskeys(img.data, :ra)  .|> ustrip ≈ 25.5:-0.1:-25.6  atol=1e-3
#         bm = beam(img)
#         @test fwhm_min(bm) == fwhm_max(bm)

#         @test_broken mod = VLBI.load(MultiComponentModel, "./data/map_stacked.fits")
#     end
# end

@testset "uvfits" begin
    @testset "simple" begin
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

        df = table(uv, pyimport)
        @test Tables.rowaccess(df)
        @test length(df) == 160560
        @test all(∈(Tables.schema(df).names), [:uv, :visibility, :if_ix, :datetime])
        @test all(isconcretetype, Tables.schema(df).types)
        @test df[1].uv === VLBI.UV(-6.4672084f7, 8.967541f7)
        target = (
            baseline = VLBIData.Baseline(1, (8, 9)),
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
        @test df[123] == target
        @test map(typeof, df[123]) == map(typeof, target)
        df_cols = Tables.columntable(df)
        @test mean(df_cols.uv_m) ≈ SVector(-174221.2u"m", 314413.6u"m")
        @test mean(df_cols.uv) === VLBI.UV(-8.926711f6, 1.6110279f7)
        @test mean(df_cols.visibility) ≈ 0.2495917 + 0.0010398296im
        # @test first(df) == first(table(uv, pyimport))
        # @test df == table(uv, pyimport)
    end

    @testset "multichannel" begin
        uv = VLBI.load(VLBI.UVData, "./data/vis_multichan.vis")
        @test length(uv.freq_windows) == 8
        df = table(uv, pyimport)
        target = (
            baseline = VLBIData.Baseline(1, (3, 5)),
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
        @test df[1234] == target
        @test map(typeof, df[1234]) == map(typeof, target)
        # @test df == table(uv, pyimport)
    end
end

# @testset "difmap model" begin
#     mod = VLBI.load("./data/difmap_model.mod")
#     @test length(components(mod)) == 4
#     @test map(flux, components(mod)) |> collect ≈ [0.522, 0.0217, 0.0176, 0.2145]u"Jy"  rtol=0.01
#     @test map(fwhm_max, components(mod)) |> collect ≈ [0.135, 0.463, 1.99, 0]u"mas"  rtol=0.01
#     @test coords(components(mod)[1]) ≈ [-0.000662, -0.00123]u"mas" rtol=0.01

#     mod_empty = VLBI.load("./data/difmap_model_empty.mod")
#     @test isempty(components(mod_empty))

#     mod_clean = VLBI.load("./data/difmap_model_clean.mod")
#     @test length(components(mod_clean)) == 631
#     @test isconcretetype(eltype(components(mod_clean)))
#     @test eltype(components(mod_clean)) <: Point

#     mod_map = VLBI.load(MultiComponentModel, "./data/map.fits")

#     tmpf = tempname()
#     VLBI.save(tmpf, mod)
#     @test VLBI.load(tmpf) == mod
#     VLBI.save(tmpf, mod_empty)
#     @test VLBI.load(tmpf) == mod_empty
#     VLBI.save(tmpf, mod_clean)
#     @test VLBI.load(tmpf) == mod_clean
#     VLBI.save(tmpf, mod_map)
#     @test VLBI.load(tmpf) ≈ mod_map  rtol=1e-4  # approx because saving to *.mod involves rounding; also Float32 vs 64
# end

# @testset "RFC" begin
#     return

#     @testset "vis" begin
#         rfci = RFC.Files()
#         @showprogress for f in rand(RFC.files(rfci, J2000="J1743-0350", suffix="vis"), 10)
#             try
#                 uv = VLBI.load(f)
#                 df = table(uv)
#             catch e
#                 @show f e
#                 rethrow()
#             end
#         end
#     end

#     @testset "maps" begin
#         rfci = RFC.Files()
#         @showprogress for f in rand(RFC.files(rfci, J2000="J1743-0350", suffix="map", extension="fits"), 10)
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
#                 beam(abspath(f))
#             catch e
#                 @show f e
#                 rethrow()
#             end
#         end
#     end
# end


# import CompatHelperLocal
# CompatHelperLocal.@check()

# import Aqua
# @testset "Aqua" begin
#     Aqua.test_ambiguities(VLBI, recursive=false)
#     Aqua.test_unbound_args(VLBI)
#     Aqua.test_undefined_exports(VLBI)
#     Aqua.test_stale_deps(VLBI)
# end
