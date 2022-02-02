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


@testset "generic loading" begin
    @test VLBI.guess_type("./data/map.fits") == VLBI.FitsImage
    @test VLBI.guess_type("./data/vis.fits") == VLBI.UVData
    @test VLBI.guess_type("./data/vis_multichan.vis") == VLBI.UVData
    @test VLBI.guess_type("./data/difmap_model.mod") == MultiComponentModel
    @test VLBI.guess_type("./data/difmap_model_empty.mod") == MultiComponentModel
end


@testset "fits image" begin
    @testset "read nothing" begin
        img = VLBI.load("./data/map.fits", read_data=false)
        @test img.data === nothing

        @test VLBI.pixel_size(img) ≈ 0.2u"mas"  rtol=1e-5
        @test VLBI.pixel_steps(img) ≈ [-0.2u"mas", 0.2u"mas"]  rtol=1e-5
        @test VLBI.pixel_area(img) ≈ 0.04u"mas^2"  rtol=1e-5
        @test frequency(img) ≈ 4.344u"GHz"

        bm = beam(img)
        @test beam("./data/map.fits") === bm
        @test effective_area(bm) ≈ 6.5519451u"mas^2"  rtol=1e-5
        @test fwhm_max(bm) ≈ 4.02887356u"mas"  rtol=1e-5
        @test fwhm_min(bm) ≈ 1.43523228u"mas"  rtol=1e-5
        @test position_angle(bm) ≈ -0.046499863  rtol=1e-5
        @test intensity(bm)(SVector(0.5u"mas", 0u"mas")) ≈ 0.714721  rtol=1e-5
    end
    
    @testset "read data" begin
        img = VLBI.load("./data/map.fits", read_data=true)
        @test size(img.data) == (512, 512)
        @test mean(img.data) ≈ 2.4069064f-5
        @test maximum(img.data) ≈ 0.021357307f0
        @test img.data[123, 456] ≈ -7.372097f-5
        @test axiskeys(img.data, :ra)  .|> ustrip ≈ 51:-0.2:-51.2  atol=1e-3
        @test axiskeys(img.data, :dec) .|> ustrip ≈ -51.2:0.2:51   atol=1e-3
        @test axiskeys(img.data, :ra) isa AbstractRange
        @test axiskeys(img.data, :dec) isa AbstractRange
    end
    
    @testset "read clean" begin
        mod = VLBI.load(MultiComponentModel, "./data/map.fits")
        @test length(components(mod)) == 361
        @test sum(flux, components(mod)) ≈ 0.038659f0u"Jy"
        @test mean(first ∘ coords, components(mod)) ≈ -0.913573508654587u"mas"
        @test mean(last ∘ coords, components(mod)) ≈ 8.145706523588233u"mas"
    end
end

@testset "uvfits" begin
    @testset "simple" begin
        uv = VLBI.load(VLBI.UVData, "./data/vis.fits")
        @test uv.header.object == "J0000+0248"
        @test Date(uv.header) === Date(2016, 1, 3)
        @test uv.header.stokes == [:RR]
        @test frequency(uv.header) ≈ 4.128u"GHz"
        @test length(uv.freq_windows) == 8
        @test length(uv.ant_arrays) == 1
        antarr = only(uv.ant_arrays)
        @test antarr.name == "VLBA Correlator"
        @test map(a -> a.id, antarr.antennas) == 1:9
        @test map(a -> a.name, antarr.antennas) == [:BR, :FD, :HN, :KP, :LA, :MK, :NL, :OV, :SC]

        df = VLBI.table(uv)
        @test Tables.rowaccess(df)
        @test length(df) == 896
        @test all(∈(Tables.schema(df).names), [:uv, :visibility, :if_ix, :datetime])
        @test all(isconcretetype, Tables.schema(df).types)
        @test df[1].uv === VLBI.UV(3.0883878f7, 1.6280775f7)
        df_cols = Tables.columntable(df)
        @test mean(df_cols.uv_m) ≈ SVector(298060.56u"m", -363224.78u"m")
        @test mean(df_cols.uv) === VLBI.UV(4.318903f6, -5.2631365f6)
        @test mean(df_cols.visibility) ≈ 0.021919968 + 0.00062215974im
    end

    @testset "multichannel" begin
        uv = VLBI.load(VLBI.UVData, "./data/vis_multichan.vis")
        @test length(uv.freq_windows) == 8
        df = VLBI.table(uv)
    end
end

@testset "difmap model" begin
    mod = VLBI.load("./data/difmap_model.mod")
    @test length(components(mod)) == 3
    @test map(flux, components(mod)) |> collect ≈ [0.64, -0.01, 1.32e9]u"Jy"  rtol=0.01
    @test map(fwhm_max, components(mod)) |> collect ≈ [2.30, 6.04, 323e3]u"mas"  rtol=0.01
    @test coords(components(mod)[1]) == [-0.09183782420814114, 0.13039751573060954]u"mas"

    mod = VLBI.load("./data/difmap_model_empty.mod")
    @test isempty(components(mod))
end

@testset "RFC" begin
    return

    @testset "vis" begin
        rfci = RFC.Files()
        @showprogress for f in RFC.files(rfci, suffix="vis") .|> abspath
            try
                uv = VLBI.load(f)
                df = VLBI.table(uv)
            catch e
                @show f e
                rethrow()
            end
        end
    end

    @testset "maps" begin
        rfci = RFC.Files()
        @showprogress for f in RFC.files(rfci, suffix="map", extension="fits") .|> abspath
            try
                VLBI.load(f)
            catch e
                @show f e
                rethrow()
            end
            try
                VLBI.load(MultiComponentModel, f)
            catch e
                @show f e
                rethrow()
            end
        end
    end
end


import CompatHelperLocal
CompatHelperLocal.@check()

import Aqua
@testset "Aqua" begin
    Aqua.test_ambiguities(VLBI, recursive=false)
    Aqua.test_unbound_args(VLBI)
    Aqua.test_undefined_exports(VLBI)
    Aqua.test_stale_deps(VLBI)
end
