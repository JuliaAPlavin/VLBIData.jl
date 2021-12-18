using Test
using Statistics
using AxisKeys
using Unitful, UnitfulAstro, MyUnitful
using Dates
import VLBIData as VLBI
using Tables


@testset "fits image" begin
    @testset "read nothing" begin
        img = VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=false, read_clean=false)
        @test img.clean_components === nothing
        @test img.data === nothing
        @test img.noise === nothing

        @test VLBI.pixel_size(img) ≈ 0.2u"mas"  rtol=1e-5
        @test VLBI.pixel_steps(img) ≈ [-0.2u"mas", 0.2u"mas"]  rtol=1e-5
        @test VLBI.pixel_area(img) ≈ 0.04u"mas^2"  rtol=1e-5
        @test VLBI.frequency(img) ≈ 4.344u"GHz"
        beam = VLBI.image_beam(img)
        @test VLBI.area(beam) ≈ 4.54146231u"mas^2"  rtol=1e-5
        @test beam.major_axis ≈ 4.02887356u"mas"  rtol=1e-5
        @test beam.minor_axis ≈ 1.43523228u"mas"  rtol=1e-5
        @test beam.pa ≈ -0.046499863  rtol=1e-5
        @test beam([0.5u"mas", 0u"mas"]) ≈ 0.714721  rtol=1e-5
    end
    
    @testset "read clean" begin
        img = VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=false, read_clean=true)
        @test img.data === nothing
        @test img.noise === nothing
        @test all(∈(Tables.schema(img.clean_components).names), [:flux, :radec])
        @test all(isconcretetype, Tables.schema(img.clean_components).types)
        @test length(img.clean_components) == 361
        @test sum(c -> c.flux, img.clean_components) ≈ 0.038659f0
        @test mean(c -> c.radec[1], img.clean_components) ≈ -0.913573508654587
        @test mean(c -> c.radec[2], img.clean_components) ≈ 8.145706523588233
    end
    
    @testset "read data" begin
        img = VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=true, read_clean=false)
        @test img.clean_components === nothing
        @test img.noise ≈ 9.5605465f-5
        @test size(img.data) == (512, 512)
        @test mean(img.data) ≈ 2.4069064f-5
        @test maximum(img.data) ≈ 0.021357307f0
        @test img.data[123, 456] ≈ -7.372097f-5
        @test axiskeys(img.data, :ra)  .|> ustrip ≈ 51:-0.2:-51.2  atol=1e-3
        @test axiskeys(img.data, :dec) .|> ustrip ≈ -51.2:0.2:51   atol=1e-3
        @test axiskeys(img.data, :ra) isa AbstractRange
        @test axiskeys(img.data, :dec) isa AbstractRange
    end
    
    @testset "read both" begin
        img = VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=true, read_clean=true)
        @test img.clean_components == VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=false, read_clean=true).clean_components
        @test img.data == VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=true, read_clean=false).data
    end
end

@testset "uvfits" begin
    @testset "simple" begin
        uv = VLBI.load(VLBI.UVData, "./data/vis.fits")
        @test uv.header.object == "J0000+0248"
        @test uv.header.date_obs === Date(2016, 1, 3)
        @test uv.header.stokes == [:RR]
        @test uv.header.frequency ≈ 4.128u"GHz"
        @test length(uv.freq_windows) == 8
        @test length(uv.ant_arrays) == 1
        antarr = only(uv.ant_arrays)
        @test antarr.name == "VLBA Correlator"
        @test map(a -> a.id, antarr.antennas) == 1:9
        @test map(a -> a.name, antarr.antennas) == [:BR, :FD, :HN, :KP, :LA, :MK, :NL, :OV, :SC]

        df = VLBI.read_data_table(uv)
        @test Tables.rowaccess(df)
        @test length(df) == 896
        @test length(Tables.schema(df).names) == 16
        @test all(∈(Tables.schema(df).names), [:u, :v, :w, :visibility, :iif])
        @test all(isconcretetype, Tables.schema(df).types)
        df_cols = Tables.columntable(df)
        @test mean(df_cols.u) ≈ 298060.56u"m"
        @test mean(df_cols.v_wl) ≈ -5.2631365e6
        @test mean(df_cols.visibility) ≈ 0.021919968 + 0.00062215974im
    end

    @testset "multichannel" begin
        uv = VLBI.load(VLBI.UVData, "./data/vis_multichan.vis")
        @test length(uv.freq_windows) == 8
    end
end

@testset "difmap_files" begin
    df = @inferred VLBI.load(VLBI.DifmapModel, "./data/difmap_model.mod")
    @test length(df) == 3
    @test map(x -> x.flux, df) ≈ [0.64, -0.01, 1.32e9]  rtol=0.01
    @test df[1].radec == [-0.09183782420814114, 0.13039751573060954]

    df_t = VLBI.load(VLBI.DifmapModel, "./data/difmap_model.mod"; radec_type=NTuple{2, Float64})

    df_e = @inferred VLBI.load(VLBI.DifmapModel, "./data/difmap_model_empty.mod")
    @test isempty(df_e)
    @test typeof(df_e) == typeof(df)
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
