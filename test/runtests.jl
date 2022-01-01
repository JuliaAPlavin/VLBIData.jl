using Test
using Statistics
using AxisKeys
using Unitful
using Dates
import VLBIData as VLBI
using Tables


import CompatHelperLocal
CompatHelperLocal.@check()


@testset "utils" begin
    @test VLBI.fname_parse("J0000+0248_C_2016_01_03_pet_map.fits") == (J2000 = "J0000+0248", band = "C", epoch = "2016_01_03", author = "pet", suffix = "map", extension = "fits", directory=nothing)
    @test VLBI.fname_parse("J0000+0248_2016_01_03_pet") == (J2000 = "J0000+0248", band = nothing, epoch = "2016_01_03", author = "pet", suffix = nothing, extension = nothing, directory=nothing)
    @test VLBI.fname_parse("/home/J0000+0248_2016_01_03_pet") == (J2000 = "J0000+0248", band = nothing, epoch = "2016_01_03", author = "pet", suffix = nothing, extension = nothing, directory="/home/")
    @test VLBI.fname_build((J2000 = "J0000+0248", band = "C", epoch = "2016_01_03", author = "pet", suffix = "map", extension = "fits")) == "J0000+0248_C_2016_01_03_pet_map.fits"
    @test VLBI.fname_build((J2000 = "J0000+0248", band = "C", epoch = "2016_01_03", author = "pet")) == "J0000+0248_C_2016_01_03_pet"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis") == "J0000+0248_C_2016_01_03_pet_vis.fits"
    @test VLBI.fname_replace("/home/J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis") == "/home/J0000+0248_C_2016_01_03_pet_vis.fits"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis", extension=nothing) == "J0000+0248_C_2016_01_03_pet_vis"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="map"=>"vis") == "J0000+0248_C_2016_01_03_pet_vis.fits"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="map"=>"vis", extension=nothing) == "J0000+0248_C_2016_01_03_pet_vis"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="map"=>"vis", extension="fits.png") == "J0000+0248_C_2016_01_03_pet_vis.fits.png"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits-1.png", suffix="map"=>"vis", extension="fits-1.png" => "fits-2.png") == "J0000+0248_C_2016_01_03_pet_vis.fits-2.png"
    @test_throws AssertionError VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis"=>"map")
end

@testset "fits image" begin
    @testset "read nothing" begin
        img = VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=false, read_clean=false)
        @test img.clean_components === nothing
        @test img.data === nothing
        @test img.noise === nothing
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
    end
    
    @testset "read both" begin
        img = VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=true, read_clean=true)
        @test img.clean_components == VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=false, read_clean=true).clean_components
        @test img.data == VLBI.load(VLBI.FitsImage, "./data/map.fits", read_data=true, read_clean=false).data
    end
end

@testset "uvfits" begin
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

@testset "difmap_files" begin
    df = VLBI.load(VLBI.DifmapModel, "./data/difmap_model.mod")
    @test length(df) == 3
    @test map(x -> x.flux, df) ≈ [0.64, -0.01, 1.32e9]  rtol=0.01
end

import Aqua
@testset "Aqua" begin
    Aqua.test_ambiguities(VLBI, recursive=false)
    Aqua.test_unbound_args(VLBI)
    Aqua.test_undefined_exports(VLBI)
    Aqua.test_stale_deps(VLBI)
end
