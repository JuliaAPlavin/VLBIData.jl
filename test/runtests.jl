using Test
using Pkg.Artifacts
using DataFrames
using Statistics
using AxisKeys
using Unitful
using Dates


import VLBIData
const VLBI = VLBIData


let
    artifact_toml = joinpath(@__DIR__, "Artifacts.toml")
    art_hash = artifact_hash("test_data", artifact_toml)
    if art_hash == nothing || !artifact_exists(art_hash)
        @info "Creating artifact"
        art_hash = create_artifact() do artifact_dir
            for url in [
                    "http://astrogeo.org/images/J0000+0248/J0000+0248_C_2016_01_03_pet_map.fits",
                    "http://astrogeo.org/images/J0000+0248/J0000+0248_C_2016_01_03_pet_vis.fits",
                ]
                download(url, joinpath(artifact_dir, basename(url)))
            end
        end
        bind_artifact!(artifact_toml, "test_data", art_hash, force=true)
    end
end


@testset "utils" begin
    @test VLBI.fname_parse("J0000+0248_C_2016_01_03_pet_map.fits") == (J2000 = "J0000+0248", band = "C", epoch = "2016_01_03", author = "pet", suffix = "map", extension = "fits")
    @test VLBI.fname_build((J2000 = "J0000+0248", band = "C", epoch = "2016_01_03", author = "pet", suffix = "map", extension = "fits")) == "J0000+0248_C_2016_01_03_pet_map.fits"
    @test VLBI.fname_build((J2000 = "J0000+0248", band = "C", epoch = "2016_01_03", author = "pet")) == "J0000+0248_C_2016_01_03_pet"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis") == "J0000+0248_C_2016_01_03_pet_vis.fits"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis", extension=nothing) == "J0000+0248_C_2016_01_03_pet_vis"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="map"=>"vis") == "J0000+0248_C_2016_01_03_pet_vis.fits"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="map"=>"vis", extension=nothing) == "J0000+0248_C_2016_01_03_pet_vis"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="map"=>"vis", extension="fits.png") == "J0000+0248_C_2016_01_03_pet_vis.fits.png"
    @test VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits-1.png", suffix="map"=>"vis", extension="fits-1.png" => "fits-2.png") == "J0000+0248_C_2016_01_03_pet_vis.fits-2.png"
    @test_throws AssertionError VLBI.fname_replace("J0000+0248_C_2016_01_03_pet_map.fits", suffix="vis"=>"map")
end

@testset "fits image" begin
    @testset "read nothing" begin
        img = VLBI.load(VLBI.FitsImage, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_map.fits"), read_data=false, read_clean=false)
        @test img.clean_components === nothing
        @test img.data === nothing
        @test img.noise === nothing
    end
    
    @testset "read clean" begin
        img = VLBI.load(VLBI.FitsImage, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_map.fits"), read_data=false, read_clean=true)
        @test img.data === nothing
        @test img.noise === nothing
        @test names(img.clean_components) == ["flux", "ra", "dec"]
        @test (<:).(eltype.(eachcol(img.clean_components)), AbstractFloat) |> all
        @test nrow(img.clean_components) == 361
        @test sum(img.clean_components[!, :flux]) ≈ 0.038659f0
        @test mean(img.clean_components[!, :ra]) ≈ -0.913573508654587
        @test mean(img.clean_components[!, :dec]) ≈ 8.145706523588233
    end
    
    @testset "read data" begin
        img = VLBI.load(VLBI.FitsImage, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_map.fits"), read_data=true, read_clean=false)
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
        img = VLBI.load(VLBI.FitsImage, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_map.fits"), read_data=true, read_clean=true)
        @test img.clean_components == VLBI.load(VLBI.FitsImage, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_map.fits"), read_data=false, read_clean=true).clean_components
        @test img.data == VLBI.load(VLBI.FitsImage, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_map.fits"), read_data=true, read_clean=false).data
    end
end

@testset "uvfits" begin
    uv = VLBI.load(VLBI.UVData, joinpath(artifact"test_data", "J0000+0248_C_2016_01_03_pet_vis.fits"))
    @test uv.header.object == "J0000+0248"
    @test uv.header.date_obs === Date(2016, 1, 3)
    @test uv.header.stokes == [:RR]
    @test uv.header.frequency ≈ 4.128u"GHz"
    @test length(uv.freq_windows) == 8
    df = VLBI.read_data_dataframe(uv)
    @test nrow(df) == 896
    @test ncol(df) == 16
    @test all(["u", "v", "w", "visibility", "iif"] .∈ Ref(names(df)))
    @test isconcretetype.(eltype.(eachcol(df))) |> all
    @test mean(df[!, :u]) ≈ 298060.56u"m"
    @test mean(df[!, :v_wl]) ≈ -5.2631365e6
    @test mean(df[!, :visibility]) ≈ 0.021919968 + 0.00062215974im
end
