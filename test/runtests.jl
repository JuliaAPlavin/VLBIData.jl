using Test
using Pkg.Artifacts
using DataFrames
using Statistics
using AxisKeys
using Unitful


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
