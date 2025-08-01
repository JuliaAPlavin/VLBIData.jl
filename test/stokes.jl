@testitem "is_parallel_hands" begin
    # Test parallel hands (same polarization)
    @test is_parallel_hands(:RR) == true
    @test is_parallel_hands(:XX) == true
    
    # Test cross hands (different polarizations)
    @test is_parallel_hands(:RL) == false
    @test is_parallel_hands(:XY) == false
    
    # Test error cases
    @test_throws ErrorException is_parallel_hands(:I)
end

@testitem "is_cross_hands" begin
    # Test parallel hands (should be false)
    @test is_cross_hands(:RR) == false
    @test is_cross_hands(:XX) == false
    
    # Test cross hands (should be true)
    @test is_cross_hands(:RL) == true
    @test is_cross_hands(:XY) == true
    
    # Test error cases
    @test_throws ErrorException is_cross_hands(:I)
end

@testitem "stokes_to_feeds" begin
    # Test circular polarization
    @test stokes_to_feeds(:RR) == (:R, :R)
    @test stokes_to_feeds(:RL) == (:R, :L)
    
    # Test linear polarization
    @test stokes_to_feeds(:XX) == (:X, :X)
    @test stokes_to_feeds(:XY) == (:X, :Y)
    
    # Test mixed polarizations
    @test stokes_to_feeds(:RX) == (:R, :X)
    @test stokes_to_feeds(:LY) == (:L, :Y)
    
    # Test error cases
    @test_throws ErrorException stokes_to_feeds(:I)
end

@testitem "Consistency between is_parallel_hands and is_cross_hands" begin
    # Test that parallel and cross hands are mutually exclusive for valid symbols
    valid_symbols = [:RR, :LL, :XX, :YY, :RL, :LR, :XY, :YX]
    for s in valid_symbols
        @test is_parallel_hands(s) != is_cross_hands(s)
    end
end
