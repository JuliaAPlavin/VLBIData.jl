@testitem "visibilties" begin
    using Unitful
    using Dates
    using Uncertain
    using VLBIData.StructArrays

    uvtbl_orig = [
        (datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 1)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 0, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 2)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 0, 20), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 3)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 20, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 4)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5)), freq_spec=123, stokes=:LL, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5)), freq_spec=124, stokes=:LL, value=1±ᵤ0.1),
    ]

    @test VLBI.average_bytime(uvtbl_orig, 0u"minute")::StructArray == [
        (freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 1))),
        (freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 0, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 2))),
        (freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 0, 20), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 3))),
        (freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 20, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 4))),
        (freq_spec=123, stokes=:LL, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5))),
        (freq_spec=124, stokes=:LL, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5))),
    ]
    @test VLBI.average_bytime(uvtbl_orig, 10u"minute")::StructArray == [
        (freq_spec = 123, stokes = :RR, count = 2, value = 1±ᵤ0.07071067811865475, datetime = DateTime("2020-01-01T00:05:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 1.5))),
        (freq_spec = 123, stokes = :RR, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:05:00"), spec = VisSpec(Baseline(1, (1, 3)), UV(0, 3.))),
        (freq_spec = 123, stokes = :RR, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 4.))),
        (freq_spec = 123, stokes = :LL, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 5.))),
        (freq_spec = 124, stokes = :LL, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 5.))),
    ]
end

@testitem "scan functionality" begin
    using Unitful
    using Dates
    using Uncertain
    using VLBIData.StructArrays
    using VLBIData.IntervalSets

    # Create test data with clear scan gaps (>1 minute)
    uvtbl_scans = StructArray([
        # Scan 1: 3 observations within 30 seconds
        (datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 1)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 0, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 2)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 0, 30), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 3)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        
        # Gap of 5 minutes (> 1 minute threshold)
        
        # Scan 2: 2 observations within 20 seconds  
        (datetime=DateTime(2020, 1, 1, 0, 5, 30), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 4)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1),
        (datetime=DateTime(2020, 1, 1, 0, 5, 50), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 5)), freq_spec=123, stokes=:LL, value=1±ᵤ0.1),
        
        # Gap of 3 minutes (> 1 minute threshold)
        
        # Scan 3: single observation
        (datetime=DateTime(2020, 1, 1, 0, 8, 50), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 6)), freq_spec=124, stokes=:LL, value=1±ᵤ0.1),
    ])

    strategy = VLBI.GapBasedScans(min_gap=1u"minute")
    
    # Test add_scan_ids
    uvtbl_with_scans = VLBI.add_scan_ids(uvtbl_scans, strategy)
    @test uvtbl_with_scans.scan_id == [1, 1, 1, 2, 2, 3]
    
    # Test that other columns remain unchanged after adding scan_ids
    @test uvtbl_with_scans.datetime == uvtbl_scans.datetime
    
    # Test scan_intervals
    intervals = VLBI.scan_intervals(uvtbl_scans, strategy)
    @test length(intervals) == 3
    @test intervals[1] == DateTime(2020, 1, 1, 0, 0, 0)..DateTime(2020, 1, 1, 0, 0, 30)
    @test intervals[2] == DateTime(2020, 1, 1, 0, 5, 30)..DateTime(2020, 1, 1, 0, 5, 50)
    @test intervals[3] == DateTime(2020, 1, 1, 0, 8, 50)..DateTime(2020, 1, 1, 0, 8, 50)
    
    # Test with different gap threshold
    uvtbl_short_gaps = VLBI.add_scan_ids(uvtbl_scans, VLBI.GapBasedScans(min_gap=15u"s"))
    # Should create more scans due to shorter threshold
    # Gap between 0:00:10 and 0:00:30 is 20s > 15s threshold, so new scan
    # Gap between 0:05:30 and 0:05:50 is 20s > 15s threshold, so new scan
    @test uvtbl_short_gaps.scan_id == [1, 1, 2, 3, 4, 5]
    
    # Test with very large gap threshold
    uvtbl_long_gaps = VLBI.add_scan_ids(uvtbl_scans, VLBI.GapBasedScans(min_gap=10u"minute"))
    # Should create only one scan since all gaps are < 10 minutes
    @test uvtbl_long_gaps.scan_id == [1, 1, 1, 1, 1, 1]
end
