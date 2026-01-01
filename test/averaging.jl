@testitem "visibilties" begin
    using Unitful
    using Dates
    using Uncertain
    using VLBIData.StructArrays

    uvtbl_orig = [
        (datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 1)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1, source="SgrA*"),
        (datetime=DateTime(2020, 1, 1, 0, 0, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 2)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1, source="SgrA*"),
        (datetime=DateTime(2020, 1, 1, 0, 0, 20), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 3)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1, source="SgrA*"),
        (datetime=DateTime(2020, 1, 1, 0, 20, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 4)), freq_spec=123, stokes=:RR, value=1±ᵤ0.1, source="SgrA*"),
        (datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5)), freq_spec=123, stokes=:LL, value=1±ᵤ0.1, source="SgrA*"),
        (datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5)), freq_spec=124, stokes=:LL, value=1±ᵤ0.1, source="SgrA*"),
    ]

    @test VLBI.average_bytime(uvtbl_orig, 0u"minute")::StructArray == [
        (source="SgrA*", freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 1))),
        (source="SgrA*", freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 0, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 2))),
        (source="SgrA*", freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 0, 20), spec=VisSpec(Baseline(1, (1, 3)), UV(0, 3))),
        (source="SgrA*", freq_spec=123, stokes=:RR, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 20, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 4))),
        (source="SgrA*", freq_spec=123, stokes=:LL, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5))),
        (source="SgrA*", freq_spec=124, stokes=:LL, count=1, value=1±ᵤ0.1, datetime=DateTime(2020, 1, 1, 0, 20, 10), spec=VisSpec(Baseline(1, (1, 2)), UV(0, 5))),
    ]
    @test VLBI.average_bytime(uvtbl_orig, 10u"minute")::StructArray == [
        (source="SgrA*", freq_spec = 123, stokes = :RR, count = 2, value = 1±ᵤ0.07071067811865475, datetime = DateTime("2020-01-01T00:05:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 1.5))),
        (source="SgrA*", freq_spec = 123, stokes = :RR, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:05:00"), spec = VisSpec(Baseline(1, (1, 3)), UV(0, 3.))),
        (source="SgrA*", freq_spec = 123, stokes = :RR, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 4.))),
        (source="SgrA*", freq_spec = 123, stokes = :LL, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 5.))),
        (source="SgrA*", freq_spec = 124, stokes = :LL, count = 1, value = 1.0 ±ᵤ 0.1, datetime = Dates.DateTime("2020-01-01T00:25:00"), spec = VisSpec(Baseline(1, (1, 2)), UV(0, 5.))),
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
    
    # Test scan_intervals with strategy=nothing on data that already has scan_id
    intervals_from_existing = VLBI.scan_intervals(uvtbl_with_scans, nothing)
    @test intervals_from_existing == intervals  # Should be same as when computed with strategy
    
    # Test that scan_intervals results are always in scan order (1, 2, 3, ...) even with unordered scan_ids
    uvtbl_unordered = StructArray([
        (datetime=DateTime(2020, 1, 1, 0, 5, 0), scan_id=2),   # Scan 2 first
        (datetime=DateTime(2020, 1, 1, 0, 0, 0), scan_id=1),   # Scan 1 second
        (datetime=DateTime(2020, 1, 1, 0, 0, 30), scan_id=1),  # More scan 1
        (datetime=DateTime(2020, 1, 1, 0, 5, 30), scan_id=2),  # More scan 2
        (datetime=DateTime(2020, 1, 1, 0, 10, 0), scan_id=3),  # Scan 3 last
    ])
    intervals_unordered = VLBI.scan_intervals(uvtbl_unordered, nothing)
    @test length(intervals_unordered) == 3
    # Should be ordered by scan_id: scan 1, scan 2, scan 3
    @test intervals_unordered[1] == DateTime(2020, 1, 1, 0, 0, 0)..DateTime(2020, 1, 1, 0, 0, 30)   # Scan 1
    @test intervals_unordered[2] == DateTime(2020, 1, 1, 0, 5, 0)..DateTime(2020, 1, 1, 0, 5, 30)   # Scan 2  
    @test intervals_unordered[3] == DateTime(2020, 1, 1, 0, 10, 0)..DateTime(2020, 1, 1, 0, 10, 0)  # Scan 3
    
    # Test that scan_intervals with strategy=nothing fails on data without scan_id
    @test_throws ErrorException VLBI.scan_intervals(uvtbl_scans, nothing)
end

@testitem "average_data by scans" begin
    using Unitful
    using Dates
    using Uncertain
    using VLBIData.StructArrays

    # Create test data with multiple observations per scan for the same baseline/freq/stokes
    uvtbl_for_avg = StructArray([
        # Scan 1: observations in chronological order
        (datetime=DateTime(2020, 1, 1, 0, 0, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(1, 1)), freq_spec=123, stokes=:RR, value=2±ᵤ0.1, source="M87"),
        (datetime=DateTime(2020, 1, 1, 0, 0, 15), spec=VisSpec(Baseline(1, (1, 3)), UV(5, 5)), freq_spec=123, stokes=:RR, value=10±ᵤ0.1, source="M87"),
        (datetime=DateTime(2020, 1, 1, 0, 0, 30), spec=VisSpec(Baseline(1, (1, 2)), UV(2, 2)), freq_spec=123, stokes=:RR, value=4±ᵤ0.1, source="M87"),
        (datetime=DateTime(2020, 1, 1, 0, 0, 45), spec=VisSpec(Baseline(1, (1, 2)), UV(6, 6)), freq_spec=123, stokes=:LL, value=12±ᵤ0.1, source="M87"),
        
        # Gap of 5 minutes (> 1 minute threshold)
        
        # Scan 2: 2 observations for same baseline/freq/stokes - should be averaged
        (datetime=DateTime(2020, 1, 1, 0, 5, 0), spec=VisSpec(Baseline(1, (1, 2)), UV(3, 3)), freq_spec=123, stokes=:RR, value=6±ᵤ0.1, source="M87"),
        (datetime=DateTime(2020, 1, 1, 0, 5, 30), spec=VisSpec(Baseline(1, (1, 2)), UV(4, 4)), freq_spec=123, stokes=:RR, value=8±ᵤ0.1, source="M87"),
    ])

    strategy = VLBI.GapBasedScans(min_gap=1u"minute")
    
    # Test average_data
    averaged = VLBI.average_data(uvtbl_for_avg, strategy)
    
    # Should have 4 rows: baseline (1,2) RR in scan 1, baseline (1,2) RR in scan 2, baseline (1,3) RR in scan 1, baseline (1,2) LL in scan 1
    @test length(averaged) == 4
    
    # Test explicit scan_id values
    @test averaged.scan_id == [1, 1, 1, 2]
    
    # Test explicit datetime values (order should be deterministic based on grouping)
    @test averaged.datetime == [
        DateTime(2020, 1, 1, 0, 0, 22, 500),
        DateTime(2020, 1, 1, 0, 0, 22, 500),
        DateTime(2020, 1, 1, 0, 0, 22, 500),
        DateTime(2020, 1, 1, 0, 5, 15, 0),
    ]
    
    # Test explicit value array (averaged values and uncertainties)
    @test U.value.(averaged.value) ≈ [
        3.0,   # (1,2) RR scan1 - weighted mean of 2±0.1 and 4±0.1
        10.0,  # (1,3) RR scan1 - single value 10±0.1
        12.0,  # (1,2) LL scan1 - single value 12±0.1
        7.0    # (1,2) RR scan2 - weighted mean of 6±0.1 and 8±0.1
    ]
    
    @test U.uncertainty.(averaged.value) ≈ [
        √2 * 0.1 / 2,      # (1,2) RR scan1 - uncertainty from weighted mean
        0.1,                # (1,2) LL scan1 - single uncertainty
        0.1,                # (1,3) RR scan1 - single uncertainty  
        √2 * 0.1 / 2       # (1,2) RR scan2 - uncertainty from weighted mean
    ]
    
    # Test stokes parameters
    @test averaged.stokes == [
        :RR,
        :RR,
        :LL,
        :RR,
    ]
    
    # Test that source column is maintained
    @test averaged.source == [
        "M87",  # (1,2) RR scan1
        "M87",  # (1,3) RR scan1
        "M87",  # (1,2) LL scan1
        "M87"   # (1,2) RR scan2
    ]
end
