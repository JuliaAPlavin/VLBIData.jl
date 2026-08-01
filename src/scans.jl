abstract type AbstractScanStrategy end

@kwdef struct GapBasedScans <: AbstractScanStrategy
    min_gap = 1u"minute"
end

@stable function scan_ids(strategy::GapBasedScans, uvtbl::StructVector)
    min_gap::typeof(1.0u"s") = strategy.min_gap
    order = sortperm(uvtbl.datetime)
    rows = view(uvtbl, order)   # rows in time order (view, no copy)
    ids = similar(order)
    id = 1
    for (k, i) in enumerate(order)
        # a new scan starts at a time gap over the threshold, or — since a scan can never span
        # two sources — at a source change (even within the gap threshold)
        newscan = k > 1 && (rows.datetime[k] - rows.datetime[k-1] > min_gap ||
                            (hasproperty(uvtbl, :source) && rows.source[k] != rows.source[k-1]))
        newscan && (id += 1)
        ids[i] = id
    end
    ids
end

function add_scan_ids(strategy::AbstractScanStrategy, uvtbl)
    uvtbl = convert(StructArray, uvtbl)  # no-op if already a StructArray — avoids copying (materializing) lazy columns
    @assert !hasproperty(uvtbl, :scan_id)
    @insert uvtbl.scan_id = scan_ids(strategy, uvtbl)
end

@stable scan_intervals(uvtbl) = scan_intervals(nothing, uvtbl)
scan_intervals(strategy, uvtbl) = scan_intervals(strategy, convert(StructArray, uvtbl))
@stable function scan_intervals(strategy::Nothing, uvtbl::StructArray)
    hasproperty(uvtbl, :scan_id) || error("No scan strategy provided and no scan_id present in the table.")
    @p let
        uvtbl
        groupview(_.scan_id; into=Vector)
        map() do gr
            Interval(extrema(gr.datetime)...)
        end
    end
end
@stable function scan_intervals(strategy::AbstractScanStrategy, uvtbl::StructArray)
    uvtbl = add_scan_ids(strategy, uvtbl)
    @p let
        uvtbl
        groupview(_.scan_id; into=Vector)
        map() do gr
            Interval(extrema(gr.datetime)...)
        end
    end
end
