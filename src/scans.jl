abstract type AbstractScanStrategy end

@kwdef struct GapBasedScans <: AbstractScanStrategy
    min_gap = 1u"minute"
end

@stable function scan_ids(strategy::GapBasedScans, uvtbl::StructVector)
    min_gap::typeof(1.0u"s") = strategy.min_gap
    scan_id_ref = Ref(1)
    prev_dt_ref = Ref(minimum(uvtbl.datetime))
    @modify($(uvtbl.datetime) |> sort) do datetimes
        map(datetimes) do dt
            gap = dt - prev_dt_ref[]
            if gap > min_gap
                scan_id_ref[] += 1
            end
            prev_dt_ref[] = dt
            scan_id_ref[]
        end
    end
end

function add_scan_ids(strategy::AbstractScanStrategy, uvtbl)
    uvtbl = StructArray(uvtbl)
    @assert !hasproperty(uvtbl, :scan_id)
    @insert uvtbl.scan_id = scan_ids(strategy, uvtbl)
end

@stable scan_intervals(uvtbl) = scan_intervals(nothing, uvtbl)
scan_intervals(strategy, uvtbl) = scan_intervals(strategy, StructArray(uvtbl))
@stable function scan_intervals(strategy::Nothing, uvtbl::StructArray)
    hasproperty(uvtbl, :scan_id) || error("No scan strategy provided and no scan_id present in the table.")
    @p let
        uvtbl
        groupview(_.scan_id; restype=Vector)
        map() do gr
            Interval(extrema(gr.datetime)...)
        end
    end
end
@stable function scan_intervals(strategy::AbstractScanStrategy, uvtbl::StructArray)
    uvtbl = add_scan_ids(strategy, uvtbl)
    @p let
        uvtbl
        groupview(_.scan_id; restype=Vector)
        map() do gr
            Interval(extrema(gr.datetime)...)
        end
    end
end
