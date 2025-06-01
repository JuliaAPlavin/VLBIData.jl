@kwdef struct GapBasedScans
    min_gap = 1u"minute"
end

function scan_ids(uvtbl::StructVector, strategy::GapBasedScans)
    uvtbl = @p uvtbl sort(by=_.datetime)
    scan_id_ref = Ref(1)
    prev_dt_ref = Ref(uvtbl.datetime[1])
    scan_id = map(uvtbl.datetime) do dt
        gap = dt - prev_dt_ref[]
        if gap > strategy.min_gap
            scan_id_ref[] += 1
        end
        prev_dt_ref[] = dt
        scan_id_ref[]
    end
    return scan_id
end

function add_scan_ids(uvtbl, strategy)
    uvtbl = StructArray(uvtbl)
    @insert uvtbl.scan_id = scan_ids(uvtbl, strategy)
end

scan_intervals(uvtbl, strategy) = @p let
    add_scan_ids(uvtbl, strategy)
    group_vg(_.scan_id)
    map() do __
        extrema(_.datetime)
        Interval(__...)
    end
end
