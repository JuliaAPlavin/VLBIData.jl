abstract type AbstractScanStrategy end

@kwdef struct GapBasedScans <: AbstractScanStrategy
    min_gap = 1u"minute"
end

function scan_ids(strategy::GapBasedScans, uvtbl::StructVector)
    @modify(uvtbl |> sort(_, by=x->x.datetime)) do uvtbl
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
end

function add_scan_ids(strategy::AbstractScanStrategy, uvtbl)
    uvtbl = StructArray(uvtbl)
    @assert !hasproperty(uvtbl, :scan_id)
    @insert uvtbl.scan_id = scan_ids(strategy, uvtbl)
end

scan_intervals(uvtbl) = scan_intervals(nothing, uvtbl)
function scan_intervals(strategy::Union{Nothing,AbstractScanStrategy}, uvtbl)
    uvtbl = StructArray(uvtbl)
    uvtbl = if hasproperty(uvtbl, :scan_id)
        @p uvtbl sort(by=_.scan_id)
    else
        isnothing(strategy) && error("No scan strategy provided and no scan_id present in the table.")
        add_scan_ids(strategy, uvtbl)
    end
    @p let
        uvtbl
        group_vg(_.scan_id)
        map() do __
            extrema(_.datetime)
            Interval(__...)
        end
    end
end
