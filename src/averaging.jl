@deprecate average_bytime(src, avg_interval; avgvals=U.weightedmean) average_data(src, FixedTimeIntervals(avg_interval); avgvals)

aggspec(bl, specs::AbstractVector{<:VisSpec}) = VisSpec(bl, mean(x->UV(x), specs) |> UV)

function average_data(uvtbl, strategy::AbstractScanStrategy; avgvals=U.weightedmean)
    uvtbl_with_scans = add_scan_ids(uvtbl, strategy)
    intervals = scan_intervals(uvtbl_with_scans)
what
    @p begin
        uvtbl_with_scans
        groupview_vg((;bl=Baseline(_), _.freq_spec, _.stokes, _.scan_id))
        map((;
            delete(key(_), @o _.bl _.scan_id)...,
            count = length(_),
            value = avgvals(_.value),
            datetime = _mean(intervals[key(_).scan_id]),
            spec = aggspec(key(_).bl, _.spec),
            scan_id = key(_).scan_id
        ))
        StructArray()
    end
end

@kwdef struct FixedTimeIntervals
    interval = 1u"minute"
end

function average_data(src_, strategy::FixedTimeIntervals; avgvals=U.weightedmean)
    src = StructArray(src_)
    mindt = minimum(src.datetime)
    avg_interval = @p let
        float(strategy.interval)
        iszero(__) ? oftype(__, 1u"ms") : __
    end
    @p begin
        src
        groupview_vg((;
            bl=Baseline(_), _.freq_spec, _.stokes,
            timestep = ((_.datetime - mindt) / avg_interval) |> upreferred |> x->trunc(Int, x),
        ))
        map((;
            delete(key(_), @o _.bl _.timestep)...,
            count = length(_),
            value = avgvals(_.value),
            datetime = mindt + (key(_).timestep + 0.5) * avg_interval,
            spec = aggspec(key(_).bl, _.spec)
        ))
        StructArray()
    end
end

_mean(i::Interval) = minimum(i) + (maximum(i) - minimum(i)) ÷ 2
