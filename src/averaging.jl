aggspec(bl, specs::AbstractVector{<:VisSpec}) = VisSpec(bl, mean(x->UV(x), specs) |> UV)

function average_data(strategy::AbstractScanStrategy, uvtbl; avgvals=U.weightedmean)
    uvtbl_with_scans = add_scan_ids(strategy, uvtbl)
    const_part = @p getproperties(uvtbl_with_scans) (@delete __[(:freq_spec, :stokes, :scan_id, :value, :spec)]) _filter(allequal) map(uniqueonly)
    return _average_data_by_scans(uvtbl_with_scans, const_part; avgvals)
end

function _average_data_by_scans(uvtbl_with_scans, const_part; avgvals)
    intervals = scan_intervals(uvtbl_with_scans)
    @p begin
        uvtbl_with_scans
        groupview_vg((;bl=Baseline(_), _.freq_spec, _.stokes, _.scan_id))
        map((;
            const_part...,
            delete(key(_), @o _.bl)...,
            count = length(_),
            value = avgvals(_.value),
            datetime = _mean(intervals[key(_).scan_id]),
            spec = aggspec(key(_).bl, _.spec),
        ))
        StructArray()
    end
end

@kwdef struct FixedTimeIntervals
    interval = 1u"minute"
end

function average_data(strategy::FixedTimeIntervals, src_; avgvals=U.weightedmean)
    src = StructArray(src_)
    const_part = @p getproperties(src) (@delete __[(:freq_spec, :stokes, :value, :spec)]) _filter(allequal) map(uniqueonly)
    return _average_data_by_time(strategy, src, const_part; avgvals)
end

function _average_data_by_time(strategy::FixedTimeIntervals, src, const_part; avgvals)
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
            const_part...,
            delete(key(_), @o _.bl _.timestep)...,
            count = length(_),
            value = avgvals(_.value),
            datetime = mindt + (key(_).timestep + 0.5) * avg_interval,
            spec = aggspec(key(_).bl, _.spec),
        ))
        StructArray()
    end
end

_mean(i::Interval) = minimum(i) + (maximum(i) - minimum(i)) ÷ 2
if VERSION ≥ v"1.11"
    _filter(pred, x) = filter(pred, x)
else
    _filter(f, xs::NamedTuple) = xs[filter(k -> f(xs[k]), keys(xs))]
end
