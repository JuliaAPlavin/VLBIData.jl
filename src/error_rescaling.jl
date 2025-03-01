@kwdef struct ConsecutiveDifferencesStandard
    maxΔt = nothing
    rayleigh_q = 0.5
    min_cnt = 50
end
ConsecutiveDifferencesStandard(maxΔt; kwargs...) = ConsecutiveDifferencesStandard(; maxΔt, kwargs...)

@kwdef struct CoherentAverageScatter
    maxΔt = nothing
    min_cnt_avg = 10
    min_cnt_after = 50
end
CoherentAverageScatter(maxΔt; kwargs...) = CoherentAverageScatter(; maxΔt, kwargs...)

struct ErrMulSame
    methods::Tuple
    rtol::Float64
end
ErrMulSame(methods...; rtol) = ErrMulSame(methods, rtol)


function find_errmul(uvtbl, m::CoherentAverageScatter)
    avg_t = @something(m.maxΔt, 15*typical_Δt(uvtbl))
    avg_t ∈ (30u"s"..16u"minute") || @warn "Strange Δt for averaging, consider specifying maxΔt explicitly" Δt=avg_t
    uvtbl_avg = @p compute_avgs(uvtbl; maxΔt=avg_t, min_cnt=m.min_cnt_avg) map(_.std / _.err)
    length(uvtbl_avg) ≥ m.min_cnt_after || return nothing
    return median(uvtbl_avg)
end

function find_errmul(uvtbl, m::ConsecutiveDifferencesStandard)
    diffs = compute_diffs(uvtbl; maxΔt=@something(m.maxΔt, 3*typical_Δt(uvtbl)))
    length(diffs) ≥ m.min_cnt || return nothing
    q = m.rayleigh_q
    @p diffs map(_.abs / _.err) quantile(__, q) / quantile(Rayleigh(), q)
end

function find_errmul(uvtbl, m::ErrMulSame)
    emuls = @p map(find_errmul(uvtbl, _), m.methods)
    if any(isnothing, emuls)
        error("""Rescaling factors for visibility errors couldn't be estimated by some methods.
        
        Obtained:
        $(@p zip(m.methods, emuls) map("$(_[1]) => $(_[2])") join(__, "\n"))
        """)
    elseif isapprox(minimum(emuls), maximum(emuls); m.rtol)
        return mean(emuls)
    else
        error("""Rescaling factors for visibility errors don't agree between methods.
        Requested rtol = $(m.rtol). Actual rel difference = $((maximum(emuls) - minimum(emuls)) / maximum(emuls)).
        
        Obtained:
        $(@p zip(m.methods, emuls) map("$(_[1]) => $(_[2])") join(__, "\n"))
        """)
    end
end


function typical_Δt(uvtbl)
    Δt = @p let
        uvtbl
        group_vg((;_.stokes, _.freq_spec, _.spec.bl))
        flatmap() do gr
            diff(gr.datetime)
        end
        map(_ |> u"minute" |> float)
        median
    end
    Δt ∈ (1u"s"..5u"minute") || @warn "Strange automatic median Δt" Δt
    Δt
end

compute_avgs(uvtbl; maxΔt, min_cnt) = @p let
    average_vis(uvtbl, maxΔt; avgvals=vals -> (err=U.uncertainty(U.weightedmean(vals)), std=std(U.value.(vals)) / √length(vals) / √2))
    # average_vis(uvtbl, avg_t; avgvals=vals -> (std=(mean(abs2∘value, vals) - abs2(mean(value, vals))) / (length(vals) - 1) / 2 |> sqrt, err=1/√(sum(v -> 1/uncertainty(v)^2, vals)) ))  # formula from difmap
    filter(_.count ≥ min_cnt)
    map(_.value)
end

compute_diffs(uvtbl; maxΔt) = @p let
    uvtbl
    group_vg((;_.stokes, _.freq_spec, _.spec.bl))
    flatmap() do gr
        Δs = @p let
            gr
            @aside @assert issorted(__, by=:datetime)
            zip(__[begin:end-1], __[begin+1:end])
            collect
            filter(_[2].datetime - _[1].datetime ≤ maxΔt)
            map(_[2].value - _[1].value)
        end
        map(Δs) do Δ
            (abs=abs(U.value(Δ)), err=U.uncertainty(Δ))
        end
    end
end

multiply_errors(uvtbl, errmul) =
    @modify(uvtbl[∗].value |> U.uncertainty) do err
        err * errmul
    end

function rescale_visibility_errors(uvtbl, m)
    mul = find_errmul(uvtbl, m)
    isnothing(mul) && error("Couldn't determine error rescaling factor using $m")
    multiply_errors(uvtbl, mul)
end
