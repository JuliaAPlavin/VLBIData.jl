"""
    ConsecutiveDifferencesStandard(; maxΔt=5u"minute", rayleigh_q=0.5, min_cnt=50)

Estimate visibility error rescaling from consecutive visibility differences on each baseline.

For each baseline, takes pairs of adjacent-in-time visibilities (within `maxΔt`) and computes
`|ΔV| / σ_ΔV` where `ΔV = V(t₂) - V(t₁)` and `σ_ΔV = √(σ₁² + σ₂²)`.
If errors are correct, `|ΔV| / σ_ΔV` follows a Rayleigh distribution.
The rescaling factor is `quantile(|ΔV|/σ_ΔV, q) / quantile(Rayleigh, q)`.

Assumes the true visibility is constant between adjacent points.
Not robust to station gain drifts between adjacent points: consider [`ClosurePhaseConsecutive`](@ref) for sparsely-sampled data.
"""
@kwdef struct ConsecutiveDifferencesStandard
    maxΔt = 5u"minute"
    rayleigh_q = 0.5
    min_cnt = 50
end
ConsecutiveDifferencesStandard(maxΔt; kwargs...) = ConsecutiveDifferencesStandard(; maxΔt, kwargs...)

"""
    CoherentAverageScatter(; maxΔt=5u"minute", min_cnt_avg=10, min_cnt_after=50)

Estimate visibility error rescaling by comparing scatter within time bins to propagated errors.

Coherently averages visibilities in time bins of length `maxΔt`, then computes
`std(V) / (√N · √2)` (empirical scatter of the mean) vs the propagated weighted-mean error.
The rescaling factor is `median(scatter / propagated_error)` across all bins.

Assumes the true visibility is constant within each bin.
Not robust to station gain drifts: consider [`ClosurePhaseConsecutive`](@ref) if gain differences within bins are significant.
"""
@kwdef struct CoherentAverageScatter
    maxΔt = 5u"minute"
    min_cnt_avg = 10
    min_cnt_after = 50
end
CoherentAverageScatter(maxΔt; kwargs...) = CoherentAverageScatter(; maxΔt, kwargs...)

"""
    ClosurePhaseConsecutive(; maxΔt=30u"minute", min_cnt=10, min_triangles=5)

Estimate visibility error rescaling from consecutive closure phase differences.

For each station triangle, computes closure phases `φ = arg(V₁₂·V₂₃·V₃₁)` and their
uncertainties `σ_φ = σ_bispectrum / |bispectrum|` at each time step.
Then takes studentized consecutive differences: `s = Δφ / √(σ₁² + σ₂²)`.
If errors are correct, `s` is standard normal, so `std(s) ≈ 1`.
The rescaling factor is `median(std(s))` across all triangles.

Robust to station-based gain drifts (closure phases cancel station gains by construction).
Based on the approach in eht-imaging (`estimate_noise_rescale_factor`), see EHT Memo 2019-CE-02.
"""
@kwdef struct ClosurePhaseConsecutive
    maxΔt = 30u"minute"
    min_cnt = 10
    min_triangles = 5
end
ClosurePhaseConsecutive(maxΔt; kwargs...) = ClosurePhaseConsecutive(; maxΔt, kwargs...)

"""
    ErrMulSame(methods...; rtol)

Run multiple error rescaling methods and return their mean, requiring they agree within `rtol`.
Errors if any method returns `nothing` or if the results disagree beyond `rtol`.
"""
struct ErrMulSame
    methods::Tuple
    rtol::Float64
end
ErrMulSame(methods...; rtol) = ErrMulSame(methods, rtol)


function find_errmul(m::CoherentAverageScatter, uvtbl)
    avg_t = min(m.maxΔt, 15*typical_Δt(uvtbl))
    30u"s" < avg_t < 16u"minute" || @warn "Strange Δt for averaging, consider specifying maxΔt explicitly" Δt=avg_t
    uvtbl_avg = @p compute_avgs(uvtbl; maxΔt=avg_t, min_cnt=m.min_cnt_avg) map(_.std / _.err)
    length(uvtbl_avg) ≥ m.min_cnt_after || return nothing
    return median(uvtbl_avg)
end

function find_errmul(m::ConsecutiveDifferencesStandard, uvtbl)
    diffs = compute_diffs(uvtbl; maxΔt=min(m.maxΔt, 3*typical_Δt(uvtbl)))
    length(diffs) ≥ m.min_cnt || return nothing
    q = m.rayleigh_q
    @p diffs map(_.abs / _.err) quantile(__, q) / quantile_rayleigh(q)
end

# same as quantile(Rayleigh(), q), but dependency-free
quantile_rayleigh(q) = √(-2 * log(1 - q))

function find_errmul(m::ClosurePhaseConsecutive, uvtbl)
    maxΔt = min(m.maxΔt, 5*typical_Δt(uvtbl))
    cps = closures_all(ClosurePhaseSpec, uvtbl)
    isempty(cps) && return nothing
    NT = intersect_nt_type(eltype(cps), NamedTuple{(:freq_spec, :stokes)})
    std_list = @p let
        cps
        group_vg((; NT(_)..., tri=antenna_names(_.spec)))
        filtermap() do gr
            @assert issorted(gr, by=:datetime)
            s_list = @p let
                gr
                map(angle(_.value))
                zip(__[begin:end-1], __[begin+1:end])
                zip(__, diff(gr.datetime))
                collect
                filter(_[2] ≤ maxΔt)
                map() do ((φ1, φ2), _)
                    Δcp = mod(φ2 - φ1 + π, 2π) - π
                    U.value(Δcp) / U.uncertainty(Δcp)
                end
            end
            length(s_list) ≥ m.min_cnt || return nothing
            std(s_list)
        end
    end
    length(std_list) ≥ m.min_triangles || return nothing
    return median(std_list)
end


function find_errmul(m::ErrMulSame, uvtbl)
    emuls = @p map(find_errmul(_, uvtbl), m.methods)
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
    NT = intersect_nt_type(eltype(uvtbl), NamedTuple{(:freq_spec, :stokes)})
    Δt = @p let
        uvtbl
        group_vg((;NT(_)..., bl=Baseline(_)))
        flatmap() do gr
            diff(gr.datetime)
        end
        map(_ |> u"minute" |> float)
        median
    end
    1u"s" < Δt < 5u"minute" || @warn "Strange automatic median Δt" Δt
    Δt
end

compute_avgs(uvtbl; maxΔt, min_cnt) = @p let
    average_data(FixedTimeIntervals(maxΔt), uvtbl; avgvals=vals -> (err=U.uncertainty(U.weightedmean(vals)), std=std(U.value.(vals)) / √length(vals) / √2))
    # average_vis(uvtbl, avg_t; avgvals=vals -> (std=(mean(abs2∘value, vals) - abs2(mean(value, vals))) / (length(vals) - 1) / 2 |> sqrt, err=1/√(sum(v -> 1/uncertainty(v)^2, vals)) ))  # formula from difmap
    filter(_.count ≥ min_cnt)
    map(_.value)
end

function compute_diffs(uvtbl; maxΔt)
    NT = intersect_nt_type(eltype(uvtbl), NamedTuple{(:freq_spec, :stokes)})
    @p let
        uvtbl
        group_vg((;NT(_)..., bl=Baseline(_)))
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
end

multiply_errors(errmul, uvtbl) =
    @modify(uvtbl[∗].value |> U.uncertainty) do err
        err * errmul
    end

function rescale_visibility_errors(m, uvtbl)
    mul = find_errmul(m, uvtbl)
    isnothing(mul) && error("Couldn't determine error rescaling factor using $m")
    multiply_errors(mul, uvtbl)
end
