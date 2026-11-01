# uv: 1/σ²-weighted when the values carry uncertainties, a plain mean otherwise or without values
aggspec(specs, vals) = set(first(specs), UV, _avg_uv(specs, vals))

_avg_uv(specs, vals::AbstractVector{<:U.Value}) =
    U.value(U.weightedmean(Iterators.map((s, v) -> U.Value(UV(s), U.uncertainty(v)), specs, vals)))
_avg_uv(specs, ::Union{AbstractVector, Nothing}) = mean(s -> UV(s), specs)

_default_avgvals(vals::AbstractVector{<:U.Value}) = U.weightedmean(vals)
_default_avgvals(vals::AbstractVector) = mean(vals)

_group_values(gr) = hasfield(eltype(gr), :value) ? gr.value : nothing

_avg_value_part(avgvals, vals) = (value = avgvals(vals),)
_avg_value_part(avgvals, ::Nothing) = (;)

"""
    Coherent()

Coherent (complex) visibility averaging mode. This is the default.
Averages complex visibilities with inverse-variance weighting, preserving phase.
"""
struct Coherent end

"""
    Incoherent()

Incoherent (RMS amplitude) visibility averaging mode.
Computes RMS amplitude with phase direction preserved from the coherent average.
Works with both complex (`VisSpec`) and amplitude-only (`VisAmpSpec`) inputs.

To debias amplitudes before averaging, call [`debias_amplitudes`](@ref) first:
```julia
avg_scans = average_data(GapBasedScans(), uvtbl)
avg_freq  = average_data(Incoherent(), ByFrequency(), debias_amplitudes(avg_scans))
```
"""
struct Incoherent end

"""
    debias_amplitudes(uvtbl)

Remove Rice distribution amplitude bias from each visibility (Thompson et al. eq. 9.86).
Returns amplitude-only data with `VisAmpSpec`.

For each visibility with amplitude ``|V|`` and per-component noise ``σ``:
debiased amplitude ``Â = √max(|V|² - 2σ², 0)``.
"""
debias_amplitudes(uvtbl::AbstractVector) = debias_amplitudes(StructArray(uvtbl))
@stable debias_amplitudes(uvtbl::StructArray) = @p let
    uvtbl
    @modify(s -> VisAmpSpec(s), __.spec[∗])
    @modify(__.value[∗]) do value
        amp = abs(U.value(value))
        σ = U.uncertainty(value)
        amp0 = √max(amp^2 - 2σ^2, zero(amp^2))
        sigma0 = _incoh_sigma(amp0, (σ,))
        U.Value(amp0, sigma0)
    end
end

average_data(::Coherent, strategy, uvtbl; kwargs...) = average_data(strategy, uvtbl; kwargs...)

average_data(::Incoherent, strategy, uvtbl) = average_data(strategy, uvtbl; avgvals=_incoh_avg)

@stable function _incoh_avg(vals::AbstractVector{<:Number})
    sigmas = U.uncertainty.(vals)
    amp0 = √mean(v -> abs2(U.value(v)), vals)
    sigma0 = _incoh_sigma(amp0, sigmas)
    # Preserve phase direction from coherent average
    coh_mean = mean(U.value, vals)
    direction = abs(coh_mean) > 0 ? coh_mean / abs(coh_mean) : one(coh_mean)
    return U.Value(amp0 * direction, sigma0)
end

# Error via variance propagation on Â² = mean(|V_i|² - 2σ_i²)
# Var(|V_i|²) = 4σ_i²(A² + σ_i²), so inc_σ = √Var(Â²)/2
# Then delta-method transform from squared-amplitude SNR to amplitude SNR
@stable function _incoh_sigma(amp0, sigmas)
    N = length(sigmas)
    inc_σ = √(amp0^2 * sum(σ -> σ^2, sigmas) + sum(σ -> σ^4, sigmas)) / N
    s = amp0^2 / (2 * inc_σ)  # SNR of squared-amplitude estimator
    if s > 0
        snr_amp = s * (1 + √(1 + 1/s))  # delta-method transform to amplitude SNR
        snr_amp > 0 ? amp0 / snr_amp : √sum(σ -> σ^2, sigmas) / N
    else
        √sum(σ -> σ^2, sigmas) / N
    end
end

# The columns with one value in every row, carried to each averaged row; averaged and grouping columns and `extra` excluded.
_const_part(tbl, extra = ()) =
    @p getproperties(tbl) (@delete __[(:source, :freq_spec, :stokes, :scan_id, :value, :spec, extra...)]) filter(allequal) map(uniqueonly)

function average_data(strategy::AbstractScanStrategy, uvtbl; avgvals=_default_avgvals)
    uvtbl_with_scans = add_scan_ids(strategy, uvtbl)
    return _average_data_by_scans(uvtbl_with_scans, _const_part(uvtbl_with_scans); avgvals)
end

function _average_data_by_scans(uvtbl_with_scans, const_part; avgvals)
    intervals = scan_intervals(uvtbl_with_scans)
    NT = intersect_nt_type(eltype(uvtbl_with_scans), NamedTuple{(:source, :freq_spec, :stokes, :scan_id)})
    @p let
        uvtbl_with_scans
        groupview_vg(ContainerOptic((; bl=@o(Baseline(_.spec)), property_lenses(NT)...)))
        map((;
            const_part...,
            delete(key(_), @o _.bl)...,
            count = length(_),
            _avg_value_part(avgvals, _group_values(_))...,
            datetime = _mean(intervals[key(_).scan_id]),
            spec = aggspec(_.spec, _group_values(_)),
        ))
        StructArray()
    end
end

@kwdef struct FixedTimeIntervals{T}
    interval::T = 1u"minute"
end

function average_data(strategy::FixedTimeIntervals, src_; avgvals=_default_avgvals)
    src = StructArray(src_)
    return _average_data_by_time(strategy, src, _const_part(src); avgvals)
end

function _average_data_by_time(strategy::FixedTimeIntervals, src, const_part; avgvals)
    mindt = minimum(src.datetime)
    avg_interval = @p let
        float(strategy.interval)
        iszero(__) ? oftype(__, 1u"ms") : __
    end
    NT = intersect_nt_type(eltype(src), NamedTuple{(:source, :freq_spec, :stokes, :scan_id)})
    @p begin
        src
        groupview_vg(ContainerOptic((;
            bl=@o(Baseline(_.spec)), property_lenses(NT)...,
            timestep=@o(trunc(Int, upreferred((_.datetime - mindt) / avg_interval))),
        )))
        map((;
            const_part...,
            delete(key(_), @o _.bl _.timestep)...,
            count = length(_),
            _avg_value_part(avgvals, _group_values(_))...,
            datetime = mindt + (key(_).timestep + 0.5) * avg_interval,
            spec = aggspec(_.spec, _group_values(_)),
        ))
        StructArray()
    end
end


struct ByFrequency end

_aggfreq(freq_specs::AbstractVector{<:Number}) = mean(freq_specs)

function average_data(::ByFrequency, uvtbl; avgvals=_default_avgvals)
    uvtbl = StructArray(uvtbl)
    merged_fs = _aggfreq(unique(uvtbl.freq_spec))
    _average_data_by_frequency(uvtbl, _const_part(uvtbl, (:count, :datetime)), merged_fs; avgvals)
end

function _average_data_by_frequency(uvtbl, const_part, merged_fs; avgvals)
    NT = intersect_nt_type(eltype(uvtbl), NamedTuple{(:source, :stokes, :scan_id, :datetime)})
    @p begin
        uvtbl
        groupview_vg(ContainerOptic((; bl=@o(Baseline(_.spec)), property_lenses(NT)...)))
        map((;
            const_part...,
            delete(key(_), @o _.bl)...,
            freq_spec=merged_fs,
            count=length(_),
            _avg_value_part(avgvals, _group_values(_))...,
            spec=aggspec(_.spec, _group_values(_)),
        ))
        StructArray()
    end
end


_mean(i::Interval) = minimum(i) + (maximum(i) - minimum(i)) ÷ 2
