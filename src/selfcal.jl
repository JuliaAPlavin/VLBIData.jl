_model_visibilities(model::AbstractVector, uvtbl) = (@assert length(model) == length(uvtbl); model)
_model_visibilities(model, uvtbl) = [complex(float(visibility(model, row.spec))) for row in uvtbl]

struct AmplitudeAnalytic end

"""
    solve_gains(AmplitudeAnalytic(), model, uvtbl)

Solve for station-based amplitude gains by comparing observed visibilities to a model.

The gain model is `|V_obs(i,j)| = g_i · g_j · |V_model(i,j)|` where `g_i > 0`.
Solved via weighted least squares in log space: `log|V_obs/V_model| = L_i + L_j`.
Returns a StructArray of `(; station, datetime, [stokes, freq_spec, ...], gain)`.

- `model`: anything supporting `visibility(model, spec)`, or a `Vector` of model visibilities
"""
solve_gains(method::AmplitudeAnalytic, model, uvtbl::AbstractVector) = solve_gains(method, model, StructArray(uvtbl))
function solve_gains(::AmplitudeAnalytic, model, uvtbl_::StructArray)
    uvtbl = @insert uvtbl_.model_value = _model_visibilities(model, uvtbl_)
    NT = intersect_nt_type(eltype(uvtbl), NamedTuple{(:datetime, :stokes, :freq_spec)})

    gain_rows = @p let
        uvtbl
        group_vg(NT(_))
        map() do gr
            k = key(gr)
            g = _solve_amp_gains(gr)
            [(; station, k..., gain) for (station, gain) in g]
        end
    end
    return StructArray(reduce(vcat, gain_rows))
end

function _solve_amp_gains(rows)
    stations = unique!(sort!([ant for row in rows for ant in antenna_names(row)]))
    st_idx = Dict(s => i for (i, s) in enumerate(stations))

    # Log-amplitude residuals and weights
    r = map(rows) do r
        V_obs = U.value(r.value)
        V_mod = r.model_value
        log(abs(V_obs / V_mod))
    end
    w = [U.nσ(row.value)^2 for row in rows]
    
    # WLS: A * L = r, where A[k,i]=1, A[k,j]=1 for baseline (i,j)
    # The system L_i + L_j is fully determined (no degeneracy) for a connected array
    N_st = length(stations)
    N_bl = length(r)
    A = zeros(N_bl, N_st)
    for (k, r) in enumerate(rows)
        a1 = st_idx[antenna_names(r)[1]]
        a2 = st_idx[antenna_names(r)[2]]
        A[k, a1] = 1.0
        A[k, a2] = 1.0
    end
    WA = w .* A
    L = (A' * WA) \ (A' * (w .* r))

    return Dict(stations[i] => exp(L[i]) for i in 1:N_st)
end

"""
    apply_gains(AmplitudeAnalytic(), gains, uvtbl; default_gain=NaN)

Apply station-based gains to visibilities.
For `AmplitudeAnalytic`: `V_cal = V_obs / (g_i · g_j)`.
When a gain is missing for a station/time, `default_gain` is used (`NaN` by default).
"""
function apply_gains(::AmplitudeAnalytic, gains, uvtbl; default_gain=NaN)
    NT = intersect_nt_type(eltype(gains), NamedTuple{(:datetime, :stokes, :freq_spec)})
    _gain_key(g) = (; g.station, NT(g)...)
    gdict = Dict(_gain_key(g) => g.gain for g in gains)
    mapset(uvtbl, value=row -> let
        ant1, ant2 = antenna_names(row)
        g1 = get(gdict, (; station=ant1, NT(row)..., datetime=row.datetime), default_gain)
        g2 = get(gdict, (; station=ant2, NT(row)..., datetime=row.datetime), default_gain)
        gg = g1 * conj(g2)
        row.value / gg
    end)
end

"""
    selfcal(AmplitudeAnalytic(), model, uvtbl)

Self-calibrate visibilities against a model. Returns `(; gains, data)`.

```julia
result = selfcal(AmplitudeAnalytic(), model, uvtbl)
result.gains  # StructArray of station gains
result.data   # calibrated visibilities
```
"""
function selfcal(method, model, uvtbl)
    gains = solve_gains(method, model, uvtbl)
    data = apply_gains(method, gains, uvtbl)
    return (; gains, data)
end
