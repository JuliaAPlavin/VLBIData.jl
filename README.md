# VLBIData.jl

Simple, composable data model for interferometric visibilities — and common operations on top of it: averaging, error rescaling, closure computation, polarization conversion, ...

Started in 2020 and actively maintained since, VLBIData is one of the earliest and most mature Julia packages for interferometric/VLBI data handling.

Visibility data is represented as an `AbstractVector{NamedTuple}`, so you can use standard Julia tools (filtering, mapping, grouping) alongside VLBIData operations. Designed for radio VLBI, but applicable to any interferometric context.

```julia
using VLBIFiles, VLBIData

uvtbl = load_uvfits("observation.uvfits")        # load with VLBIFiles
uvtbl = VLBI.rescale_visibility_errors(           # fix noise estimates
    VLBI.ConsecutiveDifferencesStandard(), uvtbl)
uvtbl = VLBI.average_data(                        # average to 30s bins
    VLBI.FixedTimeIntervals(30u"s"), uvtbl)
cps = VLBI.closures_all(ClosurePhaseSpec, uvtbl)  # compute closure phases
```

Most functions live in the `VLBI` namespace (e.g. `VLBI.average_data`), while commonly used types are also exported directly (`Antenna`, `Baseline`, `UV`, `VisSpec`, etc.).

## Visibility Table Convention

Everything revolves around the **visibility table** (`uvtable`) — a `Vector{NamedTuple}` (or `StructArray`) where each row is one visibility measurement. All VLBIData operations take and return these tables, so they compose naturally.

A row needs just two fields:

```julia
(
    value = 1.5 + 0.3im ± 0.02,                                # complex visibility (optionally with uncertainty, Uncertain.jl)
    spec  = VisSpec(Baseline((:AA, :AP)), UV(1.2e9, -0.5e9)),  # where this visibility was measured
)
```

Most real data also carries `datetime`, `freq_spec` (frequency/bandwidth/channel), `stokes` (e.g. `:RR`, `:LL`, `:RL`), `source`, and others — VLBIData operations automatically use whichever fields are present.

**[VLBIFiles.jl](https://github.com/JuliaAPlavin/VLBIFiles.jl)** reads standard VLBI file formats directly into this layout.

## Core Types

**`Antenna`** — a station, identified by `name::Symbol`, with position and metadata.

**`Baseline{T}`** — a pair of antennas (or names, or indices): `Baseline((:AA, :AP))`.

**`UV`** — 2D coordinate in the UV plane (`StaticVector`). Also works as an accessor: `UV(row)` pulls UV coordinates out of a table row or spec.

### Visibility Specifications

These describe *what* was measured — the `spec` field in each table row. All are subtypes of `AbstractSpec`:

- **`VisSpec`** — the standard spec: a baseline + UV coordinates.
- **`VisSpec0`** — UV only, when baseline info isn't available.
- **`VisAmpSpec`** — wraps a `VisSpec` for amplitude-only measurements.

### Closure Specifications

- **`ClosurePhaseSpec`** — three `VisSpec`s forming a closure triangle.
- **`ClosureAmpSpec`** — four `VisSpec`s forming a closure quadrangle.

Use `UVs(spec)` to get all UV coordinates, and `UVarea(spec)` for the triangle/quadrangle area in the UV plane.

## Operations

### Closures

Closure quantities are station-based gain-independent — useful when calibration is uncertain. Compute them from any visibility table:

- `VLBI.closures_all(ClosurePhaseSpec, uvtbl)` — closure phases for the whole table
- `VLBI.closures_all(ClosureAmpSpec, uvtbl)` — closure amplitudes for the whole table

Returns a new uvtable with closure specs, values, and metadata carried over from the input.

### Averaging

Reduce data volume while preserving signal. `VLBI.average_data(strategy, uvtbl)` supports three strategies:

- `VLBI.FixedTimeIntervals(interval)` — fixed-width time bins (default: 1 minute)
- `VLBI.GapBasedScans(min_gap)` — average within scans detected by time gaps (default gap: 1 minute)
- `VLBI.ByFrequency()` — average across frequency channels

All perform weighted averaging with proper uncertainty propagation. A `count` field is added to each output row. Pass `avgvals` for custom aggregation.

### Error Rescaling

Weights in uvfits files are often miscalibrated. These methods estimate the true noise level from the data itself:

- `VLBI.rescale_visibility_errors(method, uvtbl)` — returns a table with corrected uncertainties
- `VLBI.find_errmul(method, uvtbl)` — just get the multiplication factor

Available methods:

- `VLBI.ConsecutiveDifferencesStandard(; maxΔt, rayleigh_q=0.5, min_cnt=50)` — from consecutive visibility differences
- `VLBI.CoherentAverageScatter(; maxΔt, min_cnt_avg=10, min_cnt_after=50)` — from scatter within coherent averages
- `VLBI.ErrMulSame(method1, method2, ...; rtol)` — cross-check: require multiple methods to agree within `rtol`

### UV Data Manipulation

- `VLBI.uvshift(center, uvtbl)` — phase-shift visibilities to a new phase center
- `VLBI.add_conjvis(uvtbl)` — append conjugate visibilities (common for Stokes-I imaging and plotting)

## Polarization

Convert between polarization representations:

- `VLBI.uvtable_values_to(CoherencyMatrix, uvtbl)` — group 4 feed correlation rows into coherency matrix–valued rows (for full-polarization work)
- `VLBI.uvtable_values_to(IPol, uvtbl)` — average parallel hands into Stokes I (the common case when you just want total intensity)

Convenience utilities for working with Stokes/feed symbols:

- `is_parallel_hands(s)` / `is_cross_hands(s)` — classify correlation products
- `stokes_to_feeds(s)` — e.g. `stokes_to_feeds(:RL) == (:R, :L)`

# Integration and Interoperability

VLBIData is designed as the data layer in a composable ecosystem:

### VLBIFiles.jl

[VLBIFiles.jl](https://github.com/JuliaAPlavin/VLBIFiles.jl) — read and write VLBI file formats (uvfits, images, models, etc.). UV data loads directly into VLBIData uvtable format.

### VLBIPlots.jl

[VLBIPlots.jl](https://github.com/JuliaAPlavin/VLBIPlots.jl) — plot visibility data, UV coverage, and models. Understands VLBIData tables natively.

### InterferometricModels.jl

[InterferometricModels.jl](https://github.com/JuliaAPlavin/InterferometricModels.jl) — analytic source models (Gaussians, disks, rings, etc.) that can be evaluated on VLBIData specs via `visibility(model, spec)`.

### Comrade.jl

[Comrade.jl](https://github.com/ptiede/Comrade.jl) — Bayesian imaging and modeling framework for VLBI. Data flows both ways:

- `VLBI.uvtable(obs)` — convert a Comrade observation to a VLBIData uvtable
- `Comrade.extract_table(uvtbl; data_products=...)` — go the other direction
- `ComradeBase.visibilitymap(model, specs)` — evaluate Comrade sky models on VLBIData UV coordinates
