# VLBIData.jl

Simple, composable data model for interferometric visibilities — and common operations on top of it: averaging, error rescaling, closure computation, polarization conversion, ...

Started in 2020 and actively maintained since, VLBIData is one of the earliest and most mature Julia packages for interferometric/VLBI data handling.

Visibility data is represented as an `AbstractVector{NamedTuple}`, so you can use standard Julia tools (filtering, mapping, grouping) alongside VLBIData operations. Designed for radio VLBI, but applicable to any interferometric context.

```julia
using VLBIFiles, VLBIData, Unitful

uvtbl = VLBIFiles.load("observation.uvfits") |> uvtable  # load and tabulate
uvtbl = VLBIData.rescale_visibility_errors(              # fix noise estimates
    VLBIData.ConsecutiveDifferencesStandard(), uvtbl)
uvtbl = VLBIData.average_data(                           # average to 30s bins
    VLBIData.FixedTimeIntervals(30u"s"), uvtbl)
cps = VLBIData.closures_all(ClosurePhaseSpec, uvtbl)     # compute closure phases
```

Types and the common accessors (`Antenna`, `Baseline`, `UV`, `VisSpec`, `uvtable`, `visibility`, `frequency`, ...) are exported; the other functions are called package-qualified, e.g. `VLBIData.average_data`.

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

**`Antenna`** — an interferometric endpoint, identified by `name::Symbol`, with position and metadata.

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

Closure quantities are antenna-gain-independent — useful when calibration is uncertain. Compute them from any visibility table:

- `VLBIData.closures_all(ClosurePhaseSpec, uvtbl)` — closure phases for the whole table
- `VLBIData.closures_all(ClosureAmpSpec, uvtbl)` — closure amplitudes for the whole table

Returns a new uvtable with closure specs, values, and metadata carried over from the input.

### Averaging

Reduce data volume while preserving signal. `VLBIData.average_data(strategy, uvtbl)` supports these strategies:

- `VLBIData.FixedTimeIntervals(interval)` — fixed-width time bins
- `VLBIData.GapBasedScans(min_gap)` — average within scans, split at time gaps
- `VLBIData.IntegrationGridScans(min_gap)` — as `GapBasedScans`, also split where the integration grid jumps
- `VLBIData.ByFrequency()` — average across frequency channels

All perform weighted averaging with proper uncertainty propagation; uv coordinates are 1/σ²-weighted when the values carry uncertainties, a plain mean otherwise or without values. A `count` field is added to each output row. Pass `avgvals` for custom aggregation.

### Error Rescaling

Weights in uvfits files are often miscalibrated. These methods estimate the true noise level from the data itself:

- `VLBIData.rescale_visibility_errors(method, uvtbl)` — returns a table with corrected uncertainties
- `VLBIData.find_errmul(method, uvtbl)` — just get the multiplication factor

Methods estimate the factor from consecutive visibility differences (`VLBIData.ConsecutiveDifferencesStandard`), from the scatter within coherent averages, across frequency channels, or from consecutive closure-phase differences, which are robust to antenna gain drifts. `VLBIData.ErrMulSame(method1, method2, ...; rtol)` requires several methods to agree.

### UV Data Manipulation

- `VLBIData.uvshift(center, uvtbl)` — phase-shift visibilities to a new phase center
- `VLBIData.add_conjvis(uvtbl)` — append conjugate visibilities (common for Stokes-I imaging and plotting)

## Polarization

Convert between polarization representations:

- `VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl)` — group 4 feed correlation rows into coherency matrix–valued rows (for full-polarization work)
- `VLBIData.uvtable_values_to(IPol, uvtbl)` — average parallel hands into Stokes I (the common case when you just want total intensity)

Convenience utilities for working with Stokes/feed symbols:

- `is_parallel_hands(s)` / `is_cross_hands(s)` — classify correlation products
- `stokes_to_feeds(s)` — e.g. `stokes_to_feeds(:RL) == (:R, :L)`

# Integration and Interoperability

VLBIData is designed as the data layer in a composable ecosystem:

### VLBIFiles.jl

[VLBIFiles.jl](https://github.com/JuliaAPlavin/VLBIFiles.jl) — read and write VLBI file formats (uvfits, images, models, etc.). UV data loads directly into VLBIData uvtable format.

### VLBIPlots.jl

[VLBIPlots.jl](https://github.com/JuliaAPlavin/VLBIPlots.jl) — plot visibility data, UV coverage, and models. Understands VLBIData tables natively.

<p align="center">
  <img src="https://github.com/JuliaAPlavin/VLBIPlots.jl/raw/master/docs/images/vlbi_data_overview.png" width="44%" />
  <img src="https://github.com/JuliaAPlavin/VLBIPlots.jl/raw/master/docs/images/model_data_comparison.png" width="27%" />
  <img src="https://github.com/JuliaAPlavin/VLBIPlots.jl/raw/master/docs/images/model_visualization.png" width="27%" />
</p>

### InterferometricModels.jl

[InterferometricModels.jl](https://github.com/JuliaAPlavin/InterferometricModels.jl) — analytic source models (Gaussians etc) that can be evaluated on VLBIData specs via `visibility(model, spec)`.

### Comrade.jl

[Comrade.jl](https://github.com/ptiede/Comrade.jl) — Bayesian imaging and modeling framework for VLBI.

- `ComradeBase.visibilitymap(model, specs)` — evaluate a Comrade sky model on VLBIData UV coordinates

```julia
using VLBIFiles, Comrade, SkyCoords

uvf   = VLBIFiles.load("observation.uvfits")
uvtbl = uvtable(uvf)

# Evaluate any Comrade sky model on the uvtable's UV points
vis = visibilitymap(model, uvtbl)  # returns a uvtable with model visibilities
```
