abstract type AbstractScanStrategy end

"""
    GapBasedScans(; min_gap=1u"minute")

Scans as maximal sets of rows that share one source and are contiguous in time for every antenna
taking part.

Each antenna's own rows are sorted by time and cut where the gap exceeds `min_gap` or the source
changes. Per antenna that cut is exact — an antenna points at one source at a time — and it is the
only place the source enters. Rows are then joined into scans through the antennas they share, so a
scan is a connected set of baselines observing one source over one stretch of time.

A schedule may SUBNET: disjoint sub-arrays observing different sources at the same time. Those are
different scans whose time spans OVERLAP, so a scan is not an interval of the session's time axis and
[`scan_intervals`](@ref) may return overlapping intervals. On a schedule that observes one source at
a time this gives the same runs, numbered the same way, as cutting the globally time-sorted rows at
gaps and source changes.

Ids are numbered by `(first time, first antenna, source)` — properties of the data, so they never
depend on row order.

A table that names no antennas (neither a `baseline` column nor a `spec` one) has ONE stream and the
rule degenerates to cutting the whole time-sorted table. That is the same definition with a single
antenna set, and the one case where the difference shows — two sources at one instant — fails loud
rather than being cut by row order.
"""
@kwdef struct GapBasedScans <: AbstractScanStrategy
    min_gap = 1u"minute"
end

"""
    IntegrationGridScans(; min_gap=1u"minute")

[`GapBasedScans`](@ref) that also ends an antenna's run where its next row leaves the integration grid the
run started on: a different `int_time`, or more than [`GRID_TOLERANCE`](@ref) steps off the grid of step
`int_time` anchored at the run's first row. DiFX aligns integrations to each correlator scan's start, so
such a jump marks a correlator scan boundary. Needs an `int_time` column.
"""
@kwdef struct IntegrationGridScans <: AbstractScanStrategy
    min_gap = 1u"minute"
end

"Largest offset from an integration grid, in steps, of a sample still on it."
const GRID_TOLERANCE = 1e-2

"Position of `t`, in steps of `Δ`, on the grid anchored at `t0`: exact, integer milliseconds over the step."
grid_position(t, t0, Δ) = ustrip(u"ms", t - t0) // (1000 * Rational(ustrip(u"s", Δ)))

"Whether grid position `x` lies within [`GRID_TOLERANCE`](@ref) of a grid point."
on_grid(x::Rational) = abs(x - round(x)) ≤ GRID_TOLERANCE

# Each row's antennas. A uv table names its baseline directly (`baseline`, the wide per-record shape)
# or inside `spec` (the narrow per-datum shape); `nothing` when it names none.
function _row_antennas(uvtbl::StructVector)
    hasproperty(uvtbl, :baseline) && return map(antenna_names, uvtbl.baseline)
    hasproperty(uvtbl, :spec) && return map(r -> antenna_names(Baseline(r)), uvtbl)
    nothing
end

# The row streams the gap rule is applied within: one per antenna, or the whole table when it names
# none. Named, so that a violation of the one-source-at-a-time invariant can say where it is.
_row_streams(::Nothing, n) = [(Symbol("") => collect(1:n))]
function _row_streams(ants::AbstractVector, n)
    out = Dict{Symbol, Vector{Int}}()
    for i in 1:n, a in ants[i]
        push!(get!(Vector{Int}, out, a), i)
    end
    sort!(collect(pairs(out)); by = first)
end

_antenna_key(::Nothing, i) = Symbol("")
_antenna_key(ants::AbstractVector, i) = minimum(ants[i])
_source_of(::Nothing, i) = Symbol("")
_source_of(srcs::AbstractVector, i) = Symbol(srcs[i])
_source_changed(::Nothing, i, j) = false
_source_changed(srcs::AbstractVector, i, j) = srcs[i] != srcs[j]

function _find!(parent::Vector{Int}, x::Int)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    x
end

function _union!(parent::Vector{Int}, a::Int, b::Int)
    ra = _find!(parent, a); rb = _find!(parent, b)
    ra == rb && return nothing
    ra < rb ? (parent[rb] = ra) : (parent[ra] = rb)
    nothing
end

scan_ids(strategy::Union{GapBasedScans, IntegrationGridScans}, uvtbl::StructVector) = _scan_ids(strategy, uvtbl)

# Whether row `i` leaves the integration grid of the run anchored at row `a`.
_leaves_grid(::GapBasedScans, uvtbl, i, a) = false
_leaves_grid(::IntegrationGridScans, uvtbl, i, a) =
    uvtbl.int_time[i] != uvtbl.int_time[a] ||
    !on_grid(grid_position(uvtbl.datetime[i], uvtbl.datetime[a], uvtbl.int_time[a]))

# The per-antenna run loop shared by the scan strategies.
@stable function _scan_ids(strategy, uvtbl::StructVector)
    min_gap::typeof(1.0u"s") = strategy.min_gap
    times = uvtbl.datetime
    ants = _row_antennas(uvtbl)
    srcs = hasproperty(uvtbl, :source) ? uvtbl.source : nothing
    parent = collect(1:length(times))
    for (ant, ix) in _row_streams(ants, length(times))
        # `(time, index)` and not time alone: only rows at the same instant can tie, they join the
        # same run whichever way the tie falls, but the order still has to be fixed.
        sort!(ix; by = i -> (times[i], i))
        a = first(ix)   # the current run's first row
        for k in 2:length(ix)
            i = ix[k]; j = ix[k - 1]
            (times[i] == times[j] && _source_changed(srcs, i, j)) &&
                error("$(nameof(typeof(strategy))): $(ant == Symbol("") ? "the table" : "antenna $ant") carries " *
                      "sources $(_source_of(srcs, j)) and $(_source_of(srcs, i)) at the same instant " *
                      "$(times[i]) — one antenna observes one source at a time, so a scan is not " *
                      "defined here; a subnetting schedule needs the table to name its antennas")
            if times[i] - times[j] > min_gap || _source_changed(srcs, i, j) || _leaves_grid(strategy, uvtbl, i, a)
                a = i
                continue
            end
            _union!(parent, j, i)
        end
    end
    _numbered_by_content(parent, times, ants, srcs)
end

# Scans numbered by `(first time, first antenna, source)`: content only, so the numbering is the same
# under any permutation of the rows.
function _numbered_by_content(parent::Vector{Int}, times, ants, srcs)
    roots = [_find!(parent, i) for i in eachindex(parent)]
    firsts = Dict{Int, Tuple{eltype(times), Symbol, Symbol}}()
    for i in eachindex(roots)
        k = (times[i], _antenna_key(ants, i), _source_of(srcs, i))
        firsts[roots[i]] = haskey(firsts, roots[i]) ? min(firsts[roots[i]], k) : k
    end
    order = sort!(collect(keys(firsts)); by = r -> firsts[r])
    remap = Dict(r => n for (n, r) in enumerate(order))
    [remap[r] for r in roots]
end

function add_scan_ids(strategy::AbstractScanStrategy, uvtbl)
    uvtbl = convert(StructArray, uvtbl)  # no-op if already a StructArray — avoids copying (materializing) lazy columns
    @assert !hasproperty(uvtbl, :scan_id)
    @insert uvtbl.scan_id = scan_ids(strategy, uvtbl)
end

"""
    scan_intervals([strategy], uvtbl) -> Vector{Interval}

Each scan's `[first, last]` time, indexed by scan id. Concurrent scans on different sources give
OVERLAPPING intervals — see [`GapBasedScans`](@ref).
"""
@stable scan_intervals(uvtbl) = scan_intervals(nothing, uvtbl)
scan_intervals(strategy, uvtbl) = scan_intervals(strategy, convert(StructArray, uvtbl))
@stable function scan_intervals(strategy::Nothing, uvtbl::StructArray)
    hasproperty(uvtbl, :scan_id) || error("No scan strategy provided and no scan_id present in the table.")
    @p let
        uvtbl
        groupview(_.scan_id; into=Vector)
        map() do gr
            Interval(extrema(gr.datetime)...)
        end
    end
end
@stable function scan_intervals(strategy::AbstractScanStrategy, uvtbl::StructArray)
    uvtbl = add_scan_ids(strategy, uvtbl)
    @p let
        uvtbl
        groupview(_.scan_id; into=Vector)
        map() do gr
            Interval(extrema(gr.datetime)...)
        end
    end
end
