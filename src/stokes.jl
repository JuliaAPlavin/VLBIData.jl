@stable _reverse_stokes(s::Symbol) =
	s ∈ (:RR, :LL, :XX, :YY, :I, :total) ? s :
	s == :RL ? :LR :
	s == :LR ? :RL :
	s == :XY ? :YX :
	s == :YX ? :XY :
	s == :RX ? :XR :
	s == :XR ? :RX :
	s == :RY ? :YR :
	s == :YR ? :RY :
	s == :LX ? :XL :
	s == :XL ? :LX :
	s == :LY ? :YL :
	s == :YL ? :LY :
	error("Cannot reverse stokes $s")

@stable is_parallel_hands(s::Symbol) =
	s ∈ (:RR, :LL, :XX, :YY) ? true :
	s ∈ (:RL, :LR, :XY, :YX) ? false :
	error("Parallel vs cross hands doesn't make sense for stokes = $s")

@stable is_cross_hands(s::Symbol) =
	s ∈ (:RR, :LL, :XX, :YY) ? false :
	s ∈ (:RL, :LR, :XY, :YX) ? true :
	error("Parallel vs cross hands doesn't make sense for stokes = $s")

@stable stokes_to_feeds(s::Symbol) =
	s == :RR ? (:R, :R) :
	s == :LL ? (:L, :L) :
	s == :XX ? (:X, :X) :
	s == :YY ? (:Y, :Y) :
	s == :RL ? (:R, :L) :
	s == :LR ? (:L, :R) :
	s == :XY ? (:X, :Y) :
	s == :YX ? (:Y, :X) :
	s == :RX ? (:R, :X) :
	s == :XR ? (:X, :R) :
	s == :RY ? (:R, :Y) :
	s == :YR ? (:Y, :R) :
	s == :LX ? (:L, :X) :
	s == :XL ? (:X, :L) :
	s == :LY ? (:L, :Y) :
	s == :YL ? (:Y, :L) :
	error("Cannot convert stokes $s to feed types")

@stable feeds_to_stokes(feeds::NTuple{2,Any}) = Symbol(feeds...)
Accessors.InverseFunctions.inverse(::typeof(stokes_to_feeds)) = feeds_to_stokes
Accessors.InverseFunctions.inverse(::typeof(feeds_to_stokes)) = stokes_to_feeds

# The two feeds of a stokes axis, R/X first; single-pol axes (e.g. [:RR]) resolve too.
@stable function _feed_pair(stokes::AbstractVector{Symbol})
	circ = any(s -> s ∈ (:RR, :LL, :RL, :LR), stokes)
	lin  = any(s -> s ∈ (:XX, :YY, :XY, :YX), stokes)
	circ && !lin ? (:R, :L) :
	lin && !circ ? (:X, :Y) :
	error("stokes $stokes are not a single circular/linear basis")
end

# The four feed products (p×q) of a basis in coherency-matrix layout [RR RL; LR LL]. Precomputed
# constants — building them per call via `Symbol(f, f)` interns a fresh String every time.
const _PRODUCTS_CIR = SMatrix{2,2,Symbol}(:RR, :LR, :RL, :LL)
const _PRODUCTS_LIN = SMatrix{2,2,Symbol}(:XX, :YX, :XY, :YY)
_feed_products(stokes) = _feed_pair(stokes) == (:R, :L) ? _PRODUCTS_CIR : _PRODUCTS_LIN

# [p,q]: the index in `stokes` of antenna-1 feed p × antenna-2 feed q, or 0 if absent; slot 1 is R/X and
# slot 2 L/Y, as in a CoherencyMatrix [RR RL; LR LL].
@stable function feed_slot_map(stokes::AbstractVector{Symbol})
	targets = _feed_products(stokes)
	SMatrix{2,2,Int8}(map(t -> something(findfirst(==(t), stokes), 0), targets))
end

feed_basis(stokes::AbstractVector{Symbol}) = _feed_pair(stokes) == (:R,:L) ? PolT.CirBasis() : PolT.LinBasis()
