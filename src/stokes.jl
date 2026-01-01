_reverse_stokes(s::Symbol) =
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

is_parallel_hands(s::Symbol) =
	s ∈ (:RR, :LL, :XX, :YY) ? true :
	s ∈ (:RL, :LR, :XY, :YX) ? false :
	error("Parallel vs cross hands doesn't make sense for stokes = $s")

is_cross_hands(s::Symbol) =
	s ∈ (:RR, :LL, :XX, :YY) ? false :
	s ∈ (:RL, :LR, :XY, :YX) ? true :
	error("Parallel vs cross hands doesn't make sense for stokes = $s")

stokes_to_feeds(s::Symbol) =
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
