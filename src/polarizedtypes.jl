function uvtable_values_to(::Type{CoherencyMatrix}, uvtbl)
	grs = @p uvtbl group_vg((;_.datetime, _.freq_spec, _.spec))
	cnts = @p grs map(length) unique sort
	if cnts != [4]
		error("expected 4 Stokes parameters per group, got $cnts values")
	end

	return map(grs) do gr
		vals = map((:RR, :LR, :RL, :LL)) do stokes
			@p gr filteronly(_.stokes == stokes) __.value
		end
		coher_mat = CoherencyMatrix(vals..., PolT.CirBasis())
		@p let
			first(gr)
			@set __.value = coher_mat
			@delete __.stokes
		end
	end
end

function uvtable_values_to(::Type{IPol}, uvtbl)
	grs = @p uvtbl group_vg((;_.datetime, _.freq_spec, _.spec))
	return map(grs) do gr
		par_hands = @p gr filter(x -> is_parallel_hands(x.stokes))
		length(par_hands) ∈ (1, 2) || error("expected 1 or 2 parallel-hand Stokes per group, got $(length(par_hands))")
		val = mean(x -> x.value, par_hands)
		@p let
			first(gr)
			@set __.value = val
			@delete __.stokes
		end
	end
end
