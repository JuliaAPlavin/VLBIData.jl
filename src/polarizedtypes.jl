@stable function uvtable_values_to(::Type{CoherencyMatrix}, uvtbl)
	grs = @p uvtbl group_vg((;_.datetime, _.freq_spec, bl=Baseline(_)))
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

@stable function uvtable_values_to(::Type{IPol}, uvtbl)
	grs = @p uvtbl filter(_.stokes == :I || is_parallel_hands(_.stokes)) group_vg((;_.datetime, _.freq_spec, bl=Baseline(_)))
	return map(grs) do gr
		length(gr) ∈ (1, 2) || error("expected 1 or 2 parallel-hand Stokes per group, got $(length(gr))")
		val = mean(x -> x.value, gr)
		@p let
			first(gr)
			@set __.value = val
			@delete __.stokes
		end
	end
end
