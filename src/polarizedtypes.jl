@stable function uvtable_values_to(::Type{CoherencyMatrix}, uvtbl)
	grs = @p uvtbl group_vg((;_.datetime, _.freq_spec, bl=Baseline(_)))

	return map(grs) do gr
		nan_val = oftype(first(gr).value, NaN)
		vals = map((:RR, :LR, :RL, :LL)) do stokes
			@oget filteronly(r -> r.stokes == stokes, $gr).value nan_val
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
