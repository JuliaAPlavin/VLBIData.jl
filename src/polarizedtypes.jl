function uvtable_values_to(::Type{CoherencyMatrix}, uvtbl)
	# columnar access below needs a StructArray; normalize once (compile-time no-op when already one)
	tbl = uvtbl isa StructArray ? uvtbl : StructArray(uvtbl)
	# groupview (not group): each cell is a view, no per-row data copy
	@p tbl |>
		groupview_vg((;_.datetime, _.freq_spec, bl=Baseline(_))) |>
		map() do gr
			stk = gr.stokes                       # the cell's stokes, in row order — no temp collections
			vals = gr.value
			slot = feed_slot_map(stk)             # [p,q] → row index of that feed product, or 0
			nan_val = oftype(vals[begin], NaN)
			# value if the product occurs exactly once in the cell, else NaN (duplicates → NaN)
			coher_vals = map(slot) do i
				(i == 0 || count(==(stk[i]), stk) != 1) ? nan_val : vals[i]
			end
			coher_mat = CoherencyMatrix(coher_vals, feed_basis(stk))
			@p let
				first(gr)
				@set __.value = coher_mat
				@delete __.stokes
			end
		end |>
		StructArray()
end

@stable function uvtable_values_to(::Type{IPol}, uvtbl)
	grs = @p uvtbl filter(_.stokes == :I || is_parallel_hands(_.stokes)) group_vg((;_.datetime, _.freq_spec, bl=Baseline(_)))
	return map(grs) do gr
		length(gr) ∈ (1, 2) || error("expected 1 or 2 parallel-hand Stokes per group, got $(length(gr))")
		@p let
			first(gr)
			@set __.value = mean(x -> x.value, gr)
			@set UV(__) = mean(x -> UV(x), gr)
			modify(c -> sum(x -> x.count, gr), __, @maybe _.count)
			@delete __.stokes
		end
	end
end
