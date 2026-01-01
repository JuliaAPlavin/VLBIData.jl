function average_bytime(src_, avg_interval; avgvals=U.weightedmean)
	src = StructArray(src_)
	mindt = minimum(src.datetime)
	@p begin
		src
		groupview_vg((;
			bl=Baseline(_), _.freq_spec, _.stokes,
			timestep = ((_.datetime - mindt) / avg_interval) |> upreferred |> x->trunc(Int, x),
		))
		map((;
			delete(key(_), @o _.bl _.timestep)...,
			count = length(_),
			value = avgvals(_.value),
			datetime = mindt + (key(_).timestep + 0.5) * avg_interval,
			spec = aggspec(key(_).bl, _.spec)
		))
		StructArray()
	end
end

aggspec(bl, specs::AbstractVector{<:VisSpec}) = VisSpec(bl, mean(UV, specs) |> UV)
