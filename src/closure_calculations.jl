@stable function closures_all(::Type{T}, data::AbstractVector; rescale_redundancy::Bool=false) where {T}
	NT = intersect_nt_type(eltype(data), NamedTuple{(:freq_spec, :stokes, :datetime)})
	@p let
		data
		group_vg(NT)
		flatmap(closures_scan(T, _; rescale_redundancy))
	end
end

@stable function closures_scan(::Type{T}, data::FlexiGroups.GroupArray; rescale_redundancy::Bool=false) where {T}
	bls = @p data map(Baseline(_))
	@assert allunique(bls)
	all_ant_names = @p bls flatmap(antenna_names) unique sort
	closures = @p let
		ant_id_sets_for_closures(T, all_ant_names)
		filtermap() do ant_names
			cbls = (
				(ant_names[1], ant_names[2]),
				(ant_names[2], ant_names[3]),
				(ant_names[3], ant_names[4]),
				(ant_names[4], ant_names[1]),
			)
			curdata = map(cbls) do bl
				rs = @p data filtermap(
					antenna_names(Baseline(_)) == bl ? _ :
					antenna_names(conj(Baseline(_))) == bl ? conjvis(_) :
					nothing)
				isempty(rs) ? nothing : only(rs)
			end
			any(isnothing, curdata) && return nothing

			FlexiGroups.Group(key(data), curdata)
		end
		map(agg_to_closure(T, _))
	end
	if rescale_redundancy && !isempty(closures)
		# n_independent = B - N: gain model g_i + g_j has rank-N design matrix (unsigned incidence matrix)
		n_ind = length(bls) - length(all_ant_names)
		r = length(closures) / n_ind
		closures = map(closures) do c
			@modify(c.value |> U.uncertainty) do err
				err * √r
			end
		end
	end
	FlexiGroups.GroupArray(key(data), closures)
end

@stable function closures_scan(::Type{<:ClosurePhaseSpec}, data::FlexiGroups.GroupArray; rescale_redundancy::Bool=false)
	bls = @p data map(Baseline(_))
	@assert allunique(bls)
	all_ant_names = @p bls flatmap(antenna_names) unique sort
	closures = @p let
		Iterators.product(all_ant_names, all_ant_names, all_ant_names)
		collect
		filter(allunique ⩓ issorted)
		filtermap() do ant_names
			cbls = (
				(ant_names[1], ant_names[2]),
				(ant_names[2], ant_names[3]),
				(ant_names[3], ant_names[1]),
			)
			curdata = map(cbls) do bl
				rs = @p data filtermap(
					antenna_names(Baseline(_)) == bl ? _ :
					antenna_names(conj(Baseline(_))) == bl ? conjvis(_) :
					nothing)
				isempty(rs) ? nothing : only(rs)
			end
			any(isnothing, curdata) && return nothing

			FlexiGroups.Group(key(data), curdata)
		end
		map(agg_to_closure(ClosurePhaseSpec, _))
	end
	if rescale_redundancy && !isempty(closures)
		# n_independent = B - N + 1: gain model g_i - g_j has rank-(N-1) design matrix (signed incidence matrix)
		n_ind = length(bls) - length(all_ant_names) + 1
		r = length(closures) / n_ind
		closures = map(closures) do c
			@modify(c.value |> U.uncertainty) do err
				err * √r
			end
		end
	end
	FlexiGroups.GroupArray(key(data), closures)
end


@stable ant_id_sets_for_closures(::Type{<:ClosurePhaseSpec}, ids) = @p let
	Iterators.product(ids, ids, ids)
	collect
	filter(allunique ⩓ issorted)
end
@stable ant_id_sets_for_closures(::Type{<:ClosureAmpSpec}, ids) = @p let
	Iterators.product(ids, ids, ids, ids)
	collect
	filter(allunique)
	filter(_[1] == minimum(_))
	filter(_[2] < _[4])
end


@stable agg_to_closure(::Type{T}, gr) where {T} = @p gr (;
	key(__)...,
	spec=T(map(x->x.spec, __)),
	value=agg_value_to_closure(T, map(x->x.value, __)),
	values=map(x->x.value, __),
)
agg_value_to_closure(::Type{<:ClosurePhaseSpec}, vals) = vals[1] * vals[2] * vals[3]
agg_value_to_closure(::Type{<:ClosureAmpSpec}, vals) = vals[1] * vals[3] / (vals[2] * vals[4])
