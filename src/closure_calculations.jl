closures_all(::Type{T}, data::AbstractVector) where {T} = @p let
	data
	group_vg((;_.freq_spec, _.stokes, _.datetime))
	flatmap(closures_scan(T, _))
end

function closures_scan(::Type{T}, data::FlexiGroups.GroupArray) where {T}
	bls = @p data map(Baseline)
	@assert allunique(bls)
	array_ix = uniqueonly(map(x -> x.array_ix, bls))
	all_ant_ids = @p bls flatmap(_.ant_ids) unique sort
	@p let
		ant_id_sets_for_closures(T, all_ant_ids)
		filtermap() do ant_ids
			cbls = (
				(ant_ids[1], ant_ids[2]),
				(ant_ids[2], ant_ids[3]),
				(ant_ids[3], ant_ids[4]),
				(ant_ids[4], ant_ids[1]),
			)
			curdata = map(cbls) do bl
				rs = @p data filtermap(
					Baseline(_).ant_ids == bl ? _ :
					conj(_.spec).bl.ant_ids == bl ? conjvis(_) :
					nothing)
				isempty(rs) ? nothing : only(rs)
			end
			any(isnothing, curdata) && return nothing
	
			FlexiGroups.Group(key(data), curdata)
		end
		map(agg_to_closure(T, _))
		FlexiGroups.GroupArray(key(data), __)
	end
end

function closures_scan(::Type{<:ClosurePhaseSpec}, data::FlexiGroups.GroupArray)
	bls = @p data map(Baseline)
	@assert allunique(bls)
	array_ix = uniqueonly(map(x -> x.array_ix, bls))
	all_ant_ids = @p bls flatmap(_.ant_ids) unique sort
	@p let
		Iterators.product(all_ant_ids, all_ant_ids, all_ant_ids)
		collect
		filter(allunique ⩓ issorted)
		filtermap() do ant_ids
			cbls = (
				(ant_ids[1], ant_ids[2]),
				(ant_ids[2], ant_ids[3]),
				(ant_ids[3], ant_ids[1]),
			)
			curdata = map(cbls) do bl
				rs = @p data filtermap(
					Baseline(_).ant_ids == bl ? _ :
					conj(_.spec).bl.ant_ids == bl ? conjvis(_) :
					nothing)
				isempty(rs) ? nothing : only(rs)
			end
			any(isnothing, curdata) && return nothing
	
			FlexiGroups.Group(key(data), curdata)
		end
		map(agg_to_closure(ClosurePhaseSpec, _))
		FlexiGroups.GroupArray(key(data), __)
	end
end


ant_id_sets_for_closures(::Type{<:ClosurePhaseSpec}, ids) = @p let
	Iterators.product(ids, ids, ids)
	collect
	filter(allunique ⩓ issorted)
end
ant_id_sets_for_closures(::Type{<:ClosureAmpSpec}, ids) = @p let
	Iterators.product(ids, ids, ids, ids)
	collect
	filter(allunique)
	filter(_[1] == minimum(_))
	filter(_[2] < _[4])
end


agg_to_closure(::Type{T}, gr) where {T} = @p gr (;
	key(__)...,
	spec=T(map(x->x.spec, __)),
	value=agg_value_to_closure(T, map(x->x.value, __)),
	values=map(x->x.value, __),
)
agg_value_to_closure(::Type{<:ClosurePhaseSpec}, vals) = vals[1] * vals[2] * vals[3]
agg_value_to_closure(::Type{<:ClosureAmpSpec}, vals) = vals[1] * vals[3] / (vals[2] * vals[4])
