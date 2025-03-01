struct Baseline
    array_ix::Int8
    ant_ids::NTuple{2, Int8}
    ant_names::NTuple{2, Symbol}
end
@batteries Baseline selfconstructor=false

function Base.getproperty(b::Baseline, key::Symbol)
    key == :ants_ix && (Base.depwarn("Baseline.ants_ix is deprecated, use ant_ids instead", :ants_ix); return b.ant_ids)
    return getfield(b, key)
end

@deprecate Baseline(array_ix::Integer, ant_ids::NTuple{2, Integer}) Baseline(array_ix, ant_ids, (Symbol(:ANT, ant_ids[1]), Symbol(:ANT, ant_ids[2])))

@accessor Baseline(bl::Baseline) = bl

antennas(b::Baseline) =
    map(b.ant_ids, b.ant_names) do id, name
        Antenna(name, id, SVector(NaN, NaN, NaN))
    end
Accessors.set(b::Baseline, ::typeof(antennas), ants) = setproperties(b, ant_ids=map(a -> a.id, ants), ant_names=map(a -> a.name, ants))

AccessorsExtra.hasoptic(::Baseline, ::Type{UV}) = false
AccessorsExtra.hasoptic(::UV, ::Type{Baseline}) = false

function Base.show(io::IO, bl::Baseline)
	ants = antennas(bl)
	print(io, "Baseline ", (@p ants map(_.id) join(__, "-")), ": ", (@p ants map(_.name) join(__, " - ")))
end
