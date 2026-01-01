struct Baseline
    ant_names::NTuple{2, Symbol}
end
@batteries Baseline selfconstructor=false

Baseline(ant_ids::NTuple{2, Integer}) = Baseline(Symbol.(:ANT, ant_ids))

@accessor Baseline(bl::Baseline) = bl

antennas(b::Baseline) = map(Antenna, b.ant_names)
Accessors.set(b::Baseline, ::typeof(antennas), ants) = setproperties(b, ant_names=map(a -> a.name, ants))

AccessorsExtra.hasoptic(::Baseline, ::Type{UV}) = false
AccessorsExtra.hasoptic(::UV, ::Type{Baseline}) = false

function Base.show(io::IO, bl::Baseline)
	ants = antennas(bl)
	print(io, "Baseline ", (@p ants map(_.name) join(__, " - ")))
end
