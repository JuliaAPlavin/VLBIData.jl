struct Baseline
    ant_names::NTuple{2, Symbol}
end
@batteries Baseline selfconstructor=false

Baseline(ant_ids::NTuple{2, Integer}) = Baseline(Symbol.(:ANT, ant_ids))

@accessor Baseline(bl::Baseline) = bl
@accessor antenna_names(bl::Baseline) = bl.ant_names

AccessorsExtra.hasoptic(::Baseline, ::Type{UV}) = false
AccessorsExtra.hasoptic(::UV, ::Type{Baseline}) = false

function Base.show(io::IO, bl::Baseline)
	ants = antenna_names(bl)
	print(io, "Baseline ", join(ants, " - "))
end
