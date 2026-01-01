struct Baseline{T}
    antennas::NTuple{2, T}
end
@batteries Baseline selfconstructor=false

@accessor Baseline(bl::Baseline) = bl
@accessor antenna_names(bl::Baseline{Symbol}) = bl.antennas
@accessor antenna_names(bl::Baseline{Int}) = Symbol.(:ANT, bl.antennas)
@accessor antenna_names(bl::Baseline{Antenna}) = map(a -> a.name, bl.antennas)
antenna_names(bl::Baseline) = error("Not implemented for $bl")

Base.conj(bl::Baseline) = @modify(reverse, bl.antennas)

AccessorsExtra.hasoptic(::Baseline, ::Type{UV}) = false
AccessorsExtra.hasoptic(::UV, ::Type{Baseline}) = false

function Base.show(io::IO, bl::Baseline)
	ants = antenna_names(bl)
	print(io, "Baseline ", join(ants, " - "))
end
