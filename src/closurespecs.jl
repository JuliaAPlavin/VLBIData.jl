struct ClosurePhaseSpec{VS<:VisSpec} <: AbstractSpec
	vses::NTuple{3,VS}
end
@batteries ClosurePhaseSpec

struct ClosureAmpSpec{VS<:VisSpec} <: AbstractSpec
	vses::NTuple{4,VS}
end
@batteries ClosureAmpSpec

nom_vses(x::ClosureAmpSpec) = (x.vses[1], x.vses[3])
denom_vses(x::ClosureAmpSpec) = (x.vses[2], x.vses[4])

@accessor UVs(x::Union{ClosurePhaseSpec,ClosureAmpSpec}) = map(UV, x.vses)

antennas(x::Union{ClosurePhaseSpec,ClosureAmpSpec}) = (@assert isvalid(x); map(vs -> first(antennas(vs)), x.vses))

function isvalid(x::Union{ClosurePhaseSpec,ClosureAmpSpec})
	allequal(Baseline(vs).array_ix for vs in x.vses) || return false
	for i in eachindex(x.vses)
		antennas(x.vses[i])[2] == antennas(x.vses[mod1(i+1, length(x.vses))])[1] || return false
	end
	return true
end

visibility(visf::Function, spec::ClosureAmpSpec) =
	visf(UV(spec.vses[1]))*visf(UV(spec.vses[3])) / (visf(UV(spec.vses[2]))*visf(UV(spec.vses[4])))
