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

AccessorsExtra.hasoptic(::Union{ClosurePhaseSpec, ClosureAmpSpec}, ::Type{UV}) = false
@accessor UVs(x::Union{ClosurePhaseSpec,ClosureAmpSpec}) = map(UV, x.vses)

antenna_names(x::Union{ClosurePhaseSpec,ClosureAmpSpec}) = (@assert isvalid(x); map(vs -> first(antenna_names(vs)), x.vses))

function Base.isvalid(x::Union{ClosurePhaseSpec,ClosureAmpSpec})
	for i in eachindex(x.vses)
		antenna_names(x.vses[i])[2] == antenna_names(x.vses[mod1(i+1, length(x.vses))])[1] || return false
	end
	return true
end

visibility(visf::Function, spec::ClosurePhaseSpec) = prod(visf, UVs(spec))

visibility(visf::Function, spec::ClosureAmpSpec) =
	visf(UV(spec.vses[1]))*visf(UV(spec.vses[3])) / (visf(UV(spec.vses[2]))*visf(UV(spec.vses[4])))


UVarea(x) = _UVarea(UVs(x))

# Triangle area - compute 2D determinant directly
function _UVarea(uvs::NTuple{3,StaticVector{2}})
    v1 = uvs[2] - uvs[1]
    v2 = uvs[3] - uvs[1]
    return abs(v1[1] * v2[2] - v1[2] * v2[1]) / 2
end

# Convex quadrangle - try both diagonal splits, take the maximum
function _UVarea(uvs::NTuple{4,StaticVector{2}})
    area1 = _UVarea((uvs[1], uvs[2], uvs[3])) + _UVarea((uvs[1], uvs[3], uvs[4]))
    area2 = _UVarea((uvs[1], uvs[2], uvs[4])) + _UVarea((uvs[2], uvs[3], uvs[4]))
    return max(area1, area2)
end
