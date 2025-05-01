abstract type AbstractSpec end

visibility(model, spec::AbstractSpec) = visibility(visibility(model), spec)
visibility(visf::Function, spec::AbstractSpec) = throw(MethodError(visibility, (visf, spec)))

struct VisSpec0{TUV<:UV} <: AbstractSpec
	uv::TUV
end
@batteries VisSpec0

AccessorsExtra.hasoptic(::VisSpec0, ::Type{Baseline}) = false

struct VisSpec{TUV<:UV} <: AbstractSpec
	bl::Baseline
	uv::TUV
end
@batteries VisSpec

struct VisAmpSpec{VS<:VisSpec} <: AbstractSpec
	vs::VS
end
@batteries VisAmpSpec

@accessor VisSpec(vas::VisAmpSpec) = vas.vs
VisAmpSpec(bl::Baseline, uv::UV) = VisAmpSpec(VisSpec(bl, uv))

@accessor Baseline(vs::VisSpec) = vs.bl
@accessor Baseline(vs::VisAmpSpec) = Baseline(vs.vs)

@accessor UV(x::VisSpec0) = x.uv
@accessor UV(x::VisSpec) = x.uv
@accessor UV(x::VisAmpSpec) = UV(x.vs)

Base.conj(spec::VisSpec) = @p let
	spec
	@modify(reverse, antennas(Baseline(__)))
	@modify(-, __.uv)
end
Base.conj(spec::VisAmpSpec) = @modify(conj, spec.vs)

visibility(visf::Function, spec::Union{VisSpec,VisAmpSpec}) = visf(UV(spec))

function Base.show(io::IO, s::AbstractSpec)
	ants = antennas(s)
	print(io,
		chopsuffix(string(typeof(s).name.name), "Spec"), " ",
		(@p ants map(_.id) join(__, "-")),
		": ", (@p ants map(_.name) join(__, " - "))
	)
end

Base.show(io::IO, s::VisSpec0) = print(io, "VisSpec0: ", s.uv)
