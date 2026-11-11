module VLBIData

using DataManipulation
using AccessorsExtra
using AccessorsExtra: ContainerOptic
using Unitful
using StaticArrays
using StructArrays
using InterferometricModels
import InterferometricModels: UV, visibility, visibilities
using SkyCoords: AbstractSkyCoords, ICRSCoords
using LinearAlgebra: qr, ColumnNorm, dot, norm
using Statistics: mean, median, std
using SpecialFunctions: besselix
using Uncertain
using IntervalSets
using StructHelpers
using EnumX
using Printf
using DateFormats
import PolarizedTypes as PolT
using PolarizedTypes: CoherencyMatrix, IPol

import DispatchDoctor
macro stable(ex)
    return esc(:($(DispatchDoctor).@stable default_mode = "disable" $ex))
end


export Antenna, antenna_names, lonlatalt, Baseline,
	UV, UVs, UVarea, visibility, frequency,
	VisSpec, VisAmpSpec, ClosurePhaseSpec, ClosureAmpSpec,
	uvtable, synthesize_uvtable,
	is_parallel_hands, is_cross_hands, stokes_to_feeds,
	CoherencyMatrix, IPol


# Base has no `filter` on NamedTuples before Julia 1.11
filter(f, xs) = Base.filter(f, xs)
filter(f, xs::NamedTuple) = xs[Base.filter(k -> f(xs[k]), keys(xs))]


include("antenna.jl")
include("baseline.jl")
include("stokes.jl")
include("visspecs.jl")
include("closurespecs.jl")
include("closure_calculations.jl")
include("trivial_calculations.jl")
include("scans.jl")
include("averaging.jl")
include("error_rescaling.jl")
include("selfcal.jl")
include("reindex.jl")
include("polarizedtypes.jl")
include("likelihoods.jl")
include("simulate.jl")


function uvtable end


function frequency end
frequency(q::Unitful.Frequency) = q

function antenna_names end
# fallbacks:
@accessor antenna_names(x) = antenna_names(Baseline(x))
AccessorsExtra.hasoptic(obj, ::typeof(antenna_names)) = AccessorsExtra.hasoptic(obj, Baseline)


@accessor frequency(x::NamedTuple) = frequency(x.freq_spec)
@accessor Baseline(x::NamedTuple) = Baseline(x.spec)
@accessor UV(x::NamedTuple) = UV(x.spec)
@accessor UVs(x::NamedTuple) = UVs(x.spec)
@stable conjvis(x::NamedTuple) = @p let
	x
	modify(conj, __, @o _.spec)
	modify(conj, __, @maybe _.value)
	modify(_reverse_stokes, __, @maybe _.stokes)
end
@accessor visibility(x::NamedTuple) = x.value

# AccessorsExtra's fallback `hasoptic` reports every optic as present; these check the row's `spec`
AccessorsExtra.hasoptic(x::NamedTuple, ::Type{UV}) = hasproperty(x, :spec) && hasoptic(x.spec, UV)
AccessorsExtra.hasoptic(x::NamedTuple, ::Type{UVs}) = hasproperty(x, :spec) && hasoptic(x.spec, UVs)
AccessorsExtra.hasoptic(x::NamedTuple, ::Type{Baseline}) = hasproperty(x, :spec) && hasoptic(x.spec, Baseline)


@generated intersect_nt_type(::Type{<:NamedTuple{KS1}}, ::Type{<:NamedTuple{KS2}}) where {KS1,KS2} = NamedTuple{Tuple(KS1 ∩ KS2)}

"""Property lenses for the keys of `NT`; splat into a `ContainerOptic` so grouping reads only these columns."""
property_lenses(::Type{<:NamedTuple{KS}}) where {KS} = NamedTuple{KS}(map(k -> PropertyLens{k}(), KS))


end
