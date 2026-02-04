module VLBIData

using DataManipulation
using AccessorsExtra
using Unitful
using StaticArrays
using StructArrays
using InterferometricModels
import InterferometricModels: UV, visibility
using Statistics: mean, median, quantile, std
using Uncertain
using IntervalSets
using StructHelpers
using EnumX
using Printf
using DateFormats
using UnitfulAngles
import PolarizedTypes as PolT
using PolarizedTypes: CoherencyMatrix, IPol


export VLBI,
	Antenna, antenna_names, Baseline,
	UV, UVs, UVarea, visibility, frequency,
	VisSpec, VisAmpSpec, ClosurePhaseSpec, ClosureAmpSpec,
	uvtable,
	is_parallel_hands, is_cross_hands, stokes_to_feeds,
	CoherencyMatrix, IPol


# before Julia 1.11:
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
include("reindex.jl")
include("polarizedtypes.jl")


function uvtable end


function frequency end
function frequency(q::Quantity)
	q |> u"Hz"  # to throw an error
	q
end

function antenna_names end
# fallbacks:
@accessor antenna_names(x) = antenna_names(Baseline(x))
AccessorsExtra.hasoptic(obj, ::typeof(antenna_names)) = AccessorsExtra.hasoptic(obj, Baseline)


@accessor frequency(x::NamedTuple) = frequency(x.freq_spec)
@accessor Baseline(x::NamedTuple) = Baseline(x.spec)
@accessor UV(x::NamedTuple) = UV(x.spec)
@accessor UVs(x::NamedTuple) = UVs(x.spec)
conjvis(x::NamedTuple) = @p let
	x
	modify(conj, __, @o _.value _.spec)
	modify(_reverse_stokes, __, @maybe _.stokes)
end
@accessor visibility(x::NamedTuple) = x.value

# XXX: hasoptic should be handled by @accessor
AccessorsExtra.hasoptic(x::NamedTuple, ::Type{UV}) = hasproperty(x, :spec) && hasoptic(x.spec, UV)
AccessorsExtra.hasoptic(x::NamedTuple, ::Type{UVs}) = hasproperty(x, :spec) && hasoptic(x.spec, UVs)
AccessorsExtra.hasoptic(x::NamedTuple, ::Type{Baseline}) = hasproperty(x, :spec) && hasoptic(x.spec, Baseline)


@generated intersect_nt_type(::Type{<:NamedTuple{KS1}}, ::Type{<:NamedTuple{KS2}}) where {KS1,KS2} = NamedTuple{Tuple(KS1 ∩ KS2)}


baremodule VLBI
using Reexport
@reexport import ..VLBIData:
	Antenna, antenna_names, Baseline,
	UV, UVs, UVarea, visibility, frequency,
	AbstractSpec, VisSpec, VisSpec0, VisAmpSpec, ClosurePhaseSpec, ClosureAmpSpec,
	conjvis, add_conjvis, average_data, closures_scan, closures_all,
	uvshift,
	ConsecutiveDifferencesStandard, CoherentAverageScatter, ErrMulSame, find_errmul, rescale_visibility_errors,
	uv_reindex, uvtable, uvtable_values_to,
	GapBasedScans, FixedTimeIntervals, ByFrequency, add_scan_ids, scan_intervals,
	is_parallel_hands, is_cross_hands, stokes_to_feeds,
	CoherencyMatrix, IPol
end

end
