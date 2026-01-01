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
using StructHelpers


export VLBI,
	Antenna, antennas, Baseline,
	UV, UVs, visibility, frequency,
	VisSpec, VisAmpSpec, ClosurePhaseSpec, ClosureAmpSpec


include("antenna.jl")
include("baseline.jl")
include("visspecs.jl")
include("closurespecs.jl")
include("closure_calculations.jl")
include("trivial_calculations.jl")
include("averaging.jl")
include("error_rescaling.jl")


function frequency end

function antennas end
@accessor antennas(x) = antennas(Baseline(x))
AccessorsExtra.hasoptic(obj, ::typeof(antennas)) = AccessorsExtra.hasoptic(obj, Baseline)


@accessor frequency(x::NamedTuple) = frequency(x.freq_spec)
@accessor Baseline(x::NamedTuple) = Baseline(x.spec)
@accessor UV(x::NamedTuple) = UV(x.spec)
@accessor UVs(x::NamedTuple) = UVs(x.spec)
conjvis(x::NamedTuple) = modify(conj, x, @o _.value _.spec)
@accessor visibility(x::NamedTuple) = x.value

baremodule VLBI
import ..VLBIData:
	Antenna, antennas, Baseline,
	UV, UVs, visibility, frequency,
	VisSpec, VisAmpSpec, ClosurePhaseSpec, ClosureAmpSpec,
	conjvis, add_conjvis, average_bytime, closures_scan, closures_all,
	uvshift,
	ConsecutiveDifferencesStandard, CoherentAverageScatter, ErrMulSame, find_errmul, rescale_visibility_errors
end

end
