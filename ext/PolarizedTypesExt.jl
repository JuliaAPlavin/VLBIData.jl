module PolarizedTypesExt

using VLBIData
using VLBIData.DataManipulation
import PolarizedTypes as PolT

function VLBIData.uvtable_values_to(::Type{PolT.CoherencyMatrix}, uvtbl)
	grs = @p uvtbl group_vg((;_.datetime, _.freq_spec, _.spec))
	cnts = @p grs map(length) unique sort
	if cnts != [4]
		error("expected 4 Stokes parameters per group, got $cnts values")
	end

	return map(grs) do gr
		vals = map((:RR, :LR, :RL, :LL)) do stokes
			@p gr filteronly(_.stokes == stokes) __.value
		end
		coher_mat = PolT.CoherencyMatrix(vals..., PolT.CirBasis())
		@p let
			first(gr)
			@set __.value = coher_mat
			@delete __.stokes
		end
	end
end

end
