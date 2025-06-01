module ComradeBaseExt

using VLBIData
using VLBIData.StructArrays
using VLBIData.Unitful
using VLBIData.DataManipulation
using ComradeBase


function VLBI.uvtable(usm::ComradeBase.UnstructuredMap)
	hasproperty(usm, :U) && hasproperty(usm, :V) || error("UnstructuredMap must have U and V properties, got $(propertynames(usm))")
	nt = (;spec = UV.(usm.U, usm.V), value = Vector(usm))
	if hasproperty(usm, :Fr)
		nt = (; nt..., freq_spec = usm.Fr .* u"Hz" .|> u"GHz")
	end
	if hasproperty(usm, :Ti)
		nt = (; nt..., datetime = usm.Ti)
	end
	return StructArray(nt)
end

end