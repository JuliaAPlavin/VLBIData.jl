module VLBISkyModelsExt

using VLBIData
using VLBIData.StructArrays
using VLBIData.Unitful
using VLBIData.AccessorsExtra
using VLBISkyModels.ComradeBase
using VLBISkyModels

ComradeBase.axisdims(m::VLBISkyModels.ModifiedModel) = axisdims(m.model)

function ComradeBase.visibilitymap(model, uvtbl_spec::AbstractVector)
	gim = RectiGrid(axisdims(model))
	guv = let
		uvs = map(UV, uvtbl_spec) |> StructArray
		domain_nt = (U=uvs.u, V=uvs.v)
		if hasproperty(first(uvtbl_spec), :freq_spec)
			domain_nt = @insert domain_nt.Fr = map(r -> ustrip(u"Hz", frequency(r)), uvtbl_spec)
		end
		UnstructuredDomain(domain_nt)
	end
	gfour = FourierDualDomain(gim, guv, NFFTAlg())
	vismap = visibilitymap(model, gfour)
	return uvtable(vismap)
end

end
