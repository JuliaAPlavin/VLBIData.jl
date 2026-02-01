module VLBISkyModelsExt

using VLBIData
using VLBIData.StructArrays
using VLBISkyModels.ComradeBase
using VLBISkyModels

ComradeBase.axisdims(m::VLBISkyModels.ModifiedModel) = axisdims(m.model)

function ComradeBase.visibilitymap(model, uvtbl_spec::AbstractVector)
	gim = RectiGrid(axisdims(model))
	guv = let
		uvs = map(UV, uvtbl_spec) |> StructArray
		UnstructuredDomain((U=uvs.u, V=uvs.v)) # Fr=fill(ustrip(u"Hz", _spec.freq), length(G))))
	end
	gfour = FourierDualDomain(gim, guv, NFFTAlg())
	vismap = visibilitymap(model, gfour)
	return uvtable(vismap)
end

end
