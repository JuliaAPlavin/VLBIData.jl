### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 4eb97406-7713-11ec-00ca-8be6bef77030
begin
	using Revise
	import Pkg
	Pkg.develop(path="..")
	Pkg.resolve()
	using VLBIData
end

# ╔═╡ 5529a3e2-c874-4a21-b349-3984349b1db6
begin
	using PyPlotUtils
	pyplot_style!()
end

# ╔═╡ 7a4e96c8-8ab0-4338-82b5-e07f0afdaae5
using DataPipes

# ╔═╡ 9c6efe8c-385e-4446-acdf-bd19cffe31e2
using DisplayAs: Text as AsText

# ╔═╡ 34004188-d3a5-4468-9f23-da7d48f31ec7
using RectiGrids, StaticArrays, Unitful, MyUnitful

# ╔═╡ ee286acc-4797-481e-ab14-e91f91510751
using StatsBase

# ╔═╡ 94e40db8-6aed-4770-b126-7122aca9b495
using Dates

# ╔═╡ 01ea27f2-e96a-4b5a-b44b-f61731c3708c
using IntervalSets

# ╔═╡ 80bb2a26-89fe-4722-829a-e07e8ae6eb0b
using StructArrays

# ╔═╡ 19d92205-59d2-4592-b4dc-a1a1c3c5ea74
using SplitApplyCombine

# ╔═╡ 4d50fe5c-0d73-4ff6-bd19-937fa213e6b7
using LinearAlgebra

# ╔═╡ 70b24dbe-21ff-4375-804d-4cfb808d8f1f
using Altair

# ╔═╡ a203c15a-2473-4a4d-8914-c78ca86fcadb
using AxisKeys

# ╔═╡ 04c814f2-90b3-4ee6-8570-f14cdeeaf2c8
map(readdir("./data", join=true)) do f
	f => VLBI.guess_type(f)
end

# ╔═╡ 6bf4ae14-33d7-4faa-a664-80889ece70b7
VLBI.load("./data/map.fits")

# ╔═╡ b4e14571-a3fa-41d2-98da-cff594f88202
VLBI.load(MultiComponentModel, "./data/map.fits")

# ╔═╡ 0d973937-7a88-4910-8013-0f93c2ac74ff
convolve(VLBI.load(MultiComponentModel, "./data/map.fits"), beam("./data/map.fits"))

# ╔═╡ 406e4220-608c-45f8-8562-3b3970c306dd
VLBI.load(MultiComponentModel, "./data/map.fits")

# ╔═╡ ce74a7f5-06be-44c9-9b45-f64dc83e14a8
@p VLBI.load("./data/vis.fits") |> VLBI.table(stokes=:PI)

# ╔═╡ 5b6f8190-b3d0-403d-92a1-e504db9d6a8f
@p VLBI.load("./data/vis_multichan.vis") |> VLBI.table(stokes=:PI)

# ╔═╡ b89a56a3-b94e-4e34-af11-594d4dcab143
Hour(2)..Hour(4) |> typeof |> supertypes

# ╔═╡ 4854b4e1-131c-4aeb-a06f-c9873e1ea5f6
TimePeriod |> subtypes

# ╔═╡ 75c2edab-75d7-405a-81af-c355b32f03e1
ess(w::AbstractWeights) = sum(w)^2 / sum(w .^ 2)

# ╔═╡ 7a476b7d-713a-4be6-85b7-7a2ee637e8aa
@p VLBI.load("./data/vis_multichan.vis") |>
	VLBI.table()

# ╔═╡ 01ca3ad8-4041-4b99-a7f0-9a3994a4f192
Base.@kwdef struct AverageSpec
	time::ClosedInterval{<:TimePeriod}
	uvspread::ClosedInterval{Float64}
	ess::ClosedInterval{Float64}
end

# ╔═╡ 678f32f6-fe7f-4bcf-bbf1-c2e70b2663c6
function is_too_long(data, rng, spec)
	a, b = data[first(rng)], data[last(rng)]
	b.datetime - a.datetime > maximum(spec.time) && return true
	norm(b.uv - a.uv) > maximum(spec.uvspread) && return true
	# ess(weights(map(:weight, @view(data[rng])))) > maximum(spec.ess) && return true
	return false
end

# ╔═╡ e6004624-93d8-4b25-a2fa-8cfdcbb7972d
function is_too_short(data, rng, spec)
	isempty(rng) && return true
	a, b = data[first(rng)], data[last(rng)]
	b.datetime - a.datetime < minimum(spec.time) && return true
	norm(b.uv - a.uv) < minimum(spec.uvspread) && return true
	# ess(weights(map(:weight, @view(data[rng])))) < minimum(spec.ess) && return true
	return false
end

# ╔═╡ bad8f84e-8d3c-42f8-a8fe-1e806fe5b2a1
begin
	last_ifany(x) = isempty(x) ? nothing : last(x)
	first_ifany(x) = isempty(x) ? nothing : first(x)
end

# ╔═╡ d47f2596-d97f-42e7-9237-4b6f7984b062
function aggregate_obs(X)
	wts = weights(X.weight)
	(
		datetime=unix2datetime(mean(datetime2unix.(X.datetime), wts)),
		timerange=last(X).datetime - first(X).datetime,
		uv=mean(X.uv, wts),
		visibility=mean(X.visibility, wts),
		weight=sum(X.weight, wts),
		ess=ess(wts),
		cnt=length(X),
	)
end

# ╔═╡ b2f030ec-136c-45ee-9068-90885bf0b200
function uvaverage_singlebase(data::AbstractArray, spec::AverageSpec)
	data = sort(data; by=x -> x.datetime)
	res = []
	i_start_a = 1
	while true
		i_end_a = first_ifany(filter(i -> is_too_long(data, i_start_a:i, spec), i_start_a:length(data)))
		i_end_a = something(i_end_a, length(data) + 1) - 1
		if i_end_a == length(data)
			push!(res, aggregate_obs(view(data, i_start_a:i_end_a)))
			i_start_a = i_end_a + 1
		else
			i_end_b = first_ifany(filter(i -> is_too_long(data, (i_end_a + 1):i, spec), (i_end_a + 1):length(data)))
			i_end_b = something(i_end_b, length(data) + 1) - 1
			if is_too_short(data, (i_end_a + 1):i_end_b, spec)
				i_start_b = last_ifany(filter(i -> !is_too_short(data, (i + 1):i_end_b, spec), i_start_a:i_end_a))
				if isnothing(i_start_b)
					push!(res, aggregate_obs(view(data, i_start_a:i_end_a)))
				elseif is_too_short(data, i_start_a:(i_start_b - 1), spec)
					push!(res, aggregate_obs(view(data, i_start_a:i_end_a)))
				else
					push!(res, aggregate_obs(view(data, i_start_a:(i_start_b-1))))
					push!(res, aggregate_obs(view(data, i_start_b:i_end_b)))
				end
			else
				push!(res, aggregate_obs(view(data, i_start_a:i_end_a)))
				push!(res, aggregate_obs(view(data, (i_end_a+1):i_end_b)))
			end
			i_start_a = i_end_b + 1
		end
		i_start_a > length(data) && break
	end
	return res
end

# ╔═╡ 7fd2efb7-fbd8-4cb8-86b1-21f213808a11
function uvaverage(data::AbstractArray, spec::AverageSpec)
	@p begin
		data
		groupview((;_.baseline))
		map() do __
			uvaverage_singlebase(__, spec)
		end
		pairs
		mapmany(_[2], (;_[1]..., _2...))
		map(identity)
	end
end

# ╔═╡ 6a65ce8a-7551-4640-9391-81380e413d8f
@p VLBI.load("./data/vis_multichan.vis") |>
	VLBI.table() |>
	# filter(_.stokes == :LL && _.if_ix == 1) |>
	uvaverage(__, AverageSpec(time=Minute(1)..Minute(4), uvspread=0..Inf, ess=0..Inf))

# ╔═╡ 50668cfa-173e-43d2-bc86-53fac4db43f1
smearing(uvspread, imdist) = sin(π*uvspread*imdist) / (π*uvspread*imdist)

# ╔═╡ fa226e0a-4cd1-43f0-9f11-9ff669e27c44
uvspread(uvdist, mins) = mins * (π*uvdist)/(12*60)

# ╔═╡ a9e71043-e412-45e4-b7ba-2ab5df60887d
maxuvdist(ν) = 6400u"km" / (u"c" / ν)

# ╔═╡ 90ce174d-a85e-404b-8fd5-ef545c04f911
let
	df = @p map(grid(uvspread=[1e3, 1e4, 1e5, 3e5, 5e5, 1e6], imdist=0:50)) do p
		(;p..., smearing=smearing(p.uvspread, p.imdist*u"mas") |> upreferred)
	end |> vec
	ch = altChart(df)
	ch = ch.mark_line().encode(alt.X(:imdist), alt.Y(:smearing, scale=alt.Scale(domain=[0.9, 1]; clamp=true)), alt.Color(:uvspread; type=:nominal))
	altVLSpec(ch)
end

# ╔═╡ ca2a4a92-2d21-4622-ac15-7362d5072d80
let
	df = @p map(grid(uvdist=[1e6, 1e7, 1e8], mins=[0.5, 1, 2, 5], imdist=0:100)) do p
		(;p..., smearing=smearing(uvspread(p.uvdist, p.mins), p.imdist*u"mas") |> upreferred)
		end |> mutate(c="|UV|=$(_.uvdist), $(_.mins) min") |> vec
	ch = altChart(df)
	ch = ch.mark_line().encode(alt.X(:imdist), alt.Y(:smearing, scale=alt.Scale(domain=[0.9, 1]; clamp=true)), alt.StrokeDash(:uvdist; type=:nominal), alt.Color(:mins; type=:nominal))
	altVLSpec(ch)
end

# ╔═╡ e00a2be6-e842-4be6-aeaf-26a7a82a553a
let
	df = @p map(grid(ν=[2, 8, 22], mins=[0.5, 1, 2, 5], imdist=0:100)) do p
		(;p..., smearing=smearing(uvspread(maxuvdist(p.ν*u"GHz")/2, p.mins), p.imdist*u"mas") |> upreferred)
		end |> vec
	ch = altChart(df)
	ch = ch.mark_line().encode(alt.X(:imdist), alt.Y(:smearing, scale=alt.Scale(domain=[0.9, 1]; clamp=true)), alt.StrokeDash(:mins; type=:nominal), alt.Color(:ν; type=:nominal))
	altVLSpec(ch)
end

# ╔═╡ 498bf9cc-ef4b-475f-bab4-abf1287f6f59
smearing(1e6, 50u"mas")|>upreferred

# ╔═╡ 5848e406-18cf-4f36-922b-104d6edd2cf5
VLBI.load("./data/difmap_model.mod")

# ╔═╡ d486df81-8a1b-4785-ba8e-481384b10083
intensity_peak.(components(VLBI.load("./data/difmap_model.mod")))

# ╔═╡ 6d58f33c-c83e-4052-be6a-0e3865b2d106
VLBI.load("./data/difmap_model_empty.mod")

# ╔═╡ beace2b9-4273-4a8f-94e0-a349f3d6edb8
let
	mod = VLBI.load(MultiComponentModel, "./data/map.fits")
	b = beam("./data/map.fits")

	n = 200
	img = intensity(convolve(mod, b)).(grid(SVector, ra=range(30, -30, length=n)u"mas", dec=range(-30, 30, length=n)u"mas")) .|> ustrip
	plt.figure()
	imshow_ax(img, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=5e-4), cmap=:inferno)
	xylabels(img; inline=true)
	plt.gcf()
end

# ╔═╡ b8e07154-a609-459d-bc00-6ff61220cd4a
let
	img = VLBI.load("./data/map.fits").data
	plt.figure()
	imshow_ax(img, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=5e-4), cmap=:inferno)
	xylabels(img; inline=true)
	set_xylims((0±30)^2; inv=:x)
	plt.gcf()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Altair = "b5d8985d-ff0a-46fa-83e6-c6893fdbcf16"
AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
DataPipes = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
DisplayAs = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MyUnitful = "be63a33b-ca4d-43a5-8045-b0b8c6209429"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PyPlotUtils = "5384e752-6c47-47b3-86ac-9d091b110b31"
RectiGrids = "8ac6971d-971d-971d-971d-971d5ab1a71a"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
SplitApplyCombine = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
VLBIData = "679fc9cc-3e84-11e9-251b-cbd013bd8115"

[compat]
Altair = "~0.1.1"
AxisKeys = "~0.1.25"
DataPipes = "~0.2.5"
DisplayAs = "~0.1.2"
IntervalSets = "~0.5.3"
MyUnitful = "~0.1.0"
PyPlotUtils = "~0.1.2"
RectiGrids = "~0.1.6"
Revise = "~3.3.1"
SplitApplyCombine = "~1.2.1"
StaticArrays = "~1.3.2"
StatsBase = "~0.33.14"
StructArrays = "~0.6.4"
Unitful = "~1.10.1"
VLBIData = "~0.3.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.Accessors]]
deps = ["Compat", "CompositionsBase", "ConstructionBase", "Future", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "2e427a6196c7aad4ee35054a9a90e9cb5df5c607"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.7"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.Altair]]
deps = ["JSON", "Pandas", "PyCall", "VegaLite"]
git-tree-sha1 = "46fbac54bc17f3f0a07a808451ab6814e0b49a4d"
uuid = "b5d8985d-ff0a-46fa-83e6-c6893fdbcf16"
version = "0.1.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "ffc6588e17bcfcaa79dfa5b4f417025e755f83fc"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisKeys]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "IntervalSets", "InvertedIndices", "LazyStack", "LinearAlgebra", "NamedDims", "OffsetArrays", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "d6ff375e8229819ef8e443b25009091e0e88eb24"
uuid = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
version = "0.1.25"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "4379a2dac795014534b9895a45889aa658fca213"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.0"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg"]
git-tree-sha1 = "2fabb5fc48d185d104ca7ed7444b475705993447"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "3.49.1+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "6e39c91fb4b84dcb870813c91674bdebb9145895"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.5"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6cdc8832ba11c7695f494c9d9a1c31e90959ce0f"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.6.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CovarianceEstimation]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "a3e070133acab996660d31dcf479ea42849e368f"
uuid = "587fd27a-f159-11e8-2dae-1979310e6154"
version = "0.2.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataPipes]]
deps = ["Accessors", "SplitApplyCombine"]
git-tree-sha1 = "6da546248579d3084c06d1babd1431c64f981441"
uuid = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
version = "0.2.5"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.DateFormats]]
deps = ["Dates", "DocStringExtensions"]
git-tree-sha1 = "c3f0b713b48dab6e052384e6f4d0e077469bf296"
uuid = "44557152-fe0a-4de1-8405-416d90313ce6"
version = "0.1.12"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "66bde31636301f4d217a161cabe42536fa754ec8"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.17"

[[deps.DisplayAs]]
git-tree-sha1 = "44e8d47bc0b56ec09115056a692e5fa0976bfbff"
uuid = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
version = "0.1.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.FITSIO]]
deps = ["CFITSIO", "Printf", "Reexport", "Tables"]
git-tree-sha1 = "e6033823834ec0070125120d4d4a1234f1826a47"
uuid = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
version = "0.16.12"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "67551df041955cc6ee2ed098718c8fcd7fc7aebe"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.12.0"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InterferometricModels]]
deps = ["IntervalSets", "LinearAlgebra", "StaticArrays", "Unitful", "UnitfulAstro"]
git-tree-sha1 = "a33da56b7c4a292fb83638a26ef8b2ef9aa5d065"
uuid = "b395d269-c2ec-4df6-b679-36919ad600ca"
version = "0.1.0"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JSONSchema]]
deps = ["HTTP", "JSON", "URIs"]
git-tree-sha1 = "2f49f7f86762a0fbbeef84912265a1ae61c4ef80"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "0.3.4"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "3cbe45f4871e60fc142154252322bcf9638c2c1d"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyStack]]
deps = ["LinearAlgebra", "NamedDims", "OffsetArrays", "Test", "ZygoteRules"]
git-tree-sha1 = "a8bf67afad3f1ee59d367267adb7c44ccac7fdee"
uuid = "1fad7336-0346-5a1a-a56f-a06ba010965b"
version = "0.0.7"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "f46e8f4e38882b32dcc11c8d31c131d556063f39"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MyUnitful]]
deps = ["Unitful", "UnitfulAstro"]
git-tree-sha1 = "9f983d23f225276098e1828e8dd75719a2aafd3f"
uuid = "be63a33b-ca4d-43a5-8045-b0b8c6209429"
version = "0.1.0"

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "af6febbfede908c04e19bed954350ac687d892b2"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "0.2.45"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NodeJS]]
deps = ["Pkg"]
git-tree-sha1 = "905224bbdd4b555c69bb964514cfa387616f0d3a"
uuid = "2bd173c7-0d6d-553b-b6af-13a54713934c"
version = "1.3.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Pandas]]
deps = ["Compat", "DataValues", "Dates", "IteratorInterfaceExtensions", "Lazy", "OrderedCollections", "Pkg", "PyCall", "Statistics", "TableTraits", "TableTraitsUtils", "Tables"]
git-tree-sha1 = "beefaeb19a644d5166c7b2dff9084ee0e63934a0"
uuid = "eadc2687-ae89-51f9-a5d9-86b5a6373a9c"
version = "1.5.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "71fd4022ecd0c6d20180e23ff1b3e05a143959c2"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.93.0"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

[[deps.PyPlotUtils]]
deps = ["AxisKeys", "DomainSets", "IntervalSets", "OffsetArrays", "PyCall", "PyPlot", "StatsBase", "Unitful"]
git-tree-sha1 = "92b6cb3d131d20ed9c59804ecc3aee2b7a1c03be"
uuid = "5384e752-6c47-47b3-86ac-9d091b110b31"
version = "0.1.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RectiGrids]]
deps = ["AxisKeys", "Random"]
git-tree-sha1 = "f7e83e8dcc2b6e78d09331720d75cc37355d122f"
uuid = "8ac6971d-971d-971d-971d-971d5ab1a71a"
version = "0.1.6"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "2f9d4d6679b5f0394c52731db3794166f49d5131"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "35efd62f6f8d9142052d9c7a84e35cd1f9d2db29"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "b4912cd034cdf968e06ca5f943bb54b17b97793a"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2ae4fe21e97cd13efd857462c1869b73c9f61be3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b95e0b8a8d1b6a6c3e0b3ca393a7a285af47c264"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.10.1"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "dd21b5420bf6e9b76a8c6e56fb575319e7b1f895"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.1"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "c4e1c470a94063b911fd1b1a204cd2bb34a8cd15"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.1.1"

[[deps.Utils]]
deps = ["AxisKeys", "DocStringExtensions", "LinearAlgebra", "Match", "NamedDims", "SplitApplyCombine", "Statistics", "StatsBase"]
git-tree-sha1 = "8d2898067bb18b26d66b3e1e02ded14f6e8ebc34"
uuid = "2e7d9420-a934-49d6-aef3-a4152ec1ce09"
version = "0.4.2"

[[deps.VLBIData]]
deps = ["AxisKeys", "DataPipes", "DateFormats", "Dates", "DelimitedFiles", "FITSIO", "InterferometricModels", "MyUnitful", "PyCall", "Reexport", "SplitApplyCombine", "StaticArrays", "StatsBase", "StructArrays", "Tables", "Unitful", "UnitfulAstro", "Utils"]
path = "../../home/aplavin/.julia/dev/VLBIData"
uuid = "679fc9cc-3e84-11e9-251b-cbd013bd8115"
version = "0.3.1"

[[deps.Vega]]
deps = ["DataStructures", "DataValues", "Dates", "FileIO", "FilePaths", "IteratorInterfaceExtensions", "JSON", "JSONSchema", "MacroTools", "NodeJS", "Pkg", "REPL", "Random", "Setfield", "TableTraits", "TableTraitsUtils", "URIParser"]
git-tree-sha1 = "43f83d3119a868874d18da6bca0f4b5b6aae53f7"
uuid = "239c3e63-733f-47ad-beb7-a12fde22c578"
version = "2.3.0"

[[deps.VegaLite]]
deps = ["Base64", "DataStructures", "DataValues", "Dates", "FileIO", "FilePaths", "IteratorInterfaceExtensions", "JSON", "MacroTools", "NodeJS", "Pkg", "REPL", "Random", "TableTraits", "TableTraitsUtils", "URIParser", "Vega"]
git-tree-sha1 = "3e23f28af36da21bfb4acef08b144f92ad205660"
uuid = "112f6efa-9a02-5b7d-90c0-432ed331239a"
version = "2.6.0"

[[deps.VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═4eb97406-7713-11ec-00ca-8be6bef77030
# ╠═5529a3e2-c874-4a21-b349-3984349b1db6
# ╠═7a4e96c8-8ab0-4338-82b5-e07f0afdaae5
# ╠═9c6efe8c-385e-4446-acdf-bd19cffe31e2
# ╠═34004188-d3a5-4468-9f23-da7d48f31ec7
# ╠═04c814f2-90b3-4ee6-8570-f14cdeeaf2c8
# ╠═6bf4ae14-33d7-4faa-a664-80889ece70b7
# ╠═b4e14571-a3fa-41d2-98da-cff594f88202
# ╠═0d973937-7a88-4910-8013-0f93c2ac74ff
# ╠═406e4220-608c-45f8-8562-3b3970c306dd
# ╠═ce74a7f5-06be-44c9-9b45-f64dc83e14a8
# ╠═ee286acc-4797-481e-ab14-e91f91510751
# ╠═5b6f8190-b3d0-403d-92a1-e504db9d6a8f
# ╠═94e40db8-6aed-4770-b126-7122aca9b495
# ╠═b89a56a3-b94e-4e34-af11-594d4dcab143
# ╠═4854b4e1-131c-4aeb-a06f-c9873e1ea5f6
# ╠═01ea27f2-e96a-4b5a-b44b-f61731c3708c
# ╠═80bb2a26-89fe-4722-829a-e07e8ae6eb0b
# ╠═19d92205-59d2-4592-b4dc-a1a1c3c5ea74
# ╠═4d50fe5c-0d73-4ff6-bd19-937fa213e6b7
# ╠═75c2edab-75d7-405a-81af-c355b32f03e1
# ╠═7a476b7d-713a-4be6-85b7-7a2ee637e8aa
# ╠═6a65ce8a-7551-4640-9391-81380e413d8f
# ╠═01ca3ad8-4041-4b99-a7f0-9a3994a4f192
# ╠═678f32f6-fe7f-4bcf-bbf1-c2e70b2663c6
# ╠═e6004624-93d8-4b25-a2fa-8cfdcbb7972d
# ╠═bad8f84e-8d3c-42f8-a8fe-1e806fe5b2a1
# ╠═d47f2596-d97f-42e7-9237-4b6f7984b062
# ╠═b2f030ec-136c-45ee-9068-90885bf0b200
# ╠═7fd2efb7-fbd8-4cb8-86b1-21f213808a11
# ╠═50668cfa-173e-43d2-bc86-53fac4db43f1
# ╠═fa226e0a-4cd1-43f0-9f11-9ff669e27c44
# ╠═a9e71043-e412-45e4-b7ba-2ab5df60887d
# ╠═70b24dbe-21ff-4375-804d-4cfb808d8f1f
# ╠═90ce174d-a85e-404b-8fd5-ef545c04f911
# ╠═ca2a4a92-2d21-4622-ac15-7362d5072d80
# ╠═e00a2be6-e842-4be6-aeaf-26a7a82a553a
# ╠═498bf9cc-ef4b-475f-bab4-abf1287f6f59
# ╠═a203c15a-2473-4a4d-8914-c78ca86fcadb
# ╠═5848e406-18cf-4f36-922b-104d6edd2cf5
# ╠═d486df81-8a1b-4785-ba8e-481384b10083
# ╠═6d58f33c-c83e-4052-be6a-0e3865b2d106
# ╠═beace2b9-4273-4a8f-94e0-a349f3d6edb8
# ╠═b8e07154-a609-459d-bc00-6ff61220cd4a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
