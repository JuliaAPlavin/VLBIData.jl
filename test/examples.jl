### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 4eb97406-7713-11ec-00ca-8be6bef77030
begin
	# using Revise
	# import Pkg
	# eval(:(Pkg.develop(path="..")))
	# Pkg.resolve()
	using VLBIData
end

# ╔═╡ 5529a3e2-c874-4a21-b349-3984349b1db6
using PyPlotUtils; pyplot_style!()

# ╔═╡ 17a5e5ad-c3a7-41a0-bf56-410675e95475
using VLBIPlots

# ╔═╡ 7a4e96c8-8ab0-4338-82b5-e07f0afdaae5
using DataPipes

# ╔═╡ 9c6efe8c-385e-4446-acdf-bd19cffe31e2
using DisplayAs: Text as AsText

# ╔═╡ be208a3e-9a03-4194-8be7-05547d9871a7
using RectiGrids

# ╔═╡ 4e7278e0-0d12-4762-8ace-e7cbe8eb371c
using StaticArrays

# ╔═╡ dc179b8b-a5b7-4d2d-bf86-afc0fc5a376b
using Unitful, UnitfulAngles

# ╔═╡ f3898d40-517d-421b-b551-00d8ce1f64dd
using AxisKeys

# ╔═╡ 3ba32828-9266-4427-b396-a9c44d81f02b
using HypertextLiteral

# ╔═╡ 5491a70e-8aea-42bc-b449-a56fd120d3fd
md"""
__`VLBIData` reads a range of data formats typically used in very long baseline interferometry (VLBI).__

For convenience, it exports itself as `VLBI`: non-exported functions are accessed as `VLBI.load`.

The proper file format is guessed automatically, and typically one just needs the `VLBI.load(file)` method.

Some formats, such as "image FITS", may actually contain multiple datasets of different kinds (eg image + model). In such cases, `VLBI.load(file)` reads the main dataset (image) by default. Use `VLBI.load(type, file)` to read the dataset of `type` from the `file`: see examples below.
"""

# ╔═╡ 6cb0d385-5b13-4afe-bf8c-015938efc91d
md"""
`VLBIData` is well-tested. It is regularly confirmed to read all files in the [Astrogeo](http://astrogeo.org/vlbi_images/) database successfully.
"""

# ╔═╡ b685e8e4-efa5-453c-a27b-4453946a0839
md"""
# Image FITS
"""

# ╔═╡ 9af33a02-67f8-4233-93f5-4797ecf702d9
md"""
## Read image itself
"""

# ╔═╡ 3df0333d-734c-4d4b-bdec-49cac59f4681
md"""
Load the image data from a FITS file:
"""

# ╔═╡ 6bf4ae14-33d7-4faa-a664-80889ece70b7
fimg = VLBI.load("./data/map.fits");

# ╔═╡ 50c8dfe9-db54-4239-a14d-c8e047eab443
md"""
The returned object contains both the header information, and the image data as an array:
"""

# ╔═╡ 8b722ca5-3198-4c6c-94b1-4daf8fd2e718
fimg.header["OBJECT"]

# ╔═╡ fd52ba81-f515-4388-b677-72ad9bfb4f6b
fimg.data |> AsText

# ╔═╡ 8742d582-8bea-4fcf-a2ea-cf2ae28b8008
@htl("""
<details>
	<summary>Inspect the whole object</summary>
	$(embed_display(fimg))
</details>
""")

# ╔═╡ b104bdf4-cdee-4743-a4e7-31844f5338e7
md"""
As you see, the image data is a `KeyedArray`. It contains information about axes and their units for convenient look up, indexing, plotting, ...:
"""

# ╔═╡ dcd07de0-148d-4641-9cfc-79926e1ba3b1
fimg.data(dec=Near(5u"mas")) |> AsText

# ╔═╡ b8e07154-a609-459d-bc00-6ff61220cd4a
let
	plt.figure()
	imshow_ax(fimg.data, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=5e-4), cmap=:inferno)
	xylabels(fimg.data; inline=true)
	set_xylims((0±30)^2; inv=:x)
	plt.gcf()
end

# ╔═╡ 3075f803-d583-4764-a7d4-1bc14d8e441b
md"""
Access the so-called "image beam" - the effective point spread function:
"""

# ╔═╡ b32c9639-c2c3-4a67-8554-57aefea9d196
beam(fimg)

# ╔═╡ 6d0bf04e-51ee-43b6-8b7c-912fd051dfb5
md"""
## Read source model
"""

# ╔═╡ c4b58f30-e7db-459e-9d26-2d4bc57e5fbd
md"""
FITS image files commonly contain the source model itself, often as a set of delta functions - "CLEAN components". Components are stored in the `AIPS CC` FITS table.

To load this model instead of the image array, pass the corresponding type to `VLBI.load`:
"""

# ╔═╡ b4e14571-a3fa-41d2-98da-cff594f88202
fimg_mod = VLBI.load(MultiComponentModel, "./data/map.fits")

# ╔═╡ c8510200-f820-4b2a-8260-d21e8f7f7880
md"""
The `MultiComponentModel` type, as well as components themselves, are defined in the `InterferometricModels.jl` package. See its help for more details on how to use these models.
"""

# ╔═╡ ebcabc9f-0487-4cad-b7f7-fff448fcf203
md"""
For example, we can convolve the CLEAN model with the image beam:
"""

# ╔═╡ 7994f45f-6e01-4a62-ae3b-3031b46f2d7b
fimg_mod_convolved = convolve(fimg_mod, beam(fimg))

# ╔═╡ 9db795b2-2d79-45b6-a51c-bc09f4d05150
md"""
Convolution turned each delta component into an `EllipticGaussian` component. We can use this model to compute an image: the result should be similar to the image array stored in the FITS file, but without adding residual noise:
"""

# ╔═╡ e06baedf-8689-43e0-aa7d-4cd251879049
fimg_mod_img = intensity(fimg_mod_convolved).(
	grid(SVector,
		ra=range(30, -30, length=200)u"mas",
		dec=range(-30, 30, length=200)u"mas"
	)
);

# ╔═╡ 145c1010-c0e3-4a7a-ace5-cbec71d455f9
fimg_mod_img |> AsText

# ╔═╡ beace2b9-4273-4a8f-94e0-a349f3d6edb8
let
	plt.figure()
	imshow_ax(fimg_mod_img .|> ustrip, ColorBar(unit="Jy/b"); norm=SymLog(linthresh=5e-4), cmap=:inferno)
	xylabels(fimg_mod_img; inline=true)
	plt.gcf()
end

# ╔═╡ 7014ff5c-8975-4219-bae2-b434bbe8441f
md"""
# Visibility or UV FITS
"""

# ╔═╡ 5536507f-12d8-4f49-b163-3df60abea161
md"""
> Note: reading visibilities from UV FITS files currently relies on `astropy`. The Julia `FITSIO.jl` package doesn't support random groups data yet.
"""

# ╔═╡ b0f8bdd5-131d-43a2-845d-5785ca6274f5
md"""
Another common file format in VLBI is "uvfits": it contains interferometric visibilities.

`VLBI.load` can read such files:
"""

# ╔═╡ aad0d6b7-a58e-45f5-83b4-f73a23e5c559
uvfile = VLBI.load("./data/vis.fits");

# ╔═╡ 9ae2c18c-84fe-46e2-b19f-0abb1eccac40
md"""
The `uvfile` object only contains metadata extracted from the FITS header and from auxilary tables:
"""

# ╔═╡ d8f74ef4-fc37-4450-87ec-ef7d55075488
uvfile.ant_arrays

# ╔═╡ d971425a-785e-4424-8630-ded6bba69809
uvfile.freq_windows

# ╔═╡ 442febcb-b7de-4fc4-9a1d-53bb67204092
uvfile.header.stokes

# ╔═╡ 141f6f69-04ce-4bce-b053-3c14791347f8
md"""
... and more. Dump the object to see all fields:
"""

# ╔═╡ aa1d8bd5-e86e-4c56-8fbd-8b8f5dfc7bc8
@htl("""
<details>
	<summary>Inspect the whole object</summary>
	$(embed_display(uvfile))
</details>
""")

# ╔═╡ 546825e2-5f27-47ad-b03a-546a082cd250
md"""
The actual visibility data can be retrieved in several formats. The most useful is likely the table format:
"""

# ╔═╡ ce74a7f5-06be-44c9-9b45-f64dc83e14a8
uvtbl = VLBI.table(uvfile)

# ╔═╡ 81db478c-1185-41a8-ab16-154778c964d4
let
	_, ax = plt.subplots(2, 1, gridspec_kw=Dict(:height_ratios => [3, 1]))
	radplot([abs, angle], uvtbl; ax)
	plt.gcf()
end

# ╔═╡ faf6e8e0-b79d-4e01-94fc-433817c83fd4
md"""
This gives a `Tables.jl`-compatible table. For now it's a `StructArray`, but the concrete type may potentially change.

Fields should be self-explanatory.
"""

# ╔═╡ 77237980-eda9-4e7c-a3ef-25a0425d0149
md"""
# Source model files
"""

# ╔═╡ f946e665-fce7-4b90-a0a0-283ce60eb14b
md"""
Source models, typically consisting of several circular or elliptical Gaussians, can be stored in `.mod` files produced by software such as `difmap`.

They are loaded with `VLBI.load` as well:
"""

# ╔═╡ 5848e406-18cf-4f36-922b-104d6edd2cf5
dmod = VLBI.load("./data/difmap_model.mod")

# ╔═╡ b2a45d8a-c0d6-4e27-8f57-0ff9db75ecc8
md"""
The `MultiComponentModel` type, as well as components themselves, are defined in the `InterferometricModels.jl` package. See its help for more details on how to use these models.
"""

# ╔═╡ bd6255c9-fe17-4a29-8a5c-d9b734ae7c5c
md"""
For now we just plot amplitude and phase envelopes for each component:
"""

# ╔═╡ f2495e21-6d30-404e-810f-f9a26501eb12
let
	_, ax = plt.subplots(2, 1, gridspec_kw=Dict(:height_ratios => [3, 1]))
	for c in components(dmod)
		radplot([abs, angle], c, 0..1e9; ax)
	end
	plt.gcf()
end

# ╔═╡ ebf1c8c4-d5c4-425b-baff-42252c0d52cf


# ╔═╡ 47da1b40-3624-4dfd-a965-d2fe2286d7c1


# ╔═╡ 90b8e37c-1def-438d-981b-e488cc87a3e7


# ╔═╡ 8eb67d05-95f1-47cc-b9b5-1fced799a73a


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
DataPipes = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
DisplayAs = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PyPlotUtils = "5384e752-6c47-47b3-86ac-9d091b110b31"
RectiGrids = "8ac6971d-971d-971d-971d-971d5ab1a71a"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAngles = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
VLBIData = "679fc9cc-3e84-11e9-251b-cbd013bd8115"
VLBIPlots = "0260e397-8112-41bf-b55a-6b4577718f00"

[compat]
AxisKeys = "~0.1.25"
DataPipes = "~0.2.5"
DisplayAs = "~0.1.2"
HypertextLiteral = "~0.9.3"
PyPlotUtils = "~0.1.2"
RectiGrids = "~0.1.6"
Revise = "~3.3.1"
StaticArrays = "~1.3.4"
Unitful = "~1.10.1"
UnitfulAngles = "~0.6.2"
VLBIData = "~0.3.0"
VLBIPlots = "~0.1.0"
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
git-tree-sha1 = "0fa53d25794bf1c754f909ee4b0ac31eabff952f"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.8"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1bdcc02836402d104a46f7843b6e6730b1948264"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "4.0.2"

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
git-tree-sha1 = "8425c47db102577eefb93cb37b4480e750116b0d"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.1"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "9c91a9358de42043c3101e3a29e60883345b0b39"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "4.0.0+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

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
git-tree-sha1 = "e4b24a30eae6ad440ef909ef48d183a7950d6017"
uuid = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
version = "0.2.7"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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
git-tree-sha1 = "0cb6c7a4c30a8185cd2a67fdb0d21301bbebbaec"
uuid = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
version = "0.1.4"

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

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InterferometricModels]]
deps = ["Accessors", "IntervalSets", "LinearAlgebra", "StaticArrays", "Unitful", "UnitfulAstro"]
git-tree-sha1 = "1d0382d02cbc26e4d4dcdb7f91ecd6dd105ab1dc"
uuid = "b395d269-c2ec-4df6-b679-36919ad600ca"
version = "0.1.1"

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
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "b55aae9a2bf436fc797d9c253a900913e0e90178"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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
git-tree-sha1 = "6b0440822974cab904c8b14d79743565140567f6"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

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

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "af6febbfede908c04e19bed954350ac687d892b2"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "0.2.45"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonNegLeastSquares]]
deps = ["Distributed", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "1271344271ffae97e2855b0287356e6ea5c221cc"
uuid = "b7351bd1-99d9-5c5d-8786-f205a815c4d7"
version = "0.4.0"

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

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

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
deps = ["Accessors", "AxisKeys", "Colors", "DataPipes", "DomainSets", "IntervalSets", "LinearAlgebra", "NonNegLeastSquares", "OffsetArrays", "PyCall", "PyPlot", "StatsBase", "Unitful"]
git-tree-sha1 = "803be2153093049a12907fae6a44616b101fcbf6"
uuid = "5384e752-6c47-47b3-86ac-9d091b110b31"
version = "0.1.6"

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
git-tree-sha1 = "d4da8b728580709d736704764e55d6ef38cb7c87"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.5.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "a635a9333989a094bddc9f940c04c549cd66afcf"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.4"

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
git-tree-sha1 = "d6cfdb6ddeb388af1aea38d2b9905fa014d92d98"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "c4e1c470a94063b911fd1b1a204cd2bb34a8cd15"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.1.1"

[[deps.VLBIData]]
deps = ["AxisKeys", "DataPipes", "DateFormats", "Dates", "DelimitedFiles", "FITSIO", "InterferometricModels", "PyCall", "Reexport", "StaticArrays", "StructArrays", "Tables", "Unitful", "UnitfulAngles", "UnitfulAstro"]
path = "../../home/aplavin/.julia/dev/VLBIData"
uuid = "679fc9cc-3e84-11e9-251b-cbd013bd8115"
version = "0.3.5"

[[deps.VLBIPlots]]
deps = ["Colors", "DataPipes", "InterferometricModels", "IntervalSets", "LinearAlgebra", "PyPlotUtils", "Unitful"]
git-tree-sha1 = "f2a51d9775b53ff2378848d7c02ccbbb303f11dd"
uuid = "0260e397-8112-41bf-b55a-6b4577718f00"
version = "0.1.0"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

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
# ╟─5491a70e-8aea-42bc-b449-a56fd120d3fd
# ╟─6cb0d385-5b13-4afe-bf8c-015938efc91d
# ╟─b685e8e4-efa5-453c-a27b-4453946a0839
# ╟─9af33a02-67f8-4233-93f5-4797ecf702d9
# ╟─3df0333d-734c-4d4b-bdec-49cac59f4681
# ╠═6bf4ae14-33d7-4faa-a664-80889ece70b7
# ╟─50c8dfe9-db54-4239-a14d-c8e047eab443
# ╠═8b722ca5-3198-4c6c-94b1-4daf8fd2e718
# ╠═fd52ba81-f515-4388-b677-72ad9bfb4f6b
# ╟─8742d582-8bea-4fcf-a2ea-cf2ae28b8008
# ╟─b104bdf4-cdee-4743-a4e7-31844f5338e7
# ╠═dcd07de0-148d-4641-9cfc-79926e1ba3b1
# ╠═b8e07154-a609-459d-bc00-6ff61220cd4a
# ╟─3075f803-d583-4764-a7d4-1bc14d8e441b
# ╠═b32c9639-c2c3-4a67-8554-57aefea9d196
# ╟─6d0bf04e-51ee-43b6-8b7c-912fd051dfb5
# ╟─c4b58f30-e7db-459e-9d26-2d4bc57e5fbd
# ╠═b4e14571-a3fa-41d2-98da-cff594f88202
# ╟─c8510200-f820-4b2a-8260-d21e8f7f7880
# ╟─ebcabc9f-0487-4cad-b7f7-fff448fcf203
# ╠═7994f45f-6e01-4a62-ae3b-3031b46f2d7b
# ╟─9db795b2-2d79-45b6-a51c-bc09f4d05150
# ╟─145c1010-c0e3-4a7a-ace5-cbec71d455f9
# ╠═e06baedf-8689-43e0-aa7d-4cd251879049
# ╠═beace2b9-4273-4a8f-94e0-a349f3d6edb8
# ╟─7014ff5c-8975-4219-bae2-b434bbe8441f
# ╟─5536507f-12d8-4f49-b163-3df60abea161
# ╟─b0f8bdd5-131d-43a2-845d-5785ca6274f5
# ╠═aad0d6b7-a58e-45f5-83b4-f73a23e5c559
# ╟─9ae2c18c-84fe-46e2-b19f-0abb1eccac40
# ╠═d8f74ef4-fc37-4450-87ec-ef7d55075488
# ╠═d971425a-785e-4424-8630-ded6bba69809
# ╠═442febcb-b7de-4fc4-9a1d-53bb67204092
# ╟─141f6f69-04ce-4bce-b053-3c14791347f8
# ╟─aa1d8bd5-e86e-4c56-8fbd-8b8f5dfc7bc8
# ╟─546825e2-5f27-47ad-b03a-546a082cd250
# ╠═ce74a7f5-06be-44c9-9b45-f64dc83e14a8
# ╠═81db478c-1185-41a8-ab16-154778c964d4
# ╟─faf6e8e0-b79d-4e01-94fc-433817c83fd4
# ╟─77237980-eda9-4e7c-a3ef-25a0425d0149
# ╟─f946e665-fce7-4b90-a0a0-283ce60eb14b
# ╠═5848e406-18cf-4f36-922b-104d6edd2cf5
# ╟─b2a45d8a-c0d6-4e27-8f57-0ff9db75ecc8
# ╟─bd6255c9-fe17-4a29-8a5c-d9b734ae7c5c
# ╠═f2495e21-6d30-404e-810f-f9a26501eb12
# ╠═ebf1c8c4-d5c4-425b-baff-42252c0d52cf
# ╠═47da1b40-3624-4dfd-a965-d2fe2286d7c1
# ╠═90b8e37c-1def-438d-981b-e488cc87a3e7
# ╠═4eb97406-7713-11ec-00ca-8be6bef77030
# ╠═5529a3e2-c874-4a21-b349-3984349b1db6
# ╠═17a5e5ad-c3a7-41a0-bf56-410675e95475
# ╠═7a4e96c8-8ab0-4338-82b5-e07f0afdaae5
# ╠═9c6efe8c-385e-4446-acdf-bd19cffe31e2
# ╠═be208a3e-9a03-4194-8be7-05547d9871a7
# ╠═4e7278e0-0d12-4762-8ace-e7cbe8eb371c
# ╠═dc179b8b-a5b7-4d2d-bf86-afc0fc5a376b
# ╠═f3898d40-517d-421b-b551-00d8ce1f64dd
# ╠═3ba32828-9266-4427-b396-a9c44d81f02b
# ╠═8eb67d05-95f1-47cc-b9b5-1fced799a73a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
