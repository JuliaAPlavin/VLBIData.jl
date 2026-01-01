module VLBIData

const VLBI = VLBIData
export VLBI, table, uvtable

using Reexport
using DataManipulation
using Tables
import Tables: table
using FITSIO: FITSHeader, FITS, TableHDU, read_header, colnames, FITSIO
using Dates
using Unitful, UnitfulAstro, UnitfulAngles
using AxisKeys
using StaticArrays
using StructArrays
using AccessorsExtra
using DelimitedFiles: readdlm
using DateFormats: julian_day
using PyFormattedStrings
@reexport using InterferometricModels
import InterferometricModels: UV
using Statistics
using Uncertain


include("grouphdu.jl")

include("fitsutils.jl")
include("uvdata.jl")
include("fitsimage.jl")
include("difmap_files.jl")
include("loading.jl")

using PrecompileTools
@compile_workload begin
    VLBI.load(joinpath(@__DIR__, "../test/data/map.fits"))
    VLBI.load(joinpath(@__DIR__, "../test/data/vis.fits"))
    VLBI.load(joinpath(@__DIR__, "../test/data/difmap_model_empty.mod"))
end

end
