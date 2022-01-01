module VLBIData

using Tables
import FITSIO
using FITSIO: FITSHeader, FITS, TableHDU, read_header, colnames
using PyCall
using Dates
using Unitful, UnitfulAstro
using Parameters: @with_kw
using AxisKeys
import StatsBase: mad
using Utils  # for u"mas"
using StaticArrays
using DelimitedFiles: readdlm


include("fitstables.jl")
include("fitsutils.jl")
include("uvdata.jl")
include("fitsimage.jl")
include("difmap_files.jl")
include("rfc_files.jl")

end
