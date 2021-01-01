import CSV
using DataFrames

struct DifmapModel end

function load(::Type{DifmapModel}, path)
    df = CSV.read(path, comment="!", delim=" ", ignorerepeated=true, header=["flux", "radius", "theta", "major", "ratio", "phi", "T", "freq", "specindex"])
    df = mapcols(c -> parse.(Float64, strip.(string.(c), 'v')), df)
    select!(df,
        :theta => (x->sin(x |> deg2rad)) => :ra,
        :theta => (x->cos(x |> deg2rad)) => :dec,
    )
    return df
end
