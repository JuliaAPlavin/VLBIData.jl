import CSV
using DataFrames

struct DifmapModel end

function load(::Type{DifmapModel}, source)
    df = CSV.read(source, DataFrame, comment="!", delim=" ", ignorerepeated=true, header=["flux", "radius", "theta", "major", "ratio", "phi", "T"])
    df = mapcols(df) do col
        map(col) do x
            ismissing(x) ? missing : parse(Float64, strip(string(x), 'v'))
        end
    end
    transform!(df,
        [:radius, :theta] => ((r, t) -> r .* sind.(t)) => :ra,
        [:radius, :theta] => ((r, t) -> r .* cosd.(t)) => :dec,
    )
    return df
end
