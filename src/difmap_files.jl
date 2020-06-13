import CSV

struct DifmapModel end

function load(::Type{DifmapModel}, path)
    df = CSV.read(path, comment="!", delim=" ", ignorerepeated=true, header=["flux", "radius", "theta", "major", "ratio", "phi", "T", "freq", "specindex"])
    df = mapcols(c -> parse.(Float64, strip.(string.(c), 'v')), df)
    df = df |>
        @mutate(ra = _.radius * sin(_.theta |> deg2rad), dec = _.radius * cos(_.theta |> deg2rad)) |>
        DataFrame
    return df
end
