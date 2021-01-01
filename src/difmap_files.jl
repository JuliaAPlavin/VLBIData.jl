import CSV
using DataFrames

struct DifmapModel end

function load(::Type{DifmapModel}, source)
    df = CSV.read(source, DataFrame, comment="!", delim=" ", ignorerepeated=true, header=["flux", "radius", "theta_deg", "major", "ratio", "phi_deg", "T"])
    df = mapcols(df) do col
        map(col) do x
            ismissing(x) ? missing : parse(Float64, strip(string(x), 'v'))
        end
    end
    transform!(df,
        [:radius, :theta_deg] => ((r, t) -> r .* sind.(t)) => :ra,
        [:radius, :theta_deg] => ((r, t) -> r .* cosd.(t)) => :dec,
        :theta_deg => ByRow(deg2rad) => :theta,
        :phi_deg => ByRow(deg2rad) => :phi,
    )
    return df
end
