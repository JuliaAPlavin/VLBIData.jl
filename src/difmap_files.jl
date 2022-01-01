struct DifmapModel end

function load(::Type{DifmapModel}, source)
    mat = readdlm(source, String, comments=true, comment_char='!')
    coldata = map(eachcol(mat)) do col
        if any(x -> occursin('.', x) || occursin('v', x), col)
            parse.(Float64, rstrip.(col, 'v'))
        else
            parse.(Int, col)
        end
    end
    coltab = NamedTuple{(:flux, :radius, :theta_deg, :major, :ratio, :phi_deg, :T)}(coldata)
    rowtab = map(Tables.rows(coltab)) do r
        (;
            r.flux,
            radec=r.radius .* SVector(sincosd(r.theta_deg)),
            r.major, r.ratio, phi=deg2rad(r.phi_deg),
            r.T
        )
    end
    return rowtab
end
