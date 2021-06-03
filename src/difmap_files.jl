struct DifmapModel end

function load(::Type{DifmapModel}, source; radec_type::Type{RADEC_TYPE}=SVector) where {RADEC_TYPE}
    mat = try
        readdlm(source, String, comments=true, comment_char='!')::Matrix{String}
    catch err
        if hasfield(typeof(err), :msg) && occursin("number of rows in dims must be > 0, got 0", err.msg)
            fill("", (0, 7))  # empty string matrix with 7 columns
        else
            rethrow()
        end
    end
    coldata = map(mat) do x
        parse(Float64, rstrip(x, 'v'))
    end
    coltab = NamedTuple{(:flux, :radius, :theta_deg, :major, :ratio, :phi_deg, :T), NTuple{7, Vector{Float64}}}(eachcol(coldata))
    rowtab = map(Tables.rows(coltab)) do r
        (;
            r.flux,
            radec=r.radius .* RADEC_TYPE(sincosd(r.theta_deg)),
            r.major,
            r.ratio,
            phi=deg2rad(r.phi_deg),
            T=convert(Int, r.T),
        )
    end
    return rowtab
end
