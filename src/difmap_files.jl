struct DifmapModel end

# don't use: keeping for backwards compatibility
function load(::Type{DifmapModel}, src)
    @p begin
        load(MultiComponentModel, src)
        components()
        map() do c
            (;
                flux=flux(c),
                radec=coords(c),
                major=fwhm_max(c),
                ratio=fwhm_min(c) / fwhm_max(c),
                phi=applicable(position_angle, c) ? position_angle(c) : 0.0,
            )
        end
        collect()
    end
end

function load(::Type{MultiComponentModel}, src)
    mat = try
        readdlm(src, String, comments=true, comment_char='!')::Matrix{String}
    catch err
        if hasfield(typeof(err), :msg) && occursin("number of rows in dims must be > 0, got 0", err.msg)
            fill("", (0, 7))  # empty string matrix with 7 columns
        else
            rethrow()
        end
    end
    @p begin
        mat
        map(parse(Float64, rstrip(_, 'v')))
        eachcol()
        NamedTuple{(:flux, :radius, :theta_deg, :major, :ratio, :phi_deg, :T), NTuple{7, Vector{Float64}}}()
        Tables.rows()
        map() do r
            coords = r.radius .* SVector(sincosd(r.theta_deg))
            σ_major = InterferometricModels.fwhm_to_σ(r.major)
            if r.major == 0
                Point(; r.flux, coords)
            elseif r.ratio == 1
                CircularGaussian(; r.flux, σ=σ_major, coords)
            else
                EllipticGaussian(; r.flux, σ_major, ratio_minor_major=r.ratio, pa_major=deg2rad(r.phi_deg), coords)
            end
        end
        MultiComponentModel(Tuple(__))
    end
end
