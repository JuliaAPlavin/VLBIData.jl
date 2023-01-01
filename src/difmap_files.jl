load(T::Type{MultiComponentModel}, src) = if startswith(first(eachline(src)), "! Center")
    load(T, src, Val(:mod))
else
    load(T, src, Val(:fits))
end

function load(::Type{MultiComponentModel}, src, ::Val{:mod})
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
            flux = r.flux * u"Jy"
            coords = r.radius .* SVector(sincosd(r.theta_deg))*u"mas"
            σ_major = InterferometricModels.fwhm_to_σ(r.major)*u"mas"
            if r.major == 0
                Point(; flux, coords)
            elseif r.ratio == 1
                CircularGaussian(; flux, σ=σ_major, coords)
            else
                EllipticGaussian(; flux, σ_major, ratio_minor_major=r.ratio, pa_major=deg2rad(r.phi_deg), coords)
            end
        end
        MultiComponentModel(Tuple(__))
    end
end

function load(::Type{MultiComponentModel}, src, ::Val{:fits})
    FITS(src) do f
        @p f["AIPS CC"] |>
            columntable |>
            rowtable |>
            map() do c
                @assert c.var"TYPE OBJ" == 0
                Point(flux=c.FLUX*u"Jy", coords=SVector(c.DELTAX, c.DELTAY) * 3.6e6u"mas")
            end |>
            @aside(@assert isconcretetype(eltype(__))) |>
            MultiComponentModel
    end
end
