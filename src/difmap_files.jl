function load(T::Type{MultiComponentModel}, src)
    src = abspath(src)  # for RFC.File
    line = first(eachline(src))
    load(T, src, isvalid(line) && startswith(line, r"! \w+") ? Val(:mod) : Val(:fits))
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
    cols = @p begin
        mat
        map(isempty(_) ? nothing : parse(Float64, rstrip(_, 'v')))
        eachcol()
        NamedTuple{(:flux, :radius, :theta_deg, :major, :ratio, :phi_deg, :T, :freq, :specindex)[1:length(__)]}(__)
    end
    if length(cols) == 3
        # clean model: many components of the same type (Point)
        @p begin
            cols
            Tables.rows()
            map() do r
                flux = r.flux * u"Jy"
                coords = r.radius .* SVector(sincosd(r.theta_deg))*u"mas"
                Point(; flux, coords)
            end
            MultiComponentModel(__)
        end
    elseif length(cols) >= 7
        # arbitrary model: typically a few components, can have different types
        @p begin
            cols
            Tables.rows()
            map() do r
                flux = r.flux * u"Jy"
                coords = r.radius .* SVector(sincosd(r.theta_deg))*u"mas"
                σ_major = InterferometricModels.fwhm_to_σ(something(r.major, 0.0))*u"mas"
                if σ_major == 0
                    Point(; flux, coords)
                elseif r.ratio == 1
                    CircularGaussian(; flux, σ=σ_major, coords)
                else
                    EllipticGaussian(; flux, σ_major, ratio_minor_major=r.ratio, pa_major=deg2rad(r.phi_deg), coords)
                end
            end
            MultiComponentModel(Tuple(__))
        end
    else
        error("Unsupported: $(length(cols)) columns")
    end
end

function load(::Type{MultiComponentModel}, src, ::Val{:fits})
    mod = FITS(src) do f
        @p f["AIPS CC"] |>
            columntable |>
            rowtable |>
            map() do c
                ctype = get(c, Symbol("TYPE OBJ"), 0)
                common = (flux = c.FLUX*u"Jy", coords=SVector(c.DELTAX, c.DELTAY) * 3.6e6u"mas")
                if ctype == 0
                    Point(; common...)
                elseif ctype == 1
                    σ_major = InterferometricModels.fwhm_to_σ(c.var"MAJOR AX") * 3.6e6u"mas"
                    if c.var"MAJOR AX" == c.var"MINOR AX"
                        CircularGaussian(; common..., σ=σ_major)
                    else
                        EllipticGaussian(; common..., σ_major, ratio_minor_major=c.var"MINOR AX" / c.var"MAJOR AX", pa_major=deg2rad(c.var"POSANGLE"))
                    end
                else
                    error("Unsupported component type: `TYPE OBJ == $ctype`.")
                end
            end |>
            MultiComponentModel
    end
    if all(c -> c isa Point, components(mod))
        @assert isconcretetype(eltype(components(mod)))
    end
    return mod
end


function save(file, mod::MultiComponentModel)
    open(file, "w") do io
        write(io, "! Flux (Jy) Radius (mas)  Theta (deg)  Major FWHM (mas)  Axial ratio   Phi (deg) T \\\n! Freq (Hz)     SpecIndex\n")
        for comp in components(mod)
            @assert any(T -> comp isa T, (Point, CircularGaussian, EllipticGaussian))
            radius = hypot(coords(comp)...)
            theta° = atand(coords(comp)...)
            write(io, f"""{ustrip(u"Jy", flux(comp)):10.6g}v""")
            write(io, f""" {ustrip(u"mas", radius):11.6g}v""")
            write(io, f""" {theta°:11.6g}v""")
            if !(comp isa Point)
                write(io, f""" {ustrip(u"mas", fwhm_max(comp)):11.6g}v""")
                write(io, f""" {fwhm_min(comp) / fwhm_max(comp):11.6g}{comp isa EllipticGaussian ? "v" : " "}""")
                write(io, f""" {comp isa EllipticGaussian ? rad2deg(position_angle(comp)) : 0:10.6g} """)
                write(io, f""" {comp isa Point ? 0 : 1:d}""")
            end
            write(io, "\n")
        end
    end
end
