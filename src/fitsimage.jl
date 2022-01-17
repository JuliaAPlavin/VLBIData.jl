@with_kw struct FitsImage{TD, TCC}
    path::String
    header::FITSHeader
    clean_components::TCC
    data::TD
    noise::Union{Nothing, Float32}
end

function load(::Type{FitsImage}, path; read_data=true, read_clean=false)
    FITS(path) do f
        primary = f[1]
        header = read_header(primary)

        data = if read_data
            mul = 1u"°" |> u"mas"
            ra_vals  = axis_vals(header, "RA---SIN", zero_reference=true) * mul
            dec_vals = axis_vals(header, "DEC--SIN", zero_reference=true) * mul

            @assert header["BSCALE"] == 1 && header["BZERO"] == 0 && header["BUNIT"] == "JY/BEAM"
            data = read(primary)
            data = dropdims(data, dims=(3, 4))
            KeyedArray(data, ra=ra_vals, dec=dec_vals)
        else
            nothing
        end

        comps = if read_clean
            comps = f["AIPS CC"] |> columntable |> rowtable
            if length(Tables.schema(comps).names) > 3
                if !all(c -> c.var"MAJOR AX" == 0 && c.var"MINOR AX" == 0 && c.POSANGLE == 0 && c.var"TYPE OBJ" == 0, comps)
                    @warn "Unexpected component parameters: nonzero major or minor axes, posangle, or type."
                end
            end
            map(comps) do c
                (flux=c.FLUX, radec=(c.DELTAX*3.6e6, c.DELTAY*3.6e6))
            end
        else
            nothing
        end

        return FitsImage(
            path=path,
            header=header,
            clean_components=comps,
            data=data,
            noise=read_data ? mad(data, normalize=true, center=0) : nothing,
        )
    end
end

function pixel_size(fi::FitsImage)
    ra_step  = axis_dict(fi.header, "RA---SIN")["CDELT"] |> abs
    dec_step = axis_dict(fi.header, "DEC--SIN")["CDELT"] |> abs
    @assert ra_step == dec_step
    return ra_step * u"°" .|> u"mas"
end

function pixel_steps(fi::FitsImage)
    ra_step  = axis_dict(fi.header, "RA---SIN")["CDELT"]
    dec_step = axis_dict(fi.header, "DEC--SIN")["CDELT"]
    @assert abs(ra_step) == abs(dec_step)
    return SVector(ra_step, dec_step) .* u"°" .|> u"mas"
end

function pixel_area(fi::FitsImage)
    ra_step  = axis_dict(fi.header, "RA---SIN")["CDELT"] |> abs
    dec_step = axis_dict(fi.header, "DEC--SIN")["CDELT"] |> abs
    return (ra_step * u"°" .|> u"mas") * (dec_step * u"°" .|> u"mas")
end

InterferometricModels.beam(src::AbstractString) = FITS(src) do f
    header = read_header(f[1])
    beam(header)
end
InterferometricModels.beam(fi::FitsImage) = beam(fi.header)
InterferometricModels.beam(fh::FITSHeader) = beam(EllipticGaussian,
    σ_major = InterferometricModels.fwhm_to_σ(fh["BMAJ"] * u"°" .|> u"mas"),
    ratio_minor_major = fh["BMIN"] / fh["BMAJ"],
    pa_major = fh["BPA"] |> deg2rad,
)

frequency(fi::FitsImage) = axis_val(fi.header, "FREQ") * u"Hz" |> u"GHz"
