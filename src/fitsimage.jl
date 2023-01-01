Base.@kwdef struct FitsImage{TD}
    path::String
    header::FITSHeader
    data::TD
end

function load(::Type{FitsImage}, path; read_data=true)
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

        return FitsImage(; path, header, data)
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
