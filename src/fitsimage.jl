Base.@kwdef struct FitsImage{TD}
    path::String
    header::FITSHeader
    data::TD
end

function load(::Type{FitsImage}, path; read_data=true)
    path = abspath(path)  # for RFC.File
    FITS(path) do f
        primary = f[1]
        header = read_header(primary)

        data = if read_data
            haskey(header, "BSCALE") && @assert header["BSCALE"] == 1
            haskey(header, "BZERO" ) && @assert header["BZERO" ] == 0
            header["BUNIT"] ∈ ("JY/BEAM", "JY/PIXEL") || error("Unknown BUNIT: $(header["BUNIT"])")
            data = read(primary)
            ndims(data) ∈ (2, 4) || error("Unexpected data shape: $(size(data))")
            if ndims(data) == 4
                data = dropdims(data, dims=(3, 4))
            end
            KeyedArray(data; image_named_axiskeys(header)...)
        else
            nothing
        end

        return FitsImage(; path, header, data)
    end
end

AxisKeys.axiskeys(fi::FitsImage) = image_axiskeys(fi.header)
AxisKeys.named_axiskeys(fi::FitsImage) = image_named_axiskeys(fi.header)
image_axiskeys(fh::FITSHeader) = values(image_named_axiskeys(fh))
function image_named_axiskeys(fh::FITSHeader)
    mul = 1u"°" |> u"mas"
    ra  = axis_vals(fh, "RA---SIN", zero_reference=true) * mul
    dec = axis_vals(fh, "DEC--SIN", zero_reference=true) * mul
    (; ra, dec)
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
