using Unitful, UnitfulAstro
using FITSIO: FITSHeader, FITS, read_header
using Parameters: @with_kw
using DataFrames: select!, ncol
import StatsBase: mad
using AxisKeys
using Utils


@with_kw struct FitsImage{TD, TCC}
    path::String
    header::FITSHeader
    clean_components::TCC
    data::TD
    noise::Union{Nothing, Float32}
end

function load(::Type{FitsImage}, path; read_data = true, read_clean = true)
    FITS(path) do f
        primary = f[1]
        header = read_header(primary)

        data = if read_data
            ra_vals  = axis_vals(header, "RA---SIN", zero_reference=true) .* u"°" .|> u"mas"
            dec_vals = axis_vals(header, "DEC--SIN", zero_reference=true) .* u"°" .|> u"mas"

            @assert header["BSCALE"] == 1 && header["BZERO"] == 0 && header["BUNIT"] == "JY/BEAM"
            data = read(primary)
            data = dropdims(data, dims=(3, 4))
            KeyedArray(data, ra=ra_vals, dec=dec_vals)
        else
            nothing
        end

        comps = if read_clean
            comps = DataFrame(f["AIPS CC"])
            if ncol(comps) > 3
                @assert all(comps[!, Symbol("MAJOR AX")] .== 0) && all(comps[!, Symbol("MINOR AX")] .== 0) && all(comps[!, :POSANGLE] .== 0) && all(comps[!, Symbol("TYPE OBJ")] .== 0)
            end
            select!(comps, :FLUX => :flux, :DELTAX => (x->x*3.6e6) => :ra, :DELTAY => (x->x*3.6e6) => :dec)
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
    return ra_step * u"°" .|> mas
end

function pixel_area(fi::FitsImage)
    ra_step  = axis_dict(fi.header, "RA---SIN")["CDELT"] |> abs
    dec_step = axis_dict(fi.header, "DEC--SIN")["CDELT"] |> abs
    return (ra_step * u"°" .|> mas) * (dec_step * u"°" .|> mas)
end

@with_kw struct ImageBeam{TU}
    major_axis::TU
    minor_axis::TU
    pa::Float64
end

function image_beam(fi::FitsImage)
    return ImageBeam(
        major_axis = fi.header["BMAJ"] * u"°" .|> mas,
        minor_axis = fi.header["BMIN"] * u"°" .|> mas,
        pa  = fi.header["BPA"] |> deg2rad,
    )
end

area(beam::ImageBeam) = π * beam.major_axis/2 * beam.minor_axis/2

const mul = 4 * log(2)

function (beam::ImageBeam)((ra, dec))
    @assert beam.major_axis == beam.minor_axis
    beam_ax = beam.major_axis
    return exp(-(hypot(ra / beam_ax, dec / beam_ax))^2 * mul)
end

frequency(fi::FitsImage) = axis_val(fi.header, "FREQ") * u"Hz" |> u"GHz"
