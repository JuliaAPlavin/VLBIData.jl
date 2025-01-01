Base.@kwdef struct FitsImage{TD}
    path::String
    header::FITSHeader
    data::TD
end

AxisKeys.KeyedArray(fimg::FitsImage) = isnothing(fimg.data) ? error("Image data not loaded") : fimg.data::KeyedArray

load(::Type{KeyedArray}, path) = load(FitsImage, path).data

function load(::Type{FitsImage}, path; read_data=true)
    path = abspath(path)  # for RFC.File
    FITS(path) do f
        primary = f[1]
        header = read_header(primary)

        data = if read_data
            haskey(header, "BSCALE") && @assert header["BSCALE"] == 1
            haskey(header, "BZERO" ) && @assert header["BZERO" ] == 0
            uppercase(header["BUNIT"]) ∈ ("JY/BEAM", "JY/PIXEL") || error("Unknown BUNIT: $(header["BUNIT"])")
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
    ra  = axis_vals(fh, r"RA-+SIN", zero_reference=true) * mul
    dec = axis_vals(fh, r"DEC-+SIN", zero_reference=true) * mul
    (; ra, dec)
end

function pixel_size(fi::FitsImage)
    ra_step  = axis_dict(fi.header, r"RA-+SIN")["CDELT"] |> abs
    dec_step = axis_dict(fi.header, r"DEC-+SIN")["CDELT"] |> abs
    @assert ra_step == dec_step
    return ra_step * u"°" .|> u"mas"
end

function pixel_steps(fi::FitsImage)
    ra_step  = axis_dict(fi.header, r"RA-+SIN")["CDELT"]
    dec_step = axis_dict(fi.header, r"DEC-+SIN")["CDELT"]
    @assert abs(ra_step) == abs(dec_step)
    return SVector(ra_step, dec_step) .* u"°" .|> u"mas"
end

function pixel_area(fi::FitsImage)
    ra_step  = axis_dict(fi.header, r"RA-+SIN")["CDELT"] |> abs
    dec_step = axis_dict(fi.header, r"DEC-+SIN")["CDELT"] |> abs
    return (ra_step * u"°" .|> u"mas") * (dec_step * u"°" .|> u"mas")
end

load(::Type{Beam}, src) = FITS(abspath(src)) do f
    header = read_header(f[1])
    Beam(header)
end
InterferometricModels.beam(x::Union{FitsImage,FITSHeader}) = Beam(x)
InterferometricModels.Beam(fi::FitsImage) = Beam(fi.header)
InterferometricModels.Beam(fh::FITSHeader) = Beam(EllipticGaussian,
    σ_major = InterferometricModels.fwhm_to_σ(fh["BMAJ"] * u"°" .|> u"mas"),
    ratio_minor_major = fh["BMIN"] / fh["BMAJ"],
    pa_major = fh["BPA"] |> deg2rad,
)

frequency(fi::FitsImage) = axis_val(fi.header, "FREQ") * u"Hz" |> u"GHz"


function prepare_header(image::KeyedArray{<:Any,2}; unit::String="JY/PIXEL", freq::Union{Real,Nothing}=nothing)
    Unitful.unit(eltype(image)) == u"Jy" || error("Only images with Jy units are supported, corresponding to Jy/pix; got $(Unitful.unit(eltype(image)))")
    axks = axiskeys(image)
    pixsizes = step.(axks)
    map(NamedTuple{(:key,:value,:comment)}, [
        ("SIMPLE"  , true           , "conforms to FITS standard" ),
        ("BITPIX"  , -64            , "array data type"           ),
        ("NAXIS"   , 2              , "number of array dimensions"),
        ("NAXIS1"  , size(image, 1) , ""                          ),
        ("NAXIS2"  , size(image, 2) , ""                          ),
        ("EXTEND"  , true           , ""                          ),
        # ("OBJECT"  , source    , ""                          ),
        ("CTYPE1"  , "RA---SIN"     , ""                          ),
        ("CTYPE2"  , "DEC--SIN"     , ""                          ),
        ("CDELT1"  , ustrip(u"°", pixsizes[1]), ""                          ),
        ("CDELT2"  , ustrip(u"°", pixsizes[2]), ""                          ),
        # ("OBSRA"   , 180.0        , ""                          ),
        # ("OBSDEC"  , 0.0       , ""                          ),
        ("FREQ"    , freq         , ""                          ),
        ("CRPIX1"  , -(first(axks[1])/pixsizes[1] - 1), ""        ),
        ("CRPIX2"  , -(first(axks[2])/pixsizes[2] - 1), ""        ),
        # ("MJD"     , 0       , ""                          ),
        # ("TELESCOP", "VLBI"         , ""                          ),
        ("BUNIT"   , unit           , ""                          ),
        ("STOKES"  , "I"         , ""                          ),
    ]) |> FITSHeader
end

function save(fname, image::KeyedArray{<:Any,2}; kwargs...)
    rm(fname, force=true)
    FITS(fname, "w") do f
        try
            # img = parent(image[end:-1:1, :])  # XXX: is it needed?
            FITSIO.write(f, ustrip.(u"Jy", AxisKeys.keyless_unname(image)); header=prepare_header(image; kwargs...))
        catch e
            @error "Got exception when writing image data to FITS" e
            rethrow()
        end
    end
end
