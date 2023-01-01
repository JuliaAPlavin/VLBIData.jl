function guess_type(src)
    if startswith(first(eachline(src)), r"! \w+")
        MultiComponentModel
    else
        ctypes = try
            FITS(src) do f
                read_header(f[1]) |> axis_types
            end
        catch e
            error("Cannot read $src as a FITS file")
        end
        if "RA---SIN" ∈ ctypes
            return FitsImage
        elseif "COMPLEX" ∈ ctypes
            return UVData
        else
            error("Cannot guess data type using FITS header of $src")
        end
    end
end

load(src; kwargs...) = load(guess_type(src), src; kwargs...)
