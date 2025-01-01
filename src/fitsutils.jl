axis_types(fh::FITSHeader) = @p begin
    pairs(fh)
    filtermap() do (k, v)
        occursin(r"^CTYPE\d+$", k) ? v : nothing
    end
end

_ismatch(needle::AbstractString, haystack::AbstractString) = needle == haystack
_ismatch(needle::Regex, haystack::AbstractString) = occursin(needle, haystack)
_ismatch(needle, haystack) = false

function axis_ind(fh::FITSHeader, ctype)
    matching_card = @p pairs(fh) filteronly() do (k, v)
        _ismatch(ctype, v)
    end
    key = matching_card.first
    m = match(r"^CTYPE(\d)$", key)
    return parse(Int, m[1])
end

function axis_dict(fh::FITSHeader, ctype)
    ind = axis_ind(fh, ctype)
    re = Regex("^([A-Z]+)$(ind)\$")
    matching_cards = [match(re, k)[1] => v for (k, v) in pairs(fh) if occursin(re, k)]
    return Dict(matching_cards)
end

function axis_vals(fh::FITSHeader, ctype; zero_reference=false)
    dict = axis_dict(fh, ctype)
    if zero_reference
        return ((1:dict["NAXIS"]) .- dict["CRPIX"]) .* dict["CDELT"]
    else
        return dict["CRVAL"] .+ ((1:dict["NAXIS"]) .- dict["CRPIX"]) .* dict["CDELT"]
    end
end
    
function axis_val(fh::FITSHeader, ctype; zero_reference=false)
    vals = axis_vals(fh, ctype, zero_reference=zero_reference)
    @assert length(vals) == 1
    return first(vals)
end
