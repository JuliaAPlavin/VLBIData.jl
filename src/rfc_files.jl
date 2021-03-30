@deprecate parse_fits_fname(fname) fname_parse(fname) false
function fname_parse(fname)
    m = match(r"""
        ^
        (?<directory>.*/)?
        (?<J2000>
            J
            (?P<ra_h>\d{2})(?P<ra_m>\d{2})
            (?P<dec_sign>[+-])(?P<dec_d>\d{2})(?P<dec_m>\d{1,2})[A-Z]?
        )
        (_(?<band>[A-Z]{1,3}))?
        (_(?<epoch>\d{4}_\d{2}_\d{2}))?
        (_(?<author>[a-z]{2,3}))?
        (_(?<suffix>[\w-]*))?
        (\.(?<extension>.+))?
        $
    """x, fname)
    keys = (:J2000, :band, :epoch, :author, :suffix, :extension, :directory)
    return NamedTuple{keys, NTuple{length(keys), Union{Nothing, String}}}(map(k -> m[k], keys))
end

@deprecate build_fits_fname(fname) fname_build(fname) false
function fname_build(params)
    parts = [
        get(params, k, nothing)
        for k in [:J2000, :band, :epoch, :author, :suffix]
    ]
    fname = join(filter(!isnothing, parts), "_")
    if get(params, :extension, nothing) !== nothing
        if startswith(params.extension, ".")
            @warn "Provided extensions starts with '.' but probably shouldn't." params.extension
        end
        fname *= ".$(params.extension)"
    end
    if get(params, :directory, nothing) !== nothing
        @assert endswith(params.directory, "/")
        fname = joinpath(params.directory, fname)
    end
    return fname
end

@deprecate replace_in_fname(fname; kwargs...) fname_replace(fname; kwargs...) false
function fname_replace(fname; kwargs...)
    replace(old, (exp, new)::Pair) = (@assert old == exp; new)
    replace(old, new) = new

    parsed = fname_parse(fname)
    merged = if valtype(kwargs) <: Union{Nothing, AbstractString, Symbol}
        merge(parsed, kwargs)
    else
        for (k, v) in kwargs
            new_v = replace(get(parsed, k, nothing), v)
            parsed = merge(parsed, NamedTuple{(k,), Tuple{Union{String, Nothing}}}((new_v,)))
        end
        parsed
    end
    return fname_build(merged)
end

rfc_fname_swap_mapvis(fname::String) = if endswith(fname, "_vis.fits")
    strip_suffix(fname, "_vis.fits") * "_map.fits"
else
    strip_suffix(fname, "_map.fits") * "_vis.fits"
end
