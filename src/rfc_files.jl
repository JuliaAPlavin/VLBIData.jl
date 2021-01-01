using FormatInterp
using Utils
import SplitApplyCombine

@deprecate parse_fits_fname(fname) fname_parse(fname) false
function fname_parse(fname)
    m = match(r"""
        ^
        (?<J2000>
            J
            (?P<ra_h>\d{2})(?P<ra_m>\d{2})
            (?P<dec_sign>[+-])(?P<dec_d>\d{2})(?P<dec_m>\d{1,2})[A-Z]?
        )_
        (?<band>[A-Z])_
        (?<epoch>\d{4}_\d{2}_\d{2})_
        (?<author>[a-z]{2,3})_
        (?<suffix>\w*)
        \.(?<extension>\w+)
        $
    """x, basename(fname))
    return (; [k => m[k] for k in [:J2000, :band, :epoch, :author, :suffix, :extension]]...)
end

@deprecate build_fits_fname(fname) fname_build(fname) false
function fname_build(params)
    parts = [
        get(params, k, nothing)
        for k in [:J2000, :band, :epoch, :author, :suffix]
    ]
    filter!(!isnothing, parts)
    fname = join(parts, "_")
    if get(params, :extension, nothing) != nothing
        fname *= ".$(params.extension)"
    end
    return fname
end

@deprecate replace_in_fname(fname; kwargs...) fname_replace(fname; kwargs...) false
function fname_replace(fname; kwargs...)
    replace(old, (exp, new)::Pair) = (@assert old == exp; new)
    replace(old, new) = new

    parsed = fname_parse(fname)
    for (k, v) in kwargs
        new_v = replace(get(parsed, k, nothing), v)
        parsed = merge(parsed, (k => new_v,))
    end
    return fname_build(parsed)
end

rfc_fname_swap_mapvis(fname::String) = if endswith(fname, "_vis.fits")
    strip_suffix(fname, "_vis.fits") * "_map.fits"
else
    strip_suffix(fname, "_map.fits") * "_vis.fits"
end

function match_fits_files(
        files_a::Vector{<:AbstractString}, files_b::Vector{<:AbstractString};
        by, check_same=nothing, check_unique::Bool, check_allmatch::Bool
    )
    get_key(by::Symbol, parsed_fname::NamedTuple) = getfield(parsed_fname, by)
    get_key(by::NTuple{N, Symbol}, parsed_fname::NamedTuple) where {N} = getfield.(Ref(parsed_fname), by)
    get_key(by::Function, parsed_fname::NamedTuple) = by(parsed_fname)
    get_key(by, fname::String) = get_key(by, parse_fits_fname(fname))

    for fs in [files_a, files_b]
        parsed = parse_fits_fname.(fs)
        @assert !any(isnothing.(parsed))
        if check_same != nothing
            @assert get_key.(Ref(check_same), parsed) |> unique |> length == 1
        end
        if check_unique
            @assert allunique(get_key.(Ref(by), parsed))
        end
    end

    ab_pairs = SplitApplyCombine.innerjoin(
        f -> get_key(by, f),
        f -> get_key(by, f),
        (vf, mf) -> (vf, mf),
        files_a,
        files_b,
    )
    @info "Fits files matching:" files_a=length(files_a) files_b=length(files_b) pairs=length(ab_pairs)
    if check_allmatch
        @assert length(ab_pairs) == min(length(files_a), length(files_b)) (length(files_a), length(files_b), length(ab_pairs))
    end
    return ab_pairs
end
