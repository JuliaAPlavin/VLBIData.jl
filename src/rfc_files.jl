using FormatInterp
using Utils
import SplitApplyCombine

function parse_fits_fname(fname)
    m = match(r"""
        ^
        (?<j2000>
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
    return (; [k => m[k] for k in [:j2000, :band, :epoch, :author, :suffix, :extension]]...)
end

function build_fits_fname(params)
    res = f"{params.j2000}_{params.band}_{params.epoch}_{params.author}"
    if params.suffix != nothing
        res *= f"_{params.suffix}"
    end
    if params.extension != nothing
        res *= f".{params.extension}"
    end
    return res
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
