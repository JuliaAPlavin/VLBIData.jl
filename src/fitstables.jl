# unlike https://github.com/ajwheeler/FITSTables.jl, I eagerly load all columns immediately

const THDU = Union{FITSIO.TableHDU, FITSIO.ASCIITableHDU}

Tables.istable(::Type{<:THDU}) = true

Tables.columnaccess(::Type{<:THDU}) = true
function Tables.columns(thdu::THDU)::NamedTuple
    names = (; (Symbol(n) => n for n in colnames(thdu))...)
    data_orig = map(n -> read(thdu, n), names)
    data_reshaped = map(data_orig) do v
        if ndims(v) == 1
            return v
        else
            # see https://github.com/ajwheeler/FITSTables.jl if more dimensions are needed
            @assert ndims(v) == 2
            return [v[:, i] for i in 1:size(v, 2)]
        end
    end
    return data_reshaped
end
