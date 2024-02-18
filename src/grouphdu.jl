struct GroupedHDU <: FITSIO.HDU
    fitsfile::FITSIO.FITSFile
    ext::Int
end

function Base.read(hdu::GroupedHDU)
    FITSIO.assert_open(hdu)
    FITSIO.fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = FITSIO.fits_get_img_size(hdu.fitsfile) |> Tuple
    @assert first(sz) == 0
    h = FITSIO.read_header(hdu)
    ngroups = h["GCOUNT"]
    pcount = h["PCOUNT"]

    pnames = map(1:pcount) do i
        h["PTYPE$i"]
    end
    for i in 1:length(pnames)
        while pnames[i] in pnames[1:i-1]
            pnames[i] = "_" * pnames[i]
        end
    end

    ps_scalezero = ntuple(pcount) do i
        h["PSCAL$i"], h["PZERO$i"]
    end
    b_scalezero = h["BSCALE"], get(h, "BZERO", 0.0)
    haskey(h, "BSZERO") && @warn "FITS header contains BSZERO keyword, probably a typo"

    ptype = NTuple{pcount, Float32}
    result_grp = NamedTuple{Tuple(Symbol.(pnames))}(ntuple(_ -> Float64[], pcount))
    result_data = Array{Float32}(undef, Base.tail(sz)..., ngroups)

    T = FITSIO.type_from_bitpix(FITSIO.fits_get_img_equivtype(hdu.fitsfile))
    @assert T == Float32

    buf_grp = Array{Cfloat}(undef, pcount)
    buf_data = Array{Cfloat}(undef, Base.tail(sz))

    _fill_data!(ngroups, hdu, buf_grp, buf_data, result_grp, result_data, ptype, ps_scalezero, b_scalezero)

    return (; result_grp..., DATA=result_data)
end

function _fill_data!(ngroups, hdu, buf_grp, buf_data, result_grp, result_data, ::Type{ptype}, ps_scalezero, b_scalezero) where {ptype}
    for groupix in 1:ngroups
        status = Ref{Cint}(0)
        @ccall FITSIO.libcfitsio.ffggpe(
            hdu.fitsfile.ptr::Ptr{Cvoid},
            groupix::Clong,
            1::Clong,  # firstelem
            length(buf_grp)::Clong,
            buf_grp::Ptr{Cfloat},
            status::Ref{Cint},
        )::Cint
        FITSIO.fits_assert_ok(status[])

        anynul = Ref{Cint}(0)
        @ccall FITSIO.libcfitsio.ffgpve(
            hdu.fitsfile.ptr::Ptr{Cvoid},
            groupix::Clong,
            1::Clong,  # firstelem
            length(buf_data)::Clong,
            0.0::Cfloat,  # nulval
            buf_data::Ptr{Cfloat},
            anynul::Ref{Cint},
            status::Ref{Cint},
        )::Cint
        @assert anynul[] == 0
        FITSIO.fits_assert_ok(status[])
        buf_data .= muladd.(buf_data, b_scalezero...)

        map(Tuple(result_grp), ptype(buf_grp), ps_scalezero) do vres, val, scalezero
            push!(vres, muladd(val, scalezero...))
        end
        selectdim(result_data, ndims(result_data), groupix) .= buf_data
    end
end
