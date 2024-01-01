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
    pfuncs = map(1:pcount) do i
        x -> muladd(x, h["PSCAL$i"], h["PZERO$i"])
    end
    bfunc = x -> muladd(x, h["BSCALE"], h["BZERO"])

    pnames = map(1:pcount) do i
        h["PTYPE$i"]
    end
    for i in 1:length(pnames)
        while pnames[i] in pnames[1:i-1]
            pnames[i] = "_" * pnames[i]
        end
    end
    ptype = NamedTuple{Tuple(Symbol.(pnames)), NTuple{pcount, Float32}}
    result_grp = NamedTuple{Tuple(Symbol.(pnames))}(ntuple(_ -> Float64[], pcount))
    result_data = Array{Float32}(undef, Base.tail(sz)..., ngroups)

    T = FITSIO.type_from_bitpix(FITSIO.fits_get_img_equivtype(hdu.fitsfile))
    @assert T == Float32

    buf_grp = Array{Cfloat}(undef, pcount)
    buf_data = Array{Cfloat}(undef, Base.tail(sz))

    map(1:ngroups) do groupix
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

        @ccall FITSIO.libcfitsio.ffgpve(
            hdu.fitsfile.ptr::Ptr{Cvoid},
            groupix::Clong,
            1::Clong,  # firstelem
            length(buf_data)::Clong,
            0.0::Cfloat,  # nulval
            buf_data::Ptr{Cfloat},
            Ref{Cint}(0)::Ref{Cint},  # anynul
            status::Ref{Cint},
        )::Cint
        buf_data .= bfunc.(buf_data)

        map(result_grp, ptype(buf_grp), pfuncs) do vres, val, func
            push!(vres, func(val))
        end
        selectdim(result_data, ndims(result_data), groupix) .= buf_data
        # return (; ptype(buf_grp)..., DATA=buf_data)
    end
    return (; result_grp..., DATA=result_data)
end
