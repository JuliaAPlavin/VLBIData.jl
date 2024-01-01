module PyCallExt

using PyCall
using VLBIData: VLBIData, UVData

function VLBIData.read_data_raw(uvdata::UVData, ::typeof(pyimport))
    pyimport_conda("numpy", "numpy")
    pyimport_conda("astropy", "astropy")
    py"""
    import numpy as np
    import astropy.io.fits

    def to_native_byteorder(arr):
        dt = arr.dtype.newbyteorder('=')
        return arr.astype(dt)

    with astropy.io.fits.open(open($(uvdata.path), 'rb'), memmap=False) as f:
        raw = f[0].data
    """
    raw = Dict{Symbol,Any}(Symbol(n) => PyArray(py"to_native_byteorder(raw[$n])"o) for n in py"raw.dtype.names")
    # raw[:DATA] = permutedims(raw[:DATA], reverse(1:ndims(raw[:DATA])))
    return raw
end

end
