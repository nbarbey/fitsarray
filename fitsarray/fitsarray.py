import numpy as np
import pyfits
import copy

class InfoArray(np.ndarray):
    """An ndarray supplemented with an header object.

    There is no requirement on the nature of the object
    """
    def __new__(subtype, shape, dtype=float, buffer=None, offset=0,
                strides=None, order=None, header=None):
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, 
                                 strides, order)
        obj.header = header
        return obj
    def __array_finalize__(self, obj):
        if obj is None: return
        self.header = getattr(obj, 'header', None)

class FitsArray(InfoArray):
    """A numpy ndarray supplemented with a pyfits header
    """
    def tofits(self, filename):
        """Save FitsArray as a fits file
        """
        hdu = pyfits.PrimaryHDU(self, header=self.header)
        hdu.writeto(filename)
    def axes(self):
        axes_list = []
        for n in xrange(self.ndim):
            strn = str(n + 1)
            cdelt = self.header['cdelt' + strn]
            crpix = self.header['crpix' + strn]
            crval = self.header['crval' + strn]
            axis_len = self.shape[n]
            u = (np.arange(axis_len) - crpix) * cdelt + crval
            axes_list.append(u)
        return axes_list
    def bin(self, factor, axis=None):
        """Output binned data"""
        if type(factor) is int:
            factor = self.ndim * (factor,)
        if type(axis) is int:
            axis = (axis,)
        if axis is None:
            axis = np.arange(self.ndim)
        factors = np.ones(self.ndim)
        axes = np.arange(self.ndim)
        for f, a in zip(factor, axis):
            factors[a] = f
        out = asfitsarray(self)
        for f, a in zip(factors, axes):
            # update data
            for i in xrange(f):
                s = self.ndim * [slice(None), ]
                s[a] = slice(i, None, f)
                if i == 0:
                    tmp = copy.copy(out)
                    out = tmp[s]
                else:
                    out += tmp[s]
            # update metadata
            strn = str((a % self.ndim) + 1)
            out.header['cdelt' + strn] *= f
            out.header['crpix' + strn] /= f
        return out

def copy_header(header):
    header_dict = dict(header)
    cards = [pyfits.Card(key=k, value=header_dict[k]) for k in header_dict]
    return pyfits.Header(cards=cards)

def read_fits_array(filename, ext=0):
    """Reads a fits file and output a FitsArray
    """
    fits = pyfits.fitsopen(filename)
    header = fits[ext].header
    dtype = get_dtype(header)
    data = fits[ext].data.astype(dtype)
    fits_array = FitsArray(data.shape, header=header, dtype=dtype)
    fits_array[:] = data
    return fits_array

def asfitsarray(array):
    """Return a copy of an array casted as a FitsArray
    """
    header = copy_header(getattr(array, 'header', None))
    out = FitsArray(array.shape, header=header)
    out[:] = copy.copy(array)
    return out

def get_dtype(header, default=float):
    if header.has_key('BITPIX'):
        b = bitpix[str(header['BITPIX'])]
        return b
    else:
        return default

# fits convention for data types
bitpix = {'8':np.int8, '16':np.int16, '32':np.int32,
          '-32':np.float32,'-64':np.float64}

def infosarrays2infoarray(arrays):
    """Get a list of InfoArrays and return an info array
    as a concatenation along a new axis
    """
    N = len(arrays)
    keys = arrays[0].header.keys()
    header = dict()
    for k in keys: header[k] = []
    out = InfoArray(arrays[0].shape + (N,), header=header)
    # update values
    out[:] = np.concatenate([array.reshape(array.shape + (1,))
                             for array in arrays], axis=-1)
    # update keys
    for array in arrays:
        for k in keys:
            out.header[k] += [array.header[k],]
    return out
