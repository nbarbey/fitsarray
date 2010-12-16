import numpy as np
import pyfits
import copy

class InfoArray(np.ndarray):
    """An ndarray supplemented with an header object.

    There is no requirement on the nature of the object
    """
    def __new__(subtype, shape=None, data=None, dtype=float, buffer=None, offset=0,
                strides=None, order=None, header=None):
        if shape is not None:
            obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset,
                                     strides, order)
        elif data is not None:
            obj = np.array(data).view(subtype)
        obj.header = header
        return obj
    def __array_finalize__(self, obj):
        if obj is None: return
        self.header = getattr(obj, 'header', None)
    def copy(self):
        return asinfoarray(copy.copy(self), header=copy.copy(self.header))

class FitsArray(InfoArray):
    """A numpy ndarray supplemented with a pyfits header
    """
    def __new__(subtype, shape=None, data=None, file=None, ext=0, dtype=float, buffer=None, offset=0,
                strides=None, order=None, header=None):
        # various inputs
        if shape is not None:
            obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset,
                                     strides, order)
        elif data is not None:
            obj = np.array(data).view(subtype)
        elif file is not None:
            fits = pyfits.fitsopen(file)
            obj = fits[ext].data.view(subtype)
            header = fits[ext].header
            dtype = get_dtype(header)
        # ensure minimal header is setup
        if header is None:
            header = pyfits.PrimaryHDU(obj).header
        if isinstance(header, dict):
            header = dict2header(header)
        obj.header = header
        enforce_minimal_header(obj)
        return obj
    def update(self, key, value, **kargs):
        self.header.update(key, value, **kargs)
    def tofits(self, filename):
        """Save FitsArray as a fits file
        """
        hdu = pyfits.PrimaryHDU(self, header=self.header)
        hdu.writeto(filename)
    def axes(self):
        axes_list = []
        for n in xrange(self.ndim):
            strn = str(n + 1)
            cdelt = self.header['CDELT' + strn]
            crpix = self.header['CRPIX' + strn]
            crval = self.header['CRVAL' + strn]
            axis_len = self.shape[n]
            u = (np.arange(axis_len) - crpix) * cdelt + crval
            axes_list.append(u)
        return axes_list
    def bin(self, factor, axis=None):
        """Output binned data"""
        if np.isscalar(factor):
            factor = self.ndim * (int(factor),)
        if type(axis) is int:
            axis = (axis,)
        if axis is None:
            axis = np.arange(self.ndim, dtype=int)
        factors = np.ones(self.ndim, dtype=int)
        for f, a in zip(factor, axis):
            factors[a] = f
        out = asfitsarray(self)
        for a, f in enumerate(factors):
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
            out.header['CDELT' + strn] *= f
            out.header['CRPIX' + strn] /= f
            out.header['NAXIS' + strn] /= f
        # normalize
        out /= float(np.prod(factors))
        return out

def enforce_minimal_header(arr):
    minimal_defaults = [('SIMPLE',True,),
                        ('BITPIX', int(bitpix_inv[arr.dtype.name])),
                        ('NAXIS', arr.ndim)]
    for i in xrange(arr.ndim):
        minimal_defaults.append(('NAXIS' + str(i + 1), arr.shape[i]))
    for default in minimal_defaults[::-1]:
        arr.update(default[0], arr.header.get(default[0], default[1]), before=0)

def copy_header(header):
    header_dict = dict(header)
    dtype = get_dtype(header)
    cards = list()
    for k in header_dict:
        try:
            cards.append(pyfits.Card(key=k, value=header_dict[k]))
        except(ValueError):
            try:
                cards.append(pyfits.Card(key=k, value=float(header_dict[k])))
            except(ValueError):
                pass
    return pyfits.Header(cards=cards)

def read_fits_array(filename, ext=0):
    """Reads a fits file and output a FitsArray
    """
    fits = pyfits.fitsopen(filename)
    header = fits[ext].header
    dtype = get_dtype(header)
    data = fits[ext].data.astype(dtype)
    fits_array = FitsArray(data=data, header=header, dtype=dtype)
    return fits_array

def asfitsarray(array, header=None):
    """Returns a view of an ndarray or a subclass as a FitsArray
    """
    header = copy_header(getattr(array, 'header', header))
    if isinstance(header, dict):
        header = dict2header(header)
    return FitsArray(data=array, header=header)

def dict2header(header):
    cards = [pyfits.Card(k, header[k]) for k in header]
    return pyfits.Header(cards=cards)

def asinfoarray(array, header=None):
    """Return a view of an array casted as a FitsArray
    """
    header = copy.copy(getattr(array, 'header', header))
    out = InfoArray(data=array, header=header)
    return out

def zeros(shape, **kargs):
    out = FitsArray(shape, **kargs)
    dtype = kargs.get('dtype', float)
    order = kargs.get('order')
    out[:] = np.zeros(shape, dtype=dtype, order=order)
    return out

def ones(shape, **kargs):
    out = FitsArray(shape, **kargs)
    dtype = kargs.get('dtype', float)
    order = kargs.get('order')
    out[:] = np.ones(shape, dtype=dtype, order=order)
    return out

def get_dtype(header, default=float):
    if header.has_key('BITPIX'):
        b = bitpix[str(header['BITPIX'])]
        return b
    else:
        return default

# fits convention for data types
bitpix = {'8':np.dtype(np.int8).name,
          '16':np.dtype(np.int16).name,
          '32':np.dtype(np.int32).name,
          '-32':np.dtype(np.float32).name,
          '-64':np.dtype(np.float64).name}
bitpix_inv = dict()
for k in bitpix: bitpix_inv[bitpix[k]] = k

def infoarrays2infoarray(arrays):
    """Get a list of InfoArrays and return an info array
    as a concatenation along a new axis
    """
    headers = [a.header for a in arrays]
    [a.resize(a.shape + (1,)) for a in arrays]
    data = np.concatenate(arrays, axis=-1)
    return InfoArray(data=data, header=headers)

def fits2fitsarray(fits, ext=0):
    return hdu2fitsarray(fits[ext])

def hdu2fitsarray(hdu):
    return asfitsarray(hdu.data, header=hdu.header)

def fitsarray_from_header(header):
    if isinstance(header, dict):
        header = dict2header(header)
    shape = [header['NAXIS' + str(i + 1)] for i in xrange(header['NAXIS'])]
    dtype = bitpix[str(int(header['BITPIX']))]
    return FitsArray(shape, header=header, dtype=dtype)
