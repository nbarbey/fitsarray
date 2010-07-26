import numpy as np
import pyfits
import copy

class FitsArray(np.ndarray):
    """A numpy ndarray supplemented with a pyfits header
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
    def tofits(self, filename):
        """Save FitsArray as a fits file
        """
        hdu = pyfits.PrimaryHDU(self, header=self.header)
        hdu.writeto(filename)
    def binning(self, factor=1):
        out = copy.copy(self)
        out[:] = binning(out)
        out.header['CDELT1'] *= factor
        out.header['CDELT2'] *= factor
        return out
