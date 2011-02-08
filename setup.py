#!/usr/bin/env python
from numpy.distutils.core import Extension, setup
from numpy import get_include
from os.path import join
import sys
setup(name='fitsarray',
      version='0.2.0',
      description='An ndarray subclass with a fits header',
      author='Nicolas Barbey',
      author_email='nicolas.a.barbey@gmail.com',
      requires = ['numpy (>1.3.0)', 'pyfits'],
      packages=['fitsarray'],
      )
