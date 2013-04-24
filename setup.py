#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy

gdal_lib='/usr/lib '
gdal_inc='/usr/include/gdal'
gdal_vers='1.7.0'

warpCopyExt = Extension('pyTiff.warpCopy', ['pyTiff/warpCopy.cpp'], include_dirs=['include',numpy.get_include(),gdal_inc],library_dirs = [gdal_lib],libraries=['gdal'+gdal_vers])

setup (name = 'pyTiff',
       version = '.16',
       author = 'Reid Sawtell',
       author_email = "rwsawtel@mtu.edu",
       url = 'http://wiki.mtri.org/display/mtri/pyDepth#pyDepth',
       description = 'Convenience class for handling geotiffs through GDAL',
       packages = ['pyTiff'],
       package_dir = {'pyTiff': 'pyTiff'},
       ext_modules = [warpCopyExt],
       classifiers = [
                        "Programming Language :: Python",
                        "Programming Language :: Python :: 2.7",
                        "Operating System :: OS Independent",
                        "Development Status :: 2 - Pre-Alpha",
                        "Environment :: Console",
                        "Intended Audience :: Science/Research",
                        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
                        "Natural Language :: English",
                        "Topic :: Scientific/Engineering :: GIS",
                     ]
      )
       