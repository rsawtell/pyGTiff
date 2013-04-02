#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy

setup (name = 'pyTiff',
       version = '.15',
       author = 'Reid Sawtell',
       author_email = "rwsawtel@mtu.edu",
       url = 'http://wiki.mtri.org/display/mtri/pyDepth#pyDepth',
       description = 'Convenience class for handling geotiffs through GDAL',
       packages = ['pyTiff'],
       package_dir = {'pyTiff': 'pyTiff'},
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
       