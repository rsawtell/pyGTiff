#!/usr/bin/env python
from distutils.core import setup, Extension
from distutils.command.build import build
import numpy, sys, platform

#update to reflect your system
if platform.system() == 'Windows':
    gdal_lib='C:\OSGeo4W\lib'
    gdal_inc='C:\OSGeo4W\include'
    gdal_vers='_i'
else:
    gdal_lib='/usr/lib'
    gdal_inc='/usr/include/gdal'
    gdal_vers=''

warpCopyExt = Extension('pyGTiff.warpCopy', ['pyGTiff/warpCopy.cpp'], include_dirs=['include',numpy.get_include(),gdal_inc],library_dirs = [gdal_lib],libraries=['gdal'+gdal_vers])
shapeExt = Extension('pyGTiff.shape', ['pyGTiff/shape.cpp'], include_dirs=['include',numpy.get_include(),gdal_inc],library_dirs = [gdal_lib],libraries=['gdal'+gdal_vers])

#compile without C++ extensions
if( "-noc" in sys.argv or "--noc" in sys.argv):
    
    if '-noc' in sys.argv:
        index = sys.argv.index("-noc")
        
    if '--noc' in sys.argv:
        index = sys.argv.index("--noc")  
        
    sys.argv.pop(index)
    
    
    setup (name = 'pyGTiff',
        version = '1.0.5',
        author = 'Reid Sawtell',
        author_email = "rwsawtel@mtu.edu",
        url = 'https://bitbucket.org/rsawtell/pygtiff',
        description = 'Convenience class for handling geotiffs (and potentially other formats) through GDAL',
        packages = ['pyGTiff'],
        package_dir = {'pyGTiff': 'pyGTiff'},
        classifiers = [
                            "Programming Language :: Python",
                            "Programming Language :: Python :: 2.7",
                            "Operating System :: OS Independent",
                            "Development Status :: 5 - Production/Stable",
                            "Environment :: Console",
                            "Intended Audience :: Science/Research",
                            "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
                            "Natural Language :: English",
                            "Topic :: Scientific/Engineering :: GIS",
                        ]
        )
        
#compile with C++ extensions        
else:
    setup (name = 'pyGTiff',
        version = '1.0.4',
        author = 'Reid Sawtell',
        author_email = "rwsawtel@mtu.edu",
        url = 'https://bitbucket.org/rsawtell/pygtiff',
        description = 'Convenience class for handling geotiffs (and potentially other formats) through GDAL',
        packages = ['pyGTiff'],
        package_dir = {'pyGTiff': 'pyGTiff'},
        ext_modules = [warpCopyExt, shapeExt],
        classifiers = [
                            "Programming Language :: Python",
                            "Programming Language :: Python :: 2.7",
                            "Operating System :: OS Independent",
                            "Development Status :: 5 - Production/Stable",
                            "Environment :: Console",
                            "Intended Audience :: Science/Research",
                            "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
                            "Natural Language :: English",
                            "Topic :: Scientific/Engineering :: GIS",
                        ]
        )
       
