Developed at the Michigan Tech Research Institute (MTRI), pyGTiff is a convenience class for handling geotiffs 
(and potentially other formats) through GDAL in Python. This project is intended to facilitate scientific 
computing of geospatial information by making GDAL easier to use for scientists who are not primarily programmers.

Installation Instructions
---------------------------


UNIX:

If you haven't done so already, extract the archive contents to a temporary location or clone the repository

Open setup.py and change the following to match your setup: 

gdal_lib='/usr/lib '
gdal_inc='/usr/include/gdal'
gdal_vers='' #not usually needed, but specify if it can't find the library

Then run
python setup.py build [-noc, --noc]
sudo python setup.py install (sudo is only needed if you are installing it to the default python installation)

The noc option specifies that the supplementary C++ libraries should not be built.
This will disable some functionality (shapeIntersect, reproject and intersect)


WINDOWS:

Install osgeo4w: http://trac.osgeo.org/osgeo4w/
Download one of the binary distributions and install

Or from source:

If you haven't done so already, extract the source archive contents to a temporary location or clone the repository
Download and install Visual Studio 2008 Express
open the OSGeo4W terminal
cd to the location you put the source
python setup.py install