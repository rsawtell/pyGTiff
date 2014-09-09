//Copyright 2012 Reid Sawtell

//This file is part of pyGTiff.

//pyGTiff is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//pyGTiff is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with pyGTiff.  If not, see <http://www.gnu.org/licenses/>.

#include <Python.h>
#include <numpy/arrayobject.h>
#include "ogr_api.h"
#include "ogr_geometry.h"
#include "string"

typedef unsigned char byte;

/**
 * Recursive function takes a rectangular region of a raster and checks if the shape intersects it or not
 * If the shape intersect partially and the region contains multiple pixels it is subdivided into quads and each of those are reprocessed
 * If the shape intersects completely each pixel is marked as full area and processing stops
 * If the shape intersects partially and the region is a single pixel the proportion of overlap is recorded and processing stops
 * If the shape does not intersect processing stops
 */
static void processArea(OGRPolygon* shape,double** lat, double** lon, float* tdata, float* odata, const int width, const int height, const double pixArea, int minx, int miny, int maxx, int maxy)
{
    //compute change in each direction
    int deltax = maxx-minx;
    int deltay = maxy-miny;
    
    //short circuit if no pixels to process
    if(deltax==0 || deltay==0) return;
    
    //construct a polygon representing a cluster of pixels from the raster
    OGRLinearRing pixelRing;
    
    pixelRing.addPoint(lon[miny][minx],lat[miny][minx]);
    pixelRing.addPoint(lon[miny][maxx],lat[miny][maxx]);
    pixelRing.addPoint(lon[maxy][maxx],lat[maxy][maxx]);
    pixelRing.addPoint(lon[maxy][minx],lat[maxy][minx]);
    
    pixelRing.closeRings();
    
    OGRPolygon pixelPoly;
    
    pixelPoly.addRing(&pixelRing);
    
    //compute the area
    double regionArea = pixelPoly.get_Area();

    double intArea = 0;
    
    //compute the intersection with the polygon
    if(pixelPoly.Intersects(shape))
    {
        OGRGeometry* psint = pixelPoly.Intersection(shape);
        
        if(psint->getGeometryType()==wkbPolygon)
        {
        
            intArea = ((OGRPolygon*)psint)->get_Area();
            
        }
        else if(psint->getGeometryType()==wkbMultiPolygon)
        {
            intArea = ((OGRMultiPolygon*)psint)->get_Area();
        }
    }
    
    //if only one pixel, record the area of overlap
    if(deltax==1 && deltay==1)
    {
        odata[miny*width+minx] = intArea;
    }
    else
    {
        //if the entire set of pixels is contained within the polygon, record them all as max area per pixel
        if(regionArea==intArea)
        {
            for(int i=miny; i<maxy; i++)
            {
                for(int j=minx; j<maxx; j++)
                {
                   odata[i*width+j] = pixArea; 
                }
            }
        }
        else if(intArea==0)//if none overlap stop processing
        {
            return;
        }
        else
        {
            //split the group of pixels into 4 quandrants and recompute
            processArea(shape,lat,lon,tdata,odata,width,height,pixArea,minx,                  miny,                  (int)(minx+deltax/2.0),(int)(miny+deltay/2.0));
            processArea(shape,lat,lon,tdata,odata,width,height,pixArea,(int)(minx+deltax/2.0),miny,                  maxx,                  (int)(miny+deltay/2.0));
            processArea(shape,lat,lon,tdata,odata,width,height,pixArea,(int)(minx+deltax/2.0),(int)(miny+deltay/2.0),maxx,                  maxy);
            processArea(shape,lat,lon,tdata,odata,width,height,pixArea,minx,                  (int)(miny+deltay/2.0),(int)(minx+deltax/2.0),maxy);
        }
    }
}

/**
 * Get a the area of each pixel overlapping a polygon
 */
static PyObject * shapeSlice(PyObject *self, PyObject *args)
{
    PyArrayObject *transform;
    PyArrayObject *shapebinary;//assumed to be in the same projection as the raster
    int height, width;
    double pixelArea;
    double** lat = NULL;
    double** lon = NULL;
    
    //get arguments
    if (!PyArg_ParseTuple(args, "O!O!ii", &PyArray_Type,&transform,&PyArray_Type,&shapebinary,&width,&height))
    {
        return NULL;
    }
    
    //check for errors
    if (transform==NULL || shapebinary == NULL)
    {
        PyErr_SetString(PyExc_ValueError,"Invalid input arguments.");
        return NULL;
    }
    
    //convert wkb to polygon
    byte* sdata = (byte*)PyArray_DATA(shapebinary);
    OGRPolygon shape = OGRPolygon();
    
    if(shape.importFromWkb(sdata))
    {
        PyErr_SetString(PyExc_ValueError,"Error parsing WKB data.");
        return NULL;
    }
    
    if(!shape.IsValid())
    {
        PyErr_SetString(PyExc_ValueError,"Polygon is invalid, please fix all geometry errors before intersecting.");
        return NULL;
    }
    
    //create the output data array
    int dimensions[2] = {height,width};
    PyArrayObject* outdata = (PyArrayObject*)PyArray_FromDims(2,dimensions,PyArray_FLOAT);
    
    float* odata = (float*) PyArray_DATA(outdata);
    float* tdata = (float*) PyArray_DATA(transform);
    
    //initialize lat/lon coords for all pixel corners
    lat = new double*[height+1];
    lon = new double*[height+1];
    
    for(int i=0;i<=height;i++)
    {
        lat[i] = new double[width+1];
        lon[i] = new double[width+1];
        
        for(int j=0;j<=width;j++)
        {
            lon[i][j] = tdata[0] + tdata[1]*j + tdata[2]*i;
            lat[i][j] = tdata[3] + tdata[4]*j + tdata[5]*i;
        }
    }
    
    //compute pixel area
    OGRLinearRing pixelRing;
    
    pixelRing.addPoint(lon[0][0],lat[0][0]);
    pixelRing.addPoint(lon[0][1],lat[0][1]);
    pixelRing.addPoint(lon[1][1],lat[1][1]);
    pixelRing.addPoint(lon[1][0],lat[1][0]);
    
    pixelRing.closeRings();
    
    OGRPolygon pixelPoly;
    
    pixelPoly.addRing(&pixelRing);
    
    pixelArea = pixelPoly.get_Area();
    
    //compute intersection of pixels and polygon
    processArea(&shape,lat,lon,tdata,odata,width,height,pixelArea,0,0,width,height);
    
    //normalize areas
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            odata[i*width+j] /= pixelArea;
        }
    }
    
    //memory cleanup
    if(lat!=NULL)
    {
        for(int i=0;i<height;i++)
        {
            delete[] lat[i];
        }
        delete[] lat;
    }
    
    if(lon!=NULL)
    {
        for(int i=0;i<height;i++)
        {
            delete[] lon[i];
        }
        delete[] lon;
    }
    
    //return data
    return Py_BuildValue("N",(PyObject*)outdata); 
}

static PyObject * shapeTransform(PyObject *self, PyObject *args)
{
    PyArrayObject *shapebinary;
    const char* srcWKT;
    const char* dstWKT;
    int mode;
    
    //get arguments
    if (!PyArg_ParseTuple(args, "O!ssi", &PyArray_Type,&shapebinary,&srcWKT,&dstWKT,&mode))
    {
        return NULL;
    }
    
    //check for errors
    if (shapebinary == NULL || srcWKT==NULL || dstWKT==NULL)
    {
        PyErr_SetString(PyExc_ValueError,"Invalid input arguments");
        return NULL;
    }
    
    //convert wkb to polygon
    byte* sdata = (byte*)PyArray_DATA(shapebinary);
    OGRPolygon shape = OGRPolygon();
    
    if(shape.importFromWkb(sdata)) return NULL;
    
    OGRSpatialReference srs = OGRSpatialReference();
    
    switch(mode)
    {
        case 0://WKT
            if(srs.importFromWkt((char**)&srcWKT))
            {   
                PyErr_SetString(PyExc_ValueError,(std::string("Invalid WKT: ") + srcWKT).c_str());
                return NULL;
            }
            break;
        case 1://EPSG
            if(srs.importFromEPSG(atoi(srcWKT)))
            {
                PyErr_SetString(PyExc_ValueError,(std::string("Invalid EPSG number: ") + srcWKT).c_str());
                return NULL;
            }
            break;
        case 2://PROJ4
            if(srs.importFromProj4(srcWKT))
            {
                PyErr_SetString(PyExc_ValueError,(std::string("Invalid Proj4 String: ") + srcWKT).c_str());
                return NULL;
            }
            break;
        default:
            PyErr_SetString(PyExc_ValueError,"Mode must be one of: 0 - WKT, 1 - EPSG, 2 - Proj4");
            return NULL;
    }
    
    OGRSpatialReference dst = OGRSpatialReference();
    if(dst.importFromWkt((char**)&dstWKT))
    {
        PyErr_SetString(PyExc_ValueError,"Error determining destination spatial reference.");
        return NULL;
    }
    
    OGRCoordinateTransformation* trans = OGRCreateCoordinateTransformation(&srs,&dst);
    if(trans==NULL)
    {
        PyErr_SetString(PyExc_ValueError,"Could not compute coordinate transformation.");
        return NULL;
    }
    
    if(shape.transform(trans)) return NULL;
    
    int dimensions[1] = {shape.WkbSize()};
    PyArrayObject* outdata = (PyArrayObject*)PyArray_FromDims(1,dimensions,PyArray_BYTE);
    
    shape.exportToWkb(wkbNDR,(byte*)PyArray_DATA(outdata));
    
    //return data
    return Py_BuildValue("N",(PyObject*)outdata);
}

static PyMethodDef shapeMethods[] = 
{
    {"shapeSlice",shapeSlice, METH_VARARGS, "Get a boolean array indexing the intersection of a raster and a shape"},
    {"shapeTransform",shapeTransform,METH_VARARGS, "Translate a shape from one projection into another"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initshape(void)
{
    (void) Py_InitModule("shape", shapeMethods);
    import_array();
}