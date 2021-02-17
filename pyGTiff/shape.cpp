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
#include "list"

typedef unsigned char byte;

/**
 * Recursive function takes a rectangular region of a raster and checks if the shape intersects it or not
 * If the shape intersect partially and the region contains multiple pixels it is subdivided into quads and each of those are reprocessed
 * If the shape intersects completely each pixel is marked as full area and processing stops
 * If the shape intersects partially and the region is a single pixel the proportion of overlap is recorded and processing stops
 * If the shape does not intersect processing stops
 */
static void processArea(OGRGeometry* shape,double** lat, double** lon, PyArrayObject* tdata, PyArrayObject* odata, const int width, const int height, const double pixArea,const int minx,const int miny,const int maxx,const int maxy)
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
    
    OGRGeometry* psint = NULL;
    
    //compute the intersection with the polygon
    if(pixelPoly.Intersects(shape))
    {
        psint = pixelPoly.Intersection(shape);
        
        if(psint->getGeometryType()==wkbPolygon)
        {
            intArea = ((OGRPolygon*)psint)->get_Area(); 
        }
        else if(psint->getGeometryType()==wkbMultiPolygon)
        {
            intArea = ((OGRMultiPolygon*)psint)->get_Area();
        }
        else if(psint->getGeometryType()==wkbGeometryCollection)
        {
            intArea = ((OGRGeometryCollection*)psint)->get_Area();
            //printf("%d,%d - %d,%d\n",minx,miny,maxx,maxy);
        }
    }
    else
    {
        //printf("%d,%d - %d,%d --\n",minx,miny,maxx,maxy);
        return;
    }
    
    //printf("%d,%d - %d,%d %f\n",minx,miny,maxx,maxy,intArea/regionArea);
    
    //if only one pixel, record the area of overlap
    if(deltax==1 && deltay==1)
    {
        //odata[miny*width+minx] = intArea;
        *((float*)PyArray_GETPTR2(odata,miny,minx)) = intArea/regionArea;
        //printf(">%d,%d -> %d,%d : %f -- %f\n",miny,minx,maxy,maxx,intArea/pixArea,regionArea/pixArea);
    }
    else
    {
        //printf("%d,%d :: %d,%d R:%f I:%f  %d\n",minx,maxx,miny,maxy,regionArea/pixArea,intArea/pixArea,regionArea==intArea);
        //if the entire set of pixels is contained within the polygon, record them all as max area per pixel
        if(regionArea==intArea)
        {
            for(int i=miny; i<maxy; i++)
            {
                for(int j=minx; j<maxx; j++)
                {
                   //odata[i*width+j] = pixArea; 
                    *((float*)PyArray_GETPTR2(odata,i,j)) = 1.;//pixArea;
                    //printf("%d,%d : %f\n",i,j,pixArea/pixArea);
                }
            }
        }
        else if(intArea==0)//if none overlap stop processing
        {
            printf("%d,%d - %d,%d\n",minx,miny,maxx,maxy);
            return;
        }
        else
        {
            /*if(psint == NULL)
            {
                psint = shape;
            }*/
            if(psint != NULL)
            {
                //split the group of pixels into 4 quandrants and recompute
                processArea(psint,lat,lon,tdata,odata,width,height,pixArea,minx,                  miny,                  minx+std::max(1,(int)(deltax/2.0)),miny+std::max(1,(int)(deltay/2.0)));
                
                processArea(psint,lat,lon,tdata,odata,width,height,pixArea,(int)(minx+deltax/2.0),miny,                  maxx,                  miny+std::max(1,(int)(deltay/2.0)));
                
                processArea(psint,lat,lon,tdata,odata,width,height,pixArea,(int)(minx+deltax/2.0),(int)(miny+deltay/2.0),maxx,                  maxy);
                
                processArea(psint,lat,lon,tdata,odata,width,height,pixArea,minx,                  (int)(miny+deltay/2.0),minx+std::max(1,(int)(deltax/2.0)),maxy);

                delete psint;
            }
        }
    }
}

void TraceExterior(OGRPolygon* pg, int* buffer, double** lat, double** lon, const int height, const int width, const int si, const int sj, const int id)
{
    int ci = si;
    int cj = sj;
    int cval;
    int cval2;
    
    OGRLinearRing pixelRing;
    bool notfirst = true;
    
    int t = 0;
    int face = 0;
    int f=0;
    //printf("start: %d,%d\n",si,sj);
    do
    {
        //printf("sf: %d\n",face);
        for(f=0;f<4;f++)
        {
            //look left
            if((face + f)%4 == 0)
            {
                if(t>0 && ci==si && cj==sj)
                {
                    f += 1;
                    break;
                }
                
                //printf("left\n");
                if(cj==0)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*ci+cj-1];
                    if(ci==0)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci-1)+cj-1];
                    }
                }
                
                
                if(cval == id)
                {
                    cj -= 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci+1][cj],lat[ci+1][cj]);
                    }
                    notfirst = true;
                    pixelRing.addPoint(lon[ci][cj],lat[ci][cj]);
                    
                    if(cj!=0)
                    {
                        buffer[height*ci+cj-1] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci -= 1;
                        cj -= 1;
                        break;
                    }
                }
            }
            else if((face + f)%4 == 1)
            {
                //printf("up\n");
                //look up
                if(ci==0)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*(ci-1)+cj];
                    if(cj==width-1)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci-1)+cj+1];
                    }
                }
                
                if(cval == id)
                {
                    ci -= 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci][cj],lat[ci][cj]);
                    }
                    notfirst = true;
                    pixelRing.addPoint(lon[ci][cj+1],lat[ci][cj+1]);
                    
                    if(ci!=0)
                    {
                        buffer[height*(ci-1)+cj] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci -= 1;
                        cj += 1;
                        break;
                    }
                }
            }
            else if((face + f)%4 == 2)
            {
                //printf("right\n");
                //look right
                if(cj==width-1)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*ci+cj+1];
                    if(ci==height-1)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci+1)+cj+1];
                    }
                }
                
                if(cval == id)
                {
                    cj += 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci][cj+1],lat[ci][cj+1]);
                    }
                    notfirst = true;
                    pixelRing.addPoint(lon[ci+1][cj+1],lat[ci+1][cj+1]);
                    
                    if(cj!=width-1)
                    {
                        buffer[height*ci+cj+1] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci += 1;
                        cj += 1;
                        break;
                    }
                }
            }
            else if((face + f)%4 == 3)
            {
                //printf("down\n");
                //look down
                if(ci==height-1)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*(ci+1)+cj];
                    if(cj==0)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci+1)+cj-1];
                    }
                }
                
                if(cval == id)
                {
                    ci += 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci+1][cj+1],lat[ci+1][cj+1]);
                    }
                    pixelRing.addPoint(lon[ci+1][cj],lat[ci+1][cj]);
        
                    if(ci!=height-1)
                    {
                        buffer[height*(ci+1)+cj] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci += 1;
                        cj -= 1;
                        break;
                    }
                }
            }
        }
        
        face = (face + f +3)%4;
        
        notfirst = false;
        //printf("\n%d,%d\n",ci,cj);
        t++;
    }
    while(ci != si || cj != sj || face != 0);
    //printf("Done\n");
    
    pixelRing.closeRings();
    //printf("Closed\n");
    
    pg->addRing(&pixelRing);
    //printf("Added\n");
}

void TraceInterior(OGRPolygon* pg, int* buffer, double** lat, double** lon, const int height, const int width, const int si, const int sj, const int id, const int sf)
{
    int ci = si;
    int cj = sj;
    int cval;
    int cval2;
    
    OGRLinearRing pixelRing;
    bool notfirst = true;
    
    int t = 0;
    int face = sf;
    int f=0;
    //printf("start: %d,%d\n",si,sj);
    //printf("sf: %d\n",face);
    do
    {
        for(f=0;f<4;f++)
        {
            if(t>0 && ci==si && cj==sj && (face + f)%4 == sf)
            {
                f += 1;
                break;
            }
            
            //look left
            if((face + f)%4 == 0)
            {
                
                //printf("left\n");
                if(cj==0)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*ci+cj-1];
                    if(ci==0)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci-1)+cj-1];
                    }
                }
                
                
                if(cval == id)
                {
                    cj -= 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci+1][cj],lat[ci+1][cj]);
                    }
                    notfirst = true;
                    pixelRing.addPoint(lon[ci][cj],lat[ci][cj]);
                    
                    if(cj!=0)
                    {
                        buffer[height*ci+cj-1] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci -= 1;
                        cj -= 1;
                        break;
                    }
                }
            }
            else if((face + f)%4 == 1)
            {
                //printf("up\n");
                //look up
                if(ci==0)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*(ci-1)+cj];
                    if(cj==width-1)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci-1)+cj+1];
                    }
                }
                
                if(cval == id)
                {
                    ci -= 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci][cj],lat[ci][cj]);
                    }
                    notfirst = true;
                    pixelRing.addPoint(lon[ci][cj+1],lat[ci][cj+1]);
                    
                    if(ci!=0)
                    {
                        buffer[height*(ci-1)+cj] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci -= 1;
                        cj += 1;
                        break;
                    }
                }
            }
            else if((face + f)%4 == 2)
            {
                //printf("right\n");
                //look right
                if(cj==width-1)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*ci+cj+1];
                    if(ci==height-1)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci+1)+cj+1];
                    }
                }
                
                if(cval == id)
                {
                    cj += 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci][cj+1],lat[ci][cj+1]);
                    }
                    notfirst = true;
                    pixelRing.addPoint(lon[ci+1][cj+1],lat[ci+1][cj+1]);
                    
                    if(cj!=width-1)
                    {
                        buffer[height*ci+cj+1] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci += 1;
                        cj += 1;
                        break;
                    }
                }
            }
            else if((face + f)%4 == 3)
            {
                //printf("down\n");
                //look down
                if(ci==height-1)
                {
                    cval = 0;
                    cval2 = 0;
                }
                else
                {
                    cval = buffer[height*(ci+1)+cj];
                    if(cj==0)
                    {
                        cval2=0;
                    }
                    else
                    {
                        cval2 = buffer[height*(ci+1)+cj-1];
                    }
                }
                
                if(cval == id)
                {
                    ci += 1;
                    break;
                }
                else
                {
                    if(notfirst)
                    {
                        pixelRing.addPoint(lon[ci+1][cj+1],lat[ci+1][cj+1]);
                    }
                    pixelRing.addPoint(lon[ci+1][cj],lat[ci+1][cj]);
        
                    if(ci!=height-1)
                    {
                        buffer[height*(ci+1)+cj] = -id; 
                    }
                    
                    if(cval2 == id)
                    {
                        ci += 1;
                        cj -= 1;
                        break;
                    }
                }
            }
        }
        
        face = (face + f +3)%4;
        
        notfirst = false;
        //printf("\n%d,%d\n",ci,cj);
        t++;
    }
    while(ci != si || cj != sj || face != sf);
    //printf("Done\n");
    
    pixelRing.closeRings();
    //printf("Closed\n");
    
    pg->addRing(&pixelRing);
    //printf("Added\n");
}

void TraceInteriors(OGRPolygon* pg, int* buffer, std::list<std::pair<int,int>> &edgelist, double** lat, double** lon, const int height, const int width, const int id)
{
    int ci;
    int cj;
    
    while(edgelist.size()>0)
    {
        ci = edgelist.front().first;
        cj = edgelist.front().second;
        edgelist.pop_front();
        
        if(buffer[height*ci+cj]>-id && buffer[height*ci+cj]<=0)//not processed yet and freespace
        {
            int sf = 0;
            
            if(ci>0 && buffer[height*(ci-1)+cj] == id)
            {
                //start up
                 TraceInterior(pg,buffer,lat,lon,height,width,ci-1,cj,id,3);
            }
            else if(ci<height-1 && buffer[height*(ci+1)+cj] == id)
            {
                //start down
                TraceInterior(pg,buffer,lat,lon,height,width,ci+1,cj,id,1);
            }
            else if(cj>0 && buffer[height*ci+cj-1] == id)
            {
                //start left
                TraceInterior(pg,buffer,lat,lon,height,width,ci,cj-1,id,2);
            }
            else if(cj<width-1 && buffer[height*ci+cj+1] == id)
            {
                //start right
                TraceInterior(pg,buffer,lat,lon,height,width,ci,cj+1,id,0);
            }
            else
            {
                printf("Error: Couldn't Find Adjacent Polygon ID %d\n",id);
            }           
        }
    }
}

void FloodFill(int* buffer,std::list<std::pair<int,int>> &edgelist,const int height,const int width,const int i,const int j,const int id)
{
    int ci = i;
    int cj = j;
    int cval;
    
    std::list<std::pair<int,int>> nextpixel;
    
    nextpixel.emplace_back(ci,cj);
    buffer[height*ci+cj] = id;
    
    while(nextpixel.size()>0)
    {
        ci = nextpixel.front().first;
        cj = nextpixel.front().second;
        //printf("FF: %d,%d R:%d\n",ci,cj,nextpixel.size());
        
        nextpixel.pop_front();
        
        //up
        if(ci>0)
        {
            cval = buffer[height*(ci-1)+cj];
            if(cval == 1)
            {
                buffer[height*(ci-1)+cj] = id;
                nextpixel.emplace_back(ci-1,cj);
            }
            else if(cval == 0)
            {
                edgelist.emplace_back(ci-1,cj);
            }
        }
        
        //down
        if(ci<height-1)
        {
            cval = buffer[height*(ci+1)+cj];
            if(cval == 1)
            {
                buffer[height*(ci+1)+cj] = id;
                nextpixel.emplace_back(ci+1,cj);
            }
            else if(cval == 0)
            {
                edgelist.emplace_back(ci+1,cj);
            }
        }
        
        //left
        if(cj>0)
        {
            cval = buffer[height*ci+cj-1];
            if(cval == 1)
            {
                buffer[height*ci+cj-1] = id;
                nextpixel.emplace_back(ci,cj-1);
            }
            else if(cval == 0)
            {
                edgelist.emplace_back(ci,cj-1);
            }
        }
        
        //right
        if(cj<width-1)
        {
            cval = buffer[height*ci+cj+1];
            if(cval == 1)
            {
                buffer[height*ci+cj+1] = id;
                nextpixel.emplace_back(ci,cj+1);
            }
            else if(cval == 0)
            {
                edgelist.emplace_back(ci,cj+1);
            }
        }
    }
}

OGRGeometry* mask2polygon(double** lat, double** lon, PyArrayObject* mask, const int width, const int height)
{
    int* buffer = new int[height*width];
    
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            buffer[height*i+j] = *((bool*)PyArray_GETPTR2(mask,i,j));
        }
    }
    
    OGRMultiPolygon* mp = new OGRMultiPolygon();
    
    int id = 2;
    
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            if(buffer[height*i+j] == 1) //start a new polygon
            {
                OGRPolygon* pg = new OGRPolygon();
                std::list<std::pair<int,int>> edgelist;
                
                //floodfill to find interiors
                FloodFill(buffer,edgelist,height,width,i,j,id);
                
                //trace the exterior
                TraceExterior(pg,buffer,lat,lon,height,width,i,j,id);
                
                //printf("%d edges to check for id %d\n",edgelist.size(),id);
                
                //trace interiors
                TraceInteriors(pg,buffer,edgelist,lat,lon,height,width,id);
                
                
                id++;
                
                mp->addGeometry(pg);
                delete pg;
                
                //if(id>3)break;
            }
        }
        //break;
    }
    
    return mp;
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
    unsigned char* sdata = (unsigned char*)PyArray_DATA(shapebinary);
    OGRGeometry* shape;
    
    if(OGRGeometryFactory::createFromWkb(sdata,NULL,&shape,-1))
    {
        PyErr_SetString(PyExc_ValueError,"Error parsing WKB data.");
        return NULL;
    }
    
    //printf("%d\n",shape->getGeometryType());
    if(!(shape->getGeometryType()==wkbPolygon || shape->getGeometryType()==wkbMultiPolygon || shape->getGeometryType()==wkbGeometryCollection))
    {
        OGRGeometryFactory::destroyGeometry(shape);
        PyErr_SetString(PyExc_ValueError,"Shape must be a Polygon or Multi-Polygon.");
        return NULL;
    }
    
    if(!shape->IsValid())
    {
        OGRGeometryFactory::destroyGeometry(shape);
        PyErr_SetString(PyExc_ValueError,"Shape is invalid, please fix all geometry errors before intersecting.");
        return NULL;
    }
    
    //create the output data array
    npy_intp dimensions[2] = {height,width};
    PyArrayObject* outdata = (PyArrayObject*)PyArray_ZEROS(2,dimensions,PyArray_FLOAT,0);
    
    //float* odata = (float*) PyArray_DATA(outdata);
    //float* tdata = (float*) PyArray_DATA(transform);
    
    //initialize lat/lon coords for all pixel corners
    lat = new double*[height+1];
    lon = new double*[height+1];
    
    for(int i=0;i<=height;i++)
    {
        lat[i] = new double[width+1];
        lon[i] = new double[width+1];
        
        for(int j=0;j<=width;j++)
        {
            lon[i][j] = *((double*)PyArray_GETPTR1(transform,0)) + *((double*)PyArray_GETPTR1(transform,1))*j + *((double*)PyArray_GETPTR1(transform,2))*i;
            lat[i][j] = *((double*)PyArray_GETPTR1(transform,3)) + *((double*)PyArray_GETPTR1(transform,4))*j + *((double*)PyArray_GETPTR1(transform,5))*i;
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
    processArea(shape,lat,lon,transform,outdata,width,height,pixelArea,0,0,width,height);
    
    //normalize areas
    /*for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            //odata[i*width+j] /= pixelArea;
            *((float*)PyArray_GETPTR2(outdata,i,j)) /= pixelArea;
        }
    }*/
    
    //memory cleanup
    OGRGeometryFactory::destroyGeometry(shape);
    
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

/**
 * Convert boolean mask into polygon
 */
static PyObject * polygonize(PyObject *self, PyObject *args)
{
    PyArrayObject *transform;
    PyArrayObject *mask;
    int height, width;
    double pixelArea;
    double** lat = NULL;
    double** lon = NULL;
    
    //get arguments
    if (!PyArg_ParseTuple(args, "O!O!ii", &PyArray_Type,&transform,&PyArray_Type,&mask,&width,&height))
    {
        return NULL;
    }
    
    //check for errors
    if (transform==NULL || mask == NULL)
    {
        PyErr_SetString(PyExc_ValueError,"Invalid input arguments.");
        return NULL;
    }
    
    //initialize lat/lon coords for all pixel corners
    lat = new double*[height+1];
    lon = new double*[height+1];
    
    for(int i=0;i<=height;i++)
    {
        lat[i] = new double[width+1];
        lon[i] = new double[width+1];
        
        for(int j=0;j<=width;j++)
        {
            lon[i][j] = *((float*)PyArray_GETPTR1(transform,0)) + *((float*)PyArray_GETPTR1(transform,1))*j + *((float*)PyArray_GETPTR1(transform,2))*i;
            lat[i][j] = *((float*)PyArray_GETPTR1(transform,3)) + *((float*)PyArray_GETPTR1(transform,4))*j + *((float*)PyArray_GETPTR1(transform,5))*i;
        }
    }
    
    OGRGeometry* result = mask2polygon(lat, lon, mask, width, height);
    
    /*if(!result->IsValid())
    {
       OGRGeometry* fixed = OGRGeometry::MakeValid(result);//result->MakeValid();
       delete result;
       result = fixed;
    }*/
    
    if(result == nullptr)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }
    
    //create the output data array
    npy_intp dimensions[1] = {result->WkbSize()};
    PyArrayObject* outdata = (PyArrayObject*)PyArray_ZEROS(1,dimensions,PyArray_BYTE,0);
    
    result->exportToWkb(wkbXDR,(uint8_t*)PyArray_GETPTR1(outdata,0));
    
    delete result;
    
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
    unsigned char* sdata = (unsigned char*)PyArray_DATA(shapebinary);
    OGRGeometry* shape;
    
    if(OGRGeometryFactory::createFromWkb(sdata,NULL,&shape,-1))
    {
        PyErr_SetString(PyExc_ValueError,"Error parsing WKB data.");
        return NULL;
    }
    
    OGRSpatialReference srs = OGRSpatialReference();
    
    switch(mode)
    {
        case 0://WKT
            if(srs.importFromWkt((char**)&srcWKT))
            {   
                OGRGeometryFactory::destroyGeometry(shape);
                PyErr_SetString(PyExc_ValueError,(std::string("Invalid WKT: ") + srcWKT).c_str());
                return NULL;
            }
            break;
        case 1://EPSG
            if(srs.importFromEPSG(atoi(srcWKT)))
            {
                OGRGeometryFactory::destroyGeometry(shape);
                PyErr_SetString(PyExc_ValueError,(std::string("Invalid EPSG number: ") + srcWKT).c_str());
                return NULL;
            }
            break;
        case 2://PROJ4
            if(srs.importFromProj4(srcWKT))
            {
                OGRGeometryFactory::destroyGeometry(shape);
                PyErr_SetString(PyExc_ValueError,(std::string("Invalid Proj4 String: ") + srcWKT).c_str());
                return NULL;
            }
            break;
        default:
            OGRGeometryFactory::destroyGeometry(shape);
            PyErr_SetString(PyExc_ValueError,"Mode must be one of: 0 - WKT, 1 - EPSG, 2 - Proj4");
            return NULL;
    }
    
    OGRSpatialReference dst = OGRSpatialReference();
    if(dst.importFromWkt((char**)&dstWKT))
    {
        OGRGeometryFactory::destroyGeometry(shape);
        PyErr_SetString(PyExc_ValueError,"Error determining destination spatial reference.");
        return NULL;
    }
    
    OGRCoordinateTransformation* trans = OGRCreateCoordinateTransformation(&srs,&dst);
    if(trans==NULL)
    {
        OGRGeometryFactory::destroyGeometry(shape);
        PyErr_SetString(PyExc_ValueError,"Could not compute coordinate transformation.");
        return NULL;
    }
    
    if(shape->transform(trans)) 
    {
        OGRGeometryFactory::destroyGeometry(shape);
        return NULL;
    };
    
    int dimensions[1] = {shape->WkbSize()};
    PyArrayObject* outdata = (PyArrayObject*)PyArray_FromDims(1,dimensions,PyArray_BYTE);
    
    shape->exportToWkb(wkbNDR,(byte*)PyArray_DATA(outdata));
    
    //memory cleanup
    OGRGeometryFactory::destroyGeometry(shape);
    
    //return data
    return Py_BuildValue("N",(PyObject*)outdata);
}

struct shape_state
{
    PyObject *error;
};

#define GETSTATE(m) ((struct shape_state*)PyModule_GetState(m))

static PyMethodDef shape_methods[] = 
{
    {"shapeSlice",shapeSlice, METH_VARARGS, "Get a boolean array indexing the intersection of a raster and a shape"},
    {"shapeTransform",shapeTransform,METH_VARARGS, "Translate a shape from one projection into another"},
    {"polygonize",polygonize,METH_VARARGS,"Convert boolean mask into polygon"},
    {NULL, NULL, 0, NULL}
};

static int shape_traverse(PyObject *m, visitproc visit, void *arg)
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int shape_clear(PyObject *m)
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef shapedef = {
  PyModuleDef_HEAD_INIT,
  "shape",
  NULL,
  sizeof(struct shape_state),
  shape_methods,
  NULL,
  shape_traverse,
  shape_clear,
  NULL
};

PyMODINIT_FUNC 
PyInit_shape(void)
{
    PyObject *module = PyModule_Create(&shapedef);
    import_array();
    
    if (module == NULL)
        return NULL;
    struct shape_state *st = GETSTATE(module);

    st->error = PyErr_NewException("shape.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        return NULL;
    }
    
    return module;
}
