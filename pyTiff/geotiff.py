#!/usr/bin/env python

try:
    from osgeo import gdal
except ImportError:
    import gdal
    
import numpy as np
import tempfile, os, subprocess

def nptype2gdal(nptype):
    '''Convert numpy data type to GDAL data type'''
    if(nptype == np.bool):
        return gdal.GDT_Byte
    elif(nptype == np.uint8):
        return gdal.GDT_Byte
    elif(nptype == np.uint16):
        return gdal.GDT_UInt16
    elif(nptype == np.uint32):
        return gdal.GDT_UInt32
    elif(nptype == np.int16):
        return gdal.GDT_Int16
    elif(nptype == np.int32):
        return gdal.GDT_Int32
    elif(nptype == np.float32):
        return gdal.GDT_Float32
    elif(nptype == np.float64):
        return gdal.GDT_Float64
        
def gdaltype2np(gtype):
    '''Convert GDAL data type into numpy data type'''
    if(gtype == gdal.GDT_Byte):
        return np.uint8
    elif(gtype == gdal.GDT_UInt16):
        return np.uint16
    elif(gtype == gdal.GDT_UInt32):
        return np.uint32
    elif(gtype == gdal.GDT_Int16):
        return np.int16
    elif(gtype == gdal.GDT_Int32):
        return np.int32
    elif(gtype == gdal.GDT_Float32):
        return np.float32
    elif(gtype == gdal.GDT_Float64):
        return np.float64

class geotiff:
    '''Convenience wrapper for handling GeoTIFF files using GDAL.'''
    
    def __init__(self,inputTIF,band=None):
        '''Create a geotiff instance using either a file or data.
Instances created using data (virtual geotiff) do not assume any projection information,
this must be set through the instances geoTransform and projection properties,
these properties are the same as used by gdal.'''
        
        if not (isinstance(inputTIF,str) or isinstance(inputTIF,np.ndarray)):
            raise TypeError("inputTIF should be a string or numpy array")
        
        #file based geotiff
        if isinstance(inputTIF,str):
            
            self.inputTIF = inputTIF
            
            ds = gdal.Open(self.inputTIF, gdal.GA_ReadOnly)
            
            if ds == None:
                raise ValueError("Invalid file name: "+inputTIF)
            
            #get image information
            self.width = ds.RasterXSize
            self.height = ds.RasterYSize
            self.geoTransform = ds.GetGeoTransform()
            self.projection = ds.GetProjection()
            self.gcpProjection = ds.GetGCPProjection()
            self.GCPs = ds.GetGCPs()

            self.data = None
            
            #if the user specified a band, treat the image as if it is the only one that exists
            if not band == None:

                self.bands = 1
                self.band = band
                self.nodata = [ds.GetRasterBand(self.band+1).GetNoDataValue()]

            else:#otherwise use all bands
                
                self.bands = ds.RasterCount
                self.band = None
                self.nodata = []
                
                for x in xrange(self.bands):
                    self.nodata += [ds.GetRasterBand(x+1).GetNoDataValue()]
                
            ds = None

        #virtual geotiff, does not assume any projection information -- must be set independently if known
        elif isinstance(inputTIF,np.ndarray):
            if not inputTIF.ndim in [2,3]:
                raise ValueError("InputTIF data should be 2 (for a single band) or 3 (for multi-band)")
            
            if inputTIF.ndim==3:
                self.width = inputTIF.shape[2]
                self.height = inputTIF.shape[1]
                self.bands = inputTIF.shape[0]
            else:
                self.width = inputTIF.shape[1]
                self.height = inputTIF.shape[0]
                self.bands=1
                
            self.band = None
            self.geoTransform = None
            self.projection = ""
            self.gcpProjection = ""
            self.GCPs = ()
            self.data = inputTIF
            self.inputTIF = None
            self.nodata = [None for x in xrange(self.bands)]
            
      
    def geocopy(self, outputTIF, data=None, nodata=None, options=[], create=False):
        '''Create a new GeoTIFF file that posseses the same transform and projection as this geotiff, but with new data
        
outputTIF - name of the new GeoTIFF to be created
data - numpy array of data that will be stored in the new file.
       The array should be shaped (x,h,w) where 
       x is the number of bands
       h is the height of the image
       w is the width of the image
       Height and Width should also be the same as this geotiff or the transform copy is meaningless
nodata - list of per-band nodata values that will override the images current values
options - gdal create options
       
If this is a virtual geotiff and no input data is specified, this geotiff will instead be written to
disk and this object will be converted to a real geotiff'''
        
        #check if the user supplied data, if they didn't this instance must be virtual
        if data==None and self.data == None:
            raise ValueError("You must specify data to be copied, this is not a virtual geotiff.")
        
        #if virtual, get the data that will be written out
        if data==None:
            vrt=True
            data = self.data
        else:
            vrt=False
        
        if not isinstance(outputTIF,str):
            raise TypeError("outputTIF should be a string")
        
        if not isinstance(data,np.ndarray):
            raise TypeError("data should be a numpy array")

        #check that dimensions match
        if data.ndim==3:
            bands = data.shape[0]
            if not (data.shape[1] == self.height and data.shape[2] == self.width):
                raise ValueError("image sizes must match")
        elif data.ndim==2:
            bands = 1
            if not (data.shape[0] == self.height and data.shape[1] == self.width):
                raise ValueError("image sizes must match")
        else:
            raise ValueError("data should be 2 or 3 dimensional")
        
        if bands>20:
            print "Warning: Excessive raster bands detected. Did you format your data properly?"
        
        #create a new GeoTIFF
        driver = gdal.GetDriverByName("GTiff")
        
        if(os.path.dirname(outputTIF) != '' and not os.path.isdir(os.path.dirname(outputTIF))):
            if(create):
                os.mkdir(os.path.dirname(outputTIF))
            else:
                raise IOError("Directory '"+os.path.dirname(outputTIF)+"' does not exist.")
            
        dst_ds = driver.Create(outputTIF,self.width,self.height,bands,nptype2gdal(data.dtype),options)
            
        
        if(dst_ds==None):
            raise IOError("Error creating "+outputTIF)
                
        
        #write projection information if known
        if(not (self.geoTransform==None or self.projection=="")):
            dst_ds.SetGeoTransform(self.geoTransform)
            dst_ds.SetProjection(self.projection)
        
        #write GCP information if known
        if(not (self.GCPs == () or self.gcpProjection=="")):
            dst_ds.SetGCPs(self.GCPs,self.projection)

        #copy data into the new image and set nodata values as appropriate
        if(data.ndim==3):
            for x in xrange(0,len(data)):
                dst_dr = dst_ds.GetRasterBand(x+1)
                dst_dr.WriteArray(data[x])
                
                if(nodata != None and x<len(nodata) and nodata[x]!=None):
                    dst_dr.SetNoDataValue(nodata[x])
                elif(self.inputTIF==None and x<len(self.nodata) and self.nodata[x]!=None):
                    dst_dr.SetNoDataValue(self.nodata[x])
                    
        else:
            dst_dr = dst_ds.GetRasterBand(1)
            dst_dr.WriteArray(data)
            
            if(nodata != None and nodata[0]!=None):
                    dst_dr.SetNoDataValue(nodata[0])
            elif(self.inputTIF==None and self.nodata[0]!=None):
                dst_dr.SetNoDataValue(self.nodata[0])
            
        dst_ds = None
        
        #if this was a virtual geotiff that got written to disk, convert to regular geotiff
        if self.inputTIF==None and vrt:
            self.inputTIF = outputTIF
            self.data = None
            self.nodata = [None for x in xrange(self.bands)]
            return self
            
        else:#otherwise return the new geotiff
            return geotiff(outputTIF)
            
    def geovirt(self,data,nodata=None):
        '''Returns a virtual geotiff using data and projection information from this instance'''
        
        if data.ndim==3:
            if not (data.shape[1] == self.height and data.shape[2] == self.width):
                raise ValueError("image sizes must match")
        elif data.ndim==2:
            bands = 1
            if not (data.shape[0] == self.height and data.shape[1] == self.width):
                raise ValueError("image sizes must match")
        else:
            raise ValueError("data should be 2 or 3 dimensional")
        
        g = geotiff(data)
        
        #copy any known geolocation information
        g.geoTransform = self.geoTransform
        g.projection = self.projection
        g.gcpProjection = self.gcpProjection
        g.GCPs = self.GCPs
        
        #copy nodata values if supplied
        if(nodata!=None):
            g.nodata = nodata
        else:
            g.nodata = [None for x in xrange(g.bands)]
        
        return g
    
    def getData(self,tp=(np.float32),band=None,xoff=0,yoff=0,xsize=None,ysize=None):
        '''Open a GeoTIFF and return the raw data as a numpy array
    
tp - numpy datatype of the output array
band - if the image contains multiple bands this specifies that only a single band should be returned
       a value of -1 returns an array filled with ones with the same dimensions as a single band from this raster'''

        #set band restriction if it is in place
        if not self.band == None:
            band = self.band
        
        #setting band to -1 returns an array filled with ones with the same dimensions as a single band from this raster
        if band == -1:
            return np.ones((self.height,self.width),dtype=tp)
            
        if xsize==None:
            xsize=self.width
            
        if ysize==None:
            ysize=self.height

        #if this is a virual geotiff return data
        if not self.data == None:
            if self.data.ndim==3 and not band==None:
                na  = self.data[band,yoff:yoff+ysize,xoff:xoff+xsize]
            elif self.data.ndim==3:
                na  = self.data[:,yoff:yoff+ysize,xoff:xoff+xsize]
            else:
                na = self.data[yoff:yoff+ysize,xoff:xoff+xsize]
            
            if tp==None:
                return na
            else:
                return np.array(na,dtype=tp)
            

        #otherwise read the file to get the data
        ds = gdal.Open(self.inputTIF, gdal.GA_ReadOnly)
        
        if(tp==None):
            na = ds.ReadAsArray(xoff,yoff,xsize,ysize)
        else:
            na = np.array(ds.ReadAsArray(xoff,yoff,xsize,ysize),dtype=tp)
            
        ds = None
        
        if na.ndim==3 and na.shape[0]==1:
            return na[0,:,:]
        elif na.ndim==3 and not band==None:
            return na[band,:,:]
        else:
            return na 
            
    def getPath(self):
        '''Return the system path to the image'''
        return self.inputTIF
        
    def getName(self):
        '''Return only the name of the tiff file'''
        name = self.inputTIF[::-1]
        index = name.find(os.path.sep)
        
        if(index == -1):
            return name[::-1]
            
        return name[:index][::-1]
        
    def isVirtual(self):
        '''Return True if this is a virtual geotiff'''
        return not self.data==None
        
    def nullMask(self):
        '''return a virtual geotiff with all data values set to 255'''
        na = np.ones((self.height,self.width),dtype=np.uint8)*255
        return self.geovirt(na)
        
    def intersect(self,secondTIF,outputTIF=None,nodata = None):
        '''Intersect two geotiff datasets, the result will be a geotiff object with the same
geoTransform and projection as this dataset, but will contain warped values from the second dataset

secondTIF - geotiff to warp
outputTIF - optional name for a new dataset if one is created
nodata    - if defined, should be a list object with nodata values specified for the bands in secondTIF

returns a geotiff object with the same geoTransform and projection as this object. If this is
already true of the secondTIF parameter, this method will simply return secondTIF'''

        try:
            from warpCopy import warpCopy
            d=None
            
            if not isinstance(secondTIF,geotiff):
                raise TypeError("secondTIF must be a valid geotiff object")
            
            #check to make sure we really need to warp
            if(self.geoTransform == secondTIF.geoTransform and self.projection == secondTIF.projection and self.gcpProjection == secondTIF.gcpProjection and self.GCPs == secondTIF.GCPs):
                return secondTIF
            
            #if the second tiff has no projection information but matches the dimensions of the first, assume they are overlapping
            if(secondTIF.geoTransform == (0,1,0,0,0,1) and secondTIF.projection == "" and  secondTIF.gcpProjection == "" and secondTIF.GCPs == () and self.height == secondTIF.height and self.width == secondTIF.width):
                return secondTIF
            
            if secondTIF.isVirtual():
                bands = secondTIF.bands
                if(nodata == None):
                    nodata = [None for x in xrange(bands)]
                    
                tp = secondTIF.data.dtype
                d = tempfile.mkdtemp(prefix="pytiff")
                secondTIF = secondTIF.geocopy(os.path.join(d,"temptif.tif"),secondTIF.getData(tp=None))
                
            else:
                #get default nodata values
                ds = gdal.Open(secondTIF.getPath(),gdal.GA_ReadOnly)
                bands = ds.RasterCount

                if(nodata == None):
                    nodata = [None for x in xrange(bands)]
                for x in xrange(bands):
                    val = ds.GetRasterBand(x+1).GetNoDataValue()
                    if(nodata[x] == None):
                        nodata[x] = val
                        if(nodata[x] == None):
                            print "Warning: "+secondTIF.getPath()+" band "+repr(x)+" does not have an associated nodata value."
                
                #get the data type for the new raster, assuming the first band in the second TIF is representative of all bands
                tp = gdaltype2np(ds.GetRasterBand(1).DataType)

                
                ds = None
            
            #create the new dataset
            na = np.zeros((bands,self.height,self.width),dtype=tp)
            
            tmp = self.geovirt(na)
            tmp.nodata = nodata
            
            #create temporary file if output file is not desired
            if(outputTIF==None):
                if d==None:
                    d = tempfile.mkdtemp(prefix="pytiff")
                    
                outputTIF = os.path.join(d,'rpj.tif')
            
            tmp.geocopy(outputTIF)
            
            #perform the warp
            warpCopy(secondTIF.getPath(),tmp.getPath(),nodata)
            
            if d==None:   
                return tmp
            else:
                #create virtual geotiff
                gv = tmp.geovirt(tmp.getData(tp=None),nodata=nodata)
                tmp = None
                
                #remove temporary file
                subprocess.call(['rm','-r','-f',d])
                return gv
            
        except ImportError:
            raise NotImplementedError("You must build the supplementary C++ module to enable this method.")
        
    def __getitem__(self,slc):
        return self.getData(tp=None)[slc]
        
        
    def getCoord(self,point=None):
        '''Returns the lat/lon coords of a pixel in the image, or a tuple of numpy arrays with lat/lon for every point in the image'''
        
        #geotiff must have projection information for this to work
        if(self.projection!=None and self.geoTransform!=None):
            gt = self.geoTransform
            
            if(point != None):
                x,y = point
                
                return (gt[0]+gt[1]*x+gt[2]*y, gt[3]+gt[4]*x+gt[5]*y)
            else:
                mgx,mgy = np.meshgrid(np.arange(0,self.width,1), np.arange(0,self.height,1))
                
                lon = gt[0]+gt[1]*mgx+gt[2]*mgy
                lat = gt[3]+gt[4]*mgx+gt[5]*mgy
                
                return (lon,lat)
            
    def getXY(self, lon,lat):
        '''Returns the x,y pixel coordinates of a Geolocated point in the image'''
        
        #geotiff must have projection information for this to work
        if(self.projection!=None and self.geoTransform!=None):
            gt = self.geoTransform
            
            lon = lon-gt[0]
            lat = lat-gt[3]
            
            a = np.matrix(((gt[1:3]),(gt[4:6])))
            
            i = a.getI()
            
            vals = i*[[lon],[lat]]
            
            return int(vals[0]),int(vals[1])