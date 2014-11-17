#!/usr/bin/env python

#Copyright 2012 Reid Sawtell

#This file is part of pyGTiff.

#pyGTiff is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#pyGTiff is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with pyGTiff. If not, see <http://www.gnu.org/licenses/>.

try:
    from osgeo import gdal
except ImportError:
    import gdal
    
gdal.UseExceptions()
    
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
    elif(nptype == np.int64):
        return gdal.GDT_Int32
    elif(nptype == np.float32):
        return gdal.GDT_Float32
    elif(nptype == np.float64):
        return gdal.GDT_Float64
    elif(nptype == np.complex64):
        return gdal.GDT_CFloat32
    elif(nptype == np.complex128):
        return gdal.GDT_CFloat64
        
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
    elif(gtype == gdal.GDT_CFloat32):
        return np.complex64
    elif(gtype == gdal.GDT_CFloat64):
        return np.complex128
    
def nptype2str(nptype):
    '''Convert numpy data type to type string'''
    if(nptype == np.bool):
        return "Byte"
    elif(nptype == np.uint8):
        return "Byte"
    elif(nptype == np.uint16):
        return "UInt16"
    elif(nptype == np.uint32):
        return "UInt32"
    elif(nptype == np.int16):
        return "Int16"
    elif(nptype == np.int32):
        return "Int32"
    elif(nptype == np.int64):
        return "Int32"
    elif(nptype == np.float32):
        return "Float32"
    elif(nptype == np.float64):
        return "Float64"
    elif(nptype == np.complex64):
        return "CFloat32"
    elif(nptype == np.complex128):
        return "CFloat64"

class geotiff:
    '''Convenience wrapper for handling GeoTIFF files using GDAL.'''
    
    def __init__(self,inputTIF,band=None):
        '''Create a geotiff instance using either a file or data.
        
        Instances created using data (virtual geotiff) do not assume any projection information,
        this must be set through the instances geoTransform and projection properties,
        these properties are the same as used by gdal.
        
        Args:
            inputTIF: either a string containing the path to a geotiff file
                      or a numpy array of either 2 or 3 dimensions.
            band: optionally specify a valid band number and the tiff will be
                  treated as a single band image (2d numpy array)
                  
        Raises:
            TypeError: inputTIF is of an invalid type
            ValueError: if inputTIF is a string and does not point to a valid file.
                        if inputTIF is an array and is not either 2 or 3 dimensions.'''
        
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
           
            self.shape = (self.bands,self.height,self.width)
            ds = None

        #virtual geotiff, does not assume any projection information -- must be set independently if known
        elif isinstance(inputTIF,np.ndarray):
            if not inputTIF.ndim in [2,3]:
                raise ValueError("InputTIF data should be 2 (for a single band) or 3 (for multi-band)")
            
            if inputTIF.ndim==3:
                self.width = inputTIF.shape[2]
                self.height = inputTIF.shape[1]
                self.bands = inputTIF.shape[0]
                self.shape = (self.bands,self.height,self.width)
                self.data = inputTIF
            else:
                self.width = inputTIF.shape[1]
                self.height = inputTIF.shape[0]
                self.bands=1
                self.shape = (1, self.height,self.width)
                self.data = np.array([inputTIF])
                
            self.band = None
            self.geoTransform = None
            self.projection = ""
            self.gcpProjection = ""
            self.GCPs = ()
            
            self.inputTIF = None
            self.nodata = [None for x in xrange(self.bands)]
            
    def listFormats(self,short=False):
        nd = gdal.GetDriverCount()
        
        wdr = []
        
        for x in xrange(nd):
            drv = gdal.GetDriver(x)
            if( ('DCAP_CREATE' in drv.GetMetadata() and drv.GetMetadata()['DCAP_CREATE'] == "YES") or
                ('DCAP_CREATECOPY' in drv.GetMetadata() and drv.GetMetadata()['DCAP_CREATECOPY'] == "YES")):
                if short:
                    wdr += [drv.GetDescription()]
                else:
                    wdr += [drv.GetMetadata()['DMD_LONGNAME']+" ("+drv.GetDescription()+")"]
        return wdr
    
    def __copy_SRS(self,dst_ds,data,nodata):
            #write projection information if known
            if(not (self.geoTransform==None or self.projection=="")):
                dst_ds.SetGeoTransform(self.geoTransform)
                dst_ds.SetProjection(self.projection)
            
            #write GCP information if known
            if self.GCPs != () and self.GCPs != None:
                
                if self.gcpProjection != None and self.gcpProjection != "": #use gcpProjection if specified
                    dst_ds.SetGCPs(self.GCPs,self.gcpProjection)
                elif self.projection!=None and self.projection!="":#otherwise use projection if it is specified
                    dst_ds.SetGCPs(self.GCPs,self.projection)
                else:
                    dst_ds.SetGCPs(self.GCPs,"")#no projection available (leave it up to user to figure it out)
                    

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
                    
            return dst_ds
            
      
    def geocopy(self, outputTIF, data=None, nodata=None, options=[], create=False, compress=False, fformat="GTiff"):
        '''Create a new GeoTIFF file that posseses the same transform and projection as this geotiff, but with new data
        
        If this is a virtual geotiff and no input data is specified, this geotiff will instead be written to
        disk and this object will be converted to a real geotiff.
        
        Prints a warning if the number of raster bands exceeds 20, this does not stop you from creating a raster
        but is intended to help debugging. Most raster images have 4 or fewer bands with the exception of
        hyperspectral imagery, so if you are outputting an image in excess of 20 bands it is likely the
        input array does not have the raster bands as the first axis (ex: 640x480x3 (incorrect) instead of 3x640x480 (correct)).
            
        Args:
            outputTIF: name of the new GeoTIFF to be created.
            data: numpy array of data that will be stored in the new file.
                The array should be shaped (x,h,w) where 
                x is the number of bands
                h is the height of the image
                w is the width of the image
                Height and Width should also be the same as this geotiff or the transform copy is meaningless.
            nodata: list of per-band nodata values that will override the images current values.
            options: list gdal create options, see http://www.gdal.org/frmt_gtiff.html.
            create: set to True if the output directory should be created if it does not exist, default is False.
            compress: if True, output geotiff will use deflate compression, default is False.
            
        Returns:
            geotiff instance pointing to the newly copied file.
            
        Raises: 
            TypeError: Input arguments are of the incorrect type.
            ValueError: Image sizes do not match.
                        Data is not either 2 or 3 dimensions.
            IOError: Error occurred writing file to disk.
        '''
        
        #if virtual, get the data that will be written out
        if data==None and self.data!=None:
            vrt=True
            data = self.data
        elif data==None and self.data==None:
            vrt=False
            data = self.getData(tp=None)
        else:
            vrt=False
        
        if not isinstance(outputTIF,str):
            raise TypeError("outputTIF should be a string")
        
        if not isinstance(data,np.ndarray):
            raise TypeError("data should be a numpy array")
        
        if compress:
            options += ['COMPRESS=DEFLATE']

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
        
        #create a new file from the output format
        driver = gdal.GetDriverByName(fformat)
        
        if driver == None:
            raise ValueError('Unrecognized file format "'+fformat+'"')
        
        if not fformat in self.listFormats(short=True):
            raise IOError('"'+fformat+'" is not a writeable format')
        
        #validate output directory
        if(os.path.dirname(outputTIF) != '' and not os.path.isdir(os.path.dirname(outputTIF))):
            if(create):
                os.mkdir(os.path.dirname(outputTIF))
            else:
                raise IOError("Directory '"+os.path.dirname(outputTIF)+"' does not exist.")
            
        #use CREATE first since it's simpler
        if('DCAP_CREATE' in driver.GetMetadata() and driver.GetMetadata()['DCAP_CREATE'] == "YES"):
            
            dst_ds = driver.Create(outputTIF,self.width,self.height,bands,nptype2gdal(data.dtype),options)
                
            
            if(dst_ds==None):
                raise IOError("Error creating "+outputTIF)
                    
            
            self.__copy_SRS(dst_ds,data,nodata)
                
            dst_ds = None
            
        else: #use CREATECOPY using temporary file
            
            if not nptype2str(data.dtype) in driver.GetMetadata()['DMD_CREATIONDATATYPES']:
                raise IOError(nptype2str(data.dtype)+" is not supported by "+fformat+".\nSupported types are: "+driver.GetMetadata()['DMD_CREATIONDATATYPES'])

           
            dt = tempfile.mkdtemp(prefix="pytiff")
            temptif = self.geocopy(os.path.join(dt,"temptif.tif"),data)
            
            src_ds = gdal.Open(temptif.inputTIF, gdal.GA_ReadOnly)
            dst_ds = driver.CreateCopy(outputTIF,src_ds,1)
            
            if dst_ds == None:
                raise IOError("Failed to create output file")
            
            src_ds = None
            dst_ds = None
            
            subprocess.call(['rm','-r','-f',dt])
            
        #if this was a virtual geotiff that got written to disk, convert to regular geotiff
        if self.inputTIF==None and vrt:
            self.inputTIF = outputTIF
            self.data = None
            self.nodata = [None for x in xrange(self.bands)]
            return self
            
        else:#otherwise return the new geotiff
            return geotiff(outputTIF)
            
    def geovirt(self,data,nodata=None):
        '''Creates a virtual geotiff using data and projection information from a geotiff object
        
        Args:
            data: numpy array of data that will be stored in the new file.
                The array should be shaped (x,h,w) where 
                x is the number of bands
                h is the height of the image
                w is the width of the image
                Height and Width should also be the same as this geotiff or the transform copy is meaningless.
            nodata: list of per-band nodata values that will override the images current values
            
        Returns:
            A virtual geotiff object, the distinction being the data is stored purely in memory
            and is not yet written to disk.
            
        Raises: 
            TypeError: Input arguments are of the incorrect type.
            ValueError: Image sizes do not match.
                        Data is not either 2 or 3 dimensions.
       '''
        
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
        
        #copy nodata values if supplied and number of bands matches
        if(nodata!=None):
            g.nodata = nodata
        else:
            if len(g.nodata)<g.bands:
                g.nodata = [None]*g.bands
            else:
                g.nodata = self.nodata
        
        
        #'fix' nodata values in the data array to match the specified nodata
        if isinstance(data,np.ma.masked_array):
            if data.ndim==3:
                for b in xrange(data.shape[0]):
                    if g.nodata[b]!=None:
                        g.data[b][data[b].mask!=0]==g.nodata[b]
            elif data.ndim==2:
                if g.nodata[0]!=None:
                    g.data[0,data.mask!=0]==g.nodata[0]
        
        return g
    
    def getData(self,tp=np.float32,band=None,xoff=0,yoff=0,xsize=None,ysize=None):
        '''Open a GeoTIFF and return the raw data as a numpy array
        
        This method does not support dynamic updating, therefor every time this method is called the entirety of the file is read.
        Best practice is to use slicing notation to retrieve the largest subset of the file needed for your program, save the return
        array and reuse it to avoid additional read operations.
        
        Args:
    
            tp: numpy datatype of the output array, set to None to use the native datatype, defaults to singe precision float
            band: If the image contains multiple bands this specifies that only a single band should be returned.
                A value of -1 returns an array filled with ones with the same dimensions as a single band from this raster.
                
        Returns:
            Numpy array containing data read in from file.'''

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
            
        if band==None:
            band = range(0,self.bands)
            
        if isinstance(band,list) and len(band)==1:
            band = band[0]

        #if this is a virual geotiff return data
        if not self.data == None:
            if self.data.ndim==3 and not band==None:
                na  = self.data[band,yoff:yoff+ysize,xoff:xoff+xsize]
            elif self.data.ndim==3:
                na  = self.data[:,yoff:yoff+ysize,xoff:xoff+xsize]
            else:
                na = self.data[yoff:yoff+ysize,xoff:xoff+xsize]
            
            if tp==None:
                return self.__ndma__(na,band=band)
            else:
                return self.__ndma__(np.array(na,dtype=tp),band=band)
            

        #otherwise read the file to get the data
        ds = gdal.Open(self.inputTIF, gdal.GA_ReadOnly)
    
        if isinstance(band,int):
            na = ds.GetRasterBand(band+1).ReadAsArray(xoff,yoff,xsize,ysize)
        else:
            if(tp==None):
                na = np.array([ds.GetRasterBand(x+1).ReadAsArray(xoff,yoff,xsize,ysize) for x in band])
            else:
                na = np.array([ds.GetRasterBand(x+1).ReadAsArray(xoff,yoff,xsize,ysize) for x in band], dtype=tp)
            
        ds = None
        
        if na.ndim==3 and na.shape[0]==1:
            return self.__ndma__(na[0],band=band)
        else:
            return self.__ndma__(na,band=band)
            
    def __ndma__(self,data,band=None):
        '''Converts to numpy masked array
        
        If nodata values are specified the masked array will reflect that by masking those values.
        Note that masked arrays behave differently than normal numpy arrays for some operations,
        because of this masked arrays are returned regardless of nodata values for consistent behavior.'''
        
        if self.nodata==None:
            return np.ma.masked_array(data)

        nd = np.array([[x] for x in self.nodata])
        
        if band==None or self.band==band:
            band=np.array(np.linspace(0,nd.shape[0]-1,nd.shape[0]),dtype=np.int32)
        else:
            if isinstance(band,int):
                band = [band]
            band = np.array(band)

        if band.shape[0]>1:
            
            nd = nd[band]
            
            mask = np.zeros(data.shape,dtype=np.bool)
            
            for b in xrange(band.shape[0]):
                mask[b] = data[b]==nd[b]
            
            return np.ma.masked_array(data,mask=mask)
        else:
            return np.ma.masked_array(data,mask=data==nd[band])
        
            
    def __getitem__(self,slc):
        """Support for slicing notation
        
        Returns a numpy array containing data from the tiff file. Modifying the contents
        of that array will not affect the contents of the file itself."""

        def processList(lst,mx):
                slc2 = []
                
                for x in lst:
                    if x<0:
                        slc2 += [mx+x]
                    else:
                        slc2 += [x]
                        
                s2min = min(slc2)
                s2max = max(slc2)
                        
                if(s2min<0):
                    raise IndexError("index {0} is out of bounds for axis 0 with size {1}".format(s2min-mx,mx))
                elif(s2max>=mx):
                    raise IndexError("index {0} is out of bounds for axis 0 with size {1}".format(s2max,mx))
                
                rng = (s2min,s2max-s2min+1)
                slc = [x-s2min for x in slc2]
                
                return (rng,slc)
            
        def processSlice(slc,mx):
            
            slc2 = list(slc.indices(mx))
            #print '>>>',slc2
            
            if(slc.step>0 or slc.step==None):
                rng = (slc2[0],slc2[1]-slc2[0])
            else:
                rng = (slc2[1]+1,slc2[0]-slc2[1])
                
            #print rng
                
            slcnew = slice(None,None,slc.step)
            
            return (rng,slcnew)
        
        def process(slc,mx):
            if isinstance(slc, list):
                return processList(slc,mx)
            elif isinstance(slc,slice):
                return processSlice(slc,mx)
            elif isinstance(slc,np.ndarray):
                return (slice(None,None,None),slc)

        dim1 = slice(None,None,None)

        dim2 = (0,None)
        dim2slc = slice(None,None,None)

        dim3 = (0,None)
        dim3slc = slice(None,None,None)
        
        #only have an integer value
        if isinstance(slc, int):
            slc = ([slc],)

        #only have a list
        elif isinstance(slc, list) or isinstance(slc,slice):
            slc = (slc,)
            
        elif isinstance(slc,np.ndarray):
            return self.getData(tp=None)[slc]
        
        #slice over multiple dimensions
        if isinstance(slc,tuple):
            if len(slc)>(3-(self.band!=None or self.bands==1)):
                raise IndexError("too many indices")
                 
            for x,n in enumerate(slc):
                if isinstance(n,int):
                    n = [n]
                    
                if (self.band!=None or self.bands==1):
                    x = x+1
                
                if x==0:
                    dim1 = n
                elif x==1:
                    dim2,dim2slc = process(n,self.height)
                elif x==2:
                    dim3,dim3slc = process(n,self.width)
        
        
        
        if isinstance(dim1, slice):
            slc2 = dim1.indices(self.bands)
            dim1 = range(*slc2)
        
        data = self.getData(tp=None,band=dim1,yoff=dim2[0],ysize=dim2[1],xoff=dim3[0],xsize=dim3[1])
        
        if(data.ndim==2):
            return data[dim2slc,dim3slc]
        else:
            return data[:,dim2slc,dim2slc]
    
    def __len__(self):
        return self.bands
            
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
        '''Return a single band virtual geotiff with all data values set to 255'''
        na = np.ones((self.height,self.width),dtype=np.uint8)*255
        return self.geovirt(na)
        
    def intersect(self,secondTIF,outputTIF=None,nodata = None,resampleType=0):
        '''Intersect two geotiff datasets.
        The result will be a geotiff object with the same
        geoTransform and projection as this dataset, but will contain warped values from the second dataset.
        
        Args:
            secondTIF: geotiff to warp
            outputTIF: optional name for a new dataset if one is created
            nodata: if defined, should be a list object with nodata values specified for the bands in secondTIF
            resampleType: which resampling method should be used: one of 
                          {
                              GRA_NearestNeighbour = 0, GRA_Bilinear = 1, GRA_Cubic = 2, GRA_CubicSpline = 3,
                              GRA_Lanczos = 4, GRA_Average = 5, GRA_Mode = 6
                          }
                          defaults to 0.
                
            
        Returns:
            A geotiff object with the same geoTransform and projection as this object, with data from secondTIF.
            If secondTIF does not need to be warped then it will be returned unmodified.
            If secondTIF has no projection information but matches the image size exactly, it will be assumed to be overlapping.
            
        Raises:
            TypeError: if secondTIF is not a geotiff object
            NotImplementedError: if the C++ libraries are missing or broken'''

        try:
            from warpCopy import warpCopy
            d=None
            dt=None
            
            if not isinstance(secondTIF,geotiff):
                raise TypeError("secondTIF must be a valid geotiff object")
            
            #check to make sure we really need to warp
            if(self.geoTransform == secondTIF.geoTransform and self.projection == secondTIF.projection and self.gcpProjection == secondTIF.gcpProjection and self.GCPs == secondTIF.GCPs and self.height == secondTIF.height and self.width == secondTIF.width):
                return secondTIF
            
            #if the second tiff has no projection information but matches the dimensions of the first, assume they are overlapping
            if(secondTIF.geoTransform == (0,1,0,0,0,1) and secondTIF.projection == "" and  secondTIF.gcpProjection == "" and secondTIF.GCPs == () and self.height == secondTIF.height and self.width == secondTIF.width):
                return secondTIF
            
            if secondTIF.isVirtual():
                bands = secondTIF.bands
                if(nodata == None):
                    nodata = [None for x in xrange(bands)]
                    
                tp = secondTIF.data.dtype
                dt = tempfile.mkdtemp(prefix="pytiff")
                secondTIF = secondTIF.geocopy(os.path.join(dt,"temptif.tif"),secondTIF.getData(tp=None))
                
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
                d = tempfile.mkdtemp(prefix="pytiff")
                    
                outputTIF = os.path.join(d,'rpj.tif')
            
            tmp.geocopy(outputTIF)
            
            #perform the warp
            warpCopy(secondTIF.getPath(),tmp.getPath(),nodata,resampleType)
            
            #remove temp directory if secondTIF was virtual
            if dt!=None:
                subprocess.call(['rm','-r','-f',dt])
            
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
        
    def reproject(self,srcSRS,outputTIF=None,nodata = None,resampleType=0):
        '''Create a new reprojected geotiff according to the supplied spatial reference string
        
        Args:
            srcSRS: optionally specify the EPSG number or proj4 string for the shapes projection to have it converted
                    to match the geotiff's projection (EX: 4326, "EPSG:4326", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs").
            outputTIF: optional name for a new dataset if one is created
            nodata: if defined, should be a list object with nodata values specified for the bands in secondTIF
            resampleType: which resampling method should be used: one of 
                          {
                              GRA_NearestNeighbour = 0, GRA_Bilinear = 1, GRA_Cubic = 2, GRA_CubicSpline = 3,
                              GRA_Lanczos = 4, GRA_Average = 5, GRA_Mode = 6
                          }
                          defaults to 0.
                    
        Returns:
            A geotiff instance of the newly reprojected image
'''

        try:
            from warpCopy import warp
            d = None
            dt = None
            
            #create temporary file if virtual
            if self.isVirtual():
                bands = self.bands
                    
                tp = self.data.dtype
                dt = tempfile.mkdtemp(prefix="pytiff")
                tmp = self.geocopy(os.path.join(dt,"temptif.tif"),self.getData(tp=None))
            else:
                tmp = self
            
            if(nodata == None):
                if(tmp.nodata == None):
                    nodata = [None for x in xrange(bands)]
                else:
                    nodata = tmp.nodata
            
            for x,nd in enumerate(nodata):
                if(nd == None):
                    print "Warning: "+tmp.getPath()+" band "+repr(x)+" does not have an associated nodata value."
                    
            #create temporary file if output file is not desired
            if(outputTIF==None):
                d = tempfile.mkdtemp(prefix="pytiff")
                    
                outputTIF = os.path.join(d,'rpj.tif')
                
            #determine type of input srs    
            mode=0
            
            if(type(srcSRS)==int):
                mode=1
                srcSRS = repr(srcSRS)
            elif(srcSRS.startswith('EPSG:') or srcSRS.startswith('epsg:') or srcSRS.startswith("Epsg:")):
                srcSRS = srcSRS[5:]
                mode=1

            elif(srcSRS.startswith('+')):
                mode=2
           
            warp(tmp.getPath(),outputTIF,nodata,resampleType,srcSRS,mode)
            
            #remove temp directory if self is virtual
            if dt!=None:
                subprocess.call(['rm','-r','-f',dt])
                
            if d==None:   
                return geotiff(outputTIF)
            else:
                #create virtual geotiff
                g = geotiff(outputTIF)
                gv = g.geovirt(g.getData(tp=None),nodata=nodata)
                g = None
                
                #remove temporary file
                subprocess.call(['rm','-r','-f',d])
                return gv
            
        except ImportError:
            raise NotImplementedError("You must build the supplementary C++ module to enable this method.")
        
    def shapeIntersect(self,wkb,srcSRS=None):
        '''Rasterize a polygon.
        
        Generate an array containing the proportion of a pixel that overlaps with a polygon.
        Multiplying the return value by the area per pixel will give you the area of overlap.
        Performing a boolean operation, eg. (retVal>.5), will give you an array that can be used to index the
        geotiff data, such as for zonal statistics.
        
        Args:
            wkb: well known binary, should be in string format (For example: the shapely .wkb attribute)
                 the shape should be valid and in the same projection as the geotiff.
            srcSRS: optionally specify the EPSG number or proj4 string for the shapes projection to have it converted
                    to match the geotiff's projection (EX: 4326, "EPSG:4326", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs").
                    
        Returns:
            Numpy array with the same height,width as the tif image.
            Each value in the numpy array will contain a floating point number indicating the proportion of the
            pixel that overlaps the input shape. (Multiply by 100 to get the % overlap, so 1 is 100% overlap)
            
        Raises:
            NotImplementedError: if the C++ libraries are missing or broken
'''

        try:
            from shape import shapeSlice,shapeTransform
            
            if(self.geoTransform == None):
                gt = (0,1,0,0,0,1)
            else:
                gt = self.geoTransform
                
            sbytes = np.fromstring(wkb,dtype=np.uint8);

            #reproject shape if possible
            if self.projection!=None and srcSRS!=None:
                mode=0
                
                if(type(srcSRS)==int):
                    mode=1
                    srcSRS = repr(srcSRS)
                elif(srcSRS.startswith('EPSG:') or srcSRS.startswith('epsg:') or srcSRS.startswith("Epsg:")):
                    srcSRS = srcSRS[5:]
                    mode=1

                elif(srcSRS.startswith('+')):
                    mode=2
                    
                sbytes = shapeTransform(sbytes,srcSRS,self.projection,mode)
                
            return shapeSlice(np.array(gt,dtype=np.float32),sbytes,self.width,self.height)
            
        except ImportError:
            raise NotImplementedError("You must build the supplementary C++ module to enable this method.")
        
    def getPixelArea(self):
        '''Returns the area of the pixel in projected units (or 1 if unprojected)'''
        
        if self.geoTransform!=None:
            gt = self.geoTransform
            
            b1 = ((gt[1]+gt[2])**2 + (gt[4]+gt[5])**2)**.5
            b2 = ((gt[1]-gt[2])**2 + (gt[4]-gt[5])**2)**.5
            
            return .5*b1*b2
        else:
            return 1
        
    def getCoord(self,point=None):
        '''Get coordinates from pixel points
        
        Args: tuple containing (x,y) pixel coordinate, default is None
        
        Returns:
            If point is None this returns a tuple containing the upper-left (longitude, latitute) in projected units for every pixel in the image.
            If point is specified, a tuple containing the upper-left (longitude, latitude) for the specified point.
            Returns None if this geotiff has no projection information.'''
        
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
        '''Returns the x,y pixel coordinates of a Geolocated point in the image
        
        Args:
            lon - longitude of the x,y point in projected units
            lat - latitude of the x,y point in projected units
            
        Returns:
            A tuple containing the (x,y) pixel coordinate'''
        
        #geotiff must have projection information for this to work
        if(self.projection!=None and self.geoTransform!=None):
            gt = self.geoTransform
            
            lon = lon-gt[0]
            lat = lat-gt[3]
            
            a = np.matrix(((gt[1:3]),(gt[4:6])))
            
            i = a.getI()
            
            vals = i*[[lon],[lat]]
            
            return int(vals[0]),int(vals[1])
        
    def getBounds(self):
        '''Returns the (N,W,S,E) bounding box of the image'''
        
        ullon, ullat = self.getCoord((0,0))
        urlon, urlat = self.getCoord((self.width,0))
        lllon, lllat = self.getCoord((0,self.height))
        lrlon, lrlat = self.getCoord((self.width,self.height))
        
        lats = [ullat,urlat,lllat,lrlat]
        lons = [ullon,urlon,lllon,lrlon]
        
        return (max(lats),min(lons),min(lats),max(lons))
        