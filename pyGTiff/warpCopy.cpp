#include <Python.h>
#include <math.h>
#include "gdalwarper.h"
#include "ogr_spatialref.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

#ifdef WIN32
    #ifndef NAN
        static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
        #define NAN (*(const float *) __nan)
    #endif
#endif

/**
 * Warp one image into another
 * Copied mostly from the GDAL api tutorials 
 * http://www.gdal.org/warptut.html
 */
static PyObject * warpCopy(PyObject *self, PyObject *args)
{
    const char* inputName;
    const char* outputName;
    PyObject * nodata;
    PyObject * realObject;
    double value;
    int resampleType;
    
    if (!PyArg_ParseTuple(args, "ssO!i", &inputName,&outputName,&PyList_Type,&nodata,&resampleType))
    {
        return NULL;
    }
    
    if (inputName==NULL || outputName==NULL)
    {
        PyErr_SetString(PyExc_ValueError,"Invalid input arguments.");
        return NULL;
    }
    
    GDALDatasetH  hSrcDS, hDstDS;

    // Open input and output files. 

    GDALAllRegister();

    hSrcDS = GDALOpen( inputName, GA_ReadOnly );
    hDstDS = GDALOpen( outputName, GA_Update );

    // Setup warp options. 
    
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
    char** papszOptions = NULL;
    
    psWarpOptions->eResampleAlg = (GDALResampleAlg)resampleType;
    
    papszOptions = CSLSetNameValue(papszOptions,"INIT_DEST","NODATA");
    psWarpOptions->papszWarpOptions = papszOptions;

    psWarpOptions->hSrcDS = hSrcDS;
    psWarpOptions->hDstDS = hDstDS;

    psWarpOptions->nBandCount = PyList_Size(nodata);
    psWarpOptions->panSrcBands = 
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    psWarpOptions->panDstBands = 
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    
    //add nodata values
    psWarpOptions->padfSrcNoDataReal = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    psWarpOptions->padfSrcNoDataImag = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    psWarpOptions->padfDstNoDataReal = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    psWarpOptions->padfDstNoDataImag = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    
    for(int i=0;i<PyList_Size(nodata);i++)
    {
        psWarpOptions->panSrcBands[i] = i+1;
        psWarpOptions->panDstBands[i] = i+1;
        realObject = PyList_GetItem(nodata,i);
        value = PyFloat_AsDouble(realObject);
        
        if(PyErr_Occurred()==NULL)
        {
            psWarpOptions->padfSrcNoDataReal[i] = value;
            psWarpOptions->padfDstNoDataReal[i] = value;
        }
        else
        {
            PyErr_Clear();
            psWarpOptions->padfSrcNoDataReal[i] = NAN;
            psWarpOptions->padfDstNoDataReal[i] = NAN;
        }
        
        psWarpOptions->padfSrcNoDataImag[i] = NAN;
        psWarpOptions->padfDstNoDataImag[i] = NAN;
    }

    psWarpOptions->pfnProgress = GDALTermProgress;   

    // Establish reprojection transformer. 

    psWarpOptions->pTransformerArg = 
        GDALCreateGenImgProjTransformer( hSrcDS, 
                                         GDALGetProjectionRef(hSrcDS), 
                                         hDstDS,
                                         GDALGetProjectionRef(hDstDS), 
                                         FALSE, 0.0, 1 );
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

    // Initialize and execute the warp operation. 

    GDALWarpOperation oOperation;

    oOperation.Initialize( psWarpOptions );
    oOperation.ChunkAndWarpImage( 0, 0, 
                                  GDALGetRasterXSize( hDstDS ), 
                                  GDALGetRasterYSize( hDstDS ) );

    GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
    GDALDestroyWarpOptions( psWarpOptions );

    GDALClose( hDstDS );
    GDALClose( hSrcDS );

    return Py_BuildValue("");
}

/**
 * Reproject an image
**/
static PyObject * warp(PyObject *self, PyObject *args)
{
    const char* inputName;
    const char* outputName;
    PyObject * nodata;
    PyObject * realObject;
    double value;
    int resampleType;
    int mode;
    const char* srcWKT;
    
    GDALDataType eDT;
    GDALDriverH hDriver;
    const char *pszSrcWKT = NULL;
    char *pszDstWKT = NULL;
    OGRSpatialReference srs = OGRSpatialReference();
    void* hTransformArg;
    
    if (!PyArg_ParseTuple(args, "ssO!isi", &inputName,&outputName,&PyList_Type,&nodata,&resampleType,&srcWKT,&mode))
    {
        return NULL;
    }
    
    if (inputName==NULL || outputName==NULL)
    {
        PyErr_SetString(PyExc_ValueError,"Invalid input arguments.");
        return NULL;
    }
    
    GDALDatasetH  hSrcDS, hDstDS;

    // Open input and output files. 

    GDALAllRegister();

    hSrcDS = GDALOpen( inputName, GA_ReadOnly );
    
    eDT = GDALGetRasterDataType(GDALGetRasterBand(hSrcDS,1));
    
    hDriver = GDALGetDriverByName( "GTiff" );
    if(hDriver == NULL){PyErr_SetString(PyExc_ValueError,"Failed to get GTiff driver.");return NULL;}
    
    //get input projection
    pszSrcWKT = GDALGetProjectionRef( hSrcDS );
    if(pszSrcWKT == NULL || strlen(pszSrcWKT)==0){PyErr_SetString(PyExc_ValueError,"Could not determine source image projection."); return NULL;}
    
    
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
    
    srs.exportToWkt( &pszDstWKT);
    
    hTransformArg = 
    GDALCreateGenImgProjTransformer( hSrcDS, pszSrcWKT, NULL, pszDstWKT, 
                                        FALSE, 0, 1 );
    if(hTransformArg == NULL){PyErr_SetString(PyExc_ValueError,"Failed to generate image transformation.");return NULL;}
    
    double adfDstGeoTransform[6];
    int nPixels=0, nLines=0;
    CPLErr eErr;

    eErr = GDALSuggestedWarpOutput( hSrcDS, 
                                    GDALGenImgProjTransform, hTransformArg, 
                                    adfDstGeoTransform, &nPixels, &nLines );
    if(eErr != CE_None){PyErr_SetString(PyExc_ValueError,"Failed to generate output warp.");return NULL;}

    GDALDestroyGenImgProjTransformer( hTransformArg );

    // Create the output file.  
    hDstDS = GDALCreate( hDriver, outputName, nPixels, nLines, 
                         GDALGetRasterCount(hSrcDS), eDT, NULL );
    
    if(hDstDS == NULL){PyErr_SetString(PyExc_ValueError,(std::string("Failed to create output image: ")+outputName).c_str());return NULL;}
    
    // Write out the projection definition. 
    GDALSetProjection( hDstDS, pszDstWKT );
    GDALSetGeoTransform( hDstDS, adfDstGeoTransform );
    
    //color table check
    GDALColorTableH hCT;

    hCT = GDALGetRasterColorTable( GDALGetRasterBand(hSrcDS,1) );
    if( hCT != NULL )
        GDALSetRasterColorTable( GDALGetRasterBand(hDstDS,1), hCT );
    

    // Setup warp options. 
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
    char** papszOptions = NULL;
    
    psWarpOptions->eResampleAlg = (GDALResampleAlg)resampleType;
    
    papszOptions = CSLSetNameValue(papszOptions,"INIT_DEST","NODATA");
    psWarpOptions->papszWarpOptions = papszOptions;

    psWarpOptions->hSrcDS = hSrcDS;
    psWarpOptions->hDstDS = hDstDS;

    psWarpOptions->nBandCount = PyList_Size(nodata);
    psWarpOptions->panSrcBands = 
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    psWarpOptions->panDstBands = 
        (int *) CPLMalloc(sizeof(int) * psWarpOptions->nBandCount );
    
    //add nodata values
    psWarpOptions->padfSrcNoDataReal = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    psWarpOptions->padfSrcNoDataImag = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    psWarpOptions->padfDstNoDataReal = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    psWarpOptions->padfDstNoDataImag = (double *) CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    
    for(int i=0;i<PyList_Size(nodata);i++)
    {
        psWarpOptions->panSrcBands[i] = i+1;
        psWarpOptions->panDstBands[i] = i+1;
        realObject = PyList_GetItem(nodata,i);
        value = PyFloat_AsDouble(realObject);
        
        if(PyErr_Occurred()==NULL)
        {
            psWarpOptions->padfSrcNoDataReal[i] = value;
            psWarpOptions->padfDstNoDataReal[i] = value;
        }
        else
        {
            PyErr_Clear();
            psWarpOptions->padfSrcNoDataReal[i] = NAN;
            psWarpOptions->padfDstNoDataReal[i] = NAN;
        }
        
        psWarpOptions->padfSrcNoDataImag[i] = NAN;
        psWarpOptions->padfDstNoDataImag[i] = NAN;
    }

    psWarpOptions->pfnProgress = GDALTermProgress;   

    // Establish reprojection transformer. 
    psWarpOptions->pTransformerArg = 
        GDALCreateGenImgProjTransformer( hSrcDS, 
                                         GDALGetProjectionRef(hSrcDS), 
                                         hDstDS,
                                         GDALGetProjectionRef(hDstDS), 
                                         FALSE, 0.0, 1 );
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

    // Initialize and execute the warp operation. 
    GDALWarpOperation oOperation;

    oOperation.Initialize( psWarpOptions );
    oOperation.ChunkAndWarpImage( 0, 0, 
                                  GDALGetRasterXSize( hDstDS ), 
                                  GDALGetRasterYSize( hDstDS ) );

    GDALDestroyGenImgProjTransformer( psWarpOptions->pTransformerArg );
    GDALDestroyWarpOptions( psWarpOptions );

    GDALClose( hDstDS );
    GDALClose( hSrcDS );
    
    CPLFree(pszDstWKT);

    return Py_BuildValue("");
}

static PyMethodDef warpCopyMethods[] = 
{
    {"warpCopy",warpCopy, METH_VARARGS, "reproject an image into another existing image"},
    {"warp",warp, METH_VARARGS, "reproject an image by specifying an output SRS"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initwarpCopy(void)
{
    (void) Py_InitModule("warpCopy", warpCopyMethods);
}