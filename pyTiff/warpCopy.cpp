#include <Python.h>
#include <math.h>
#include "gdalwarper.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

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
    
    if (!PyArg_ParseTuple(args, "ssO!", &inputName,&outputName,&PyList_Type,&nodata))
    {
        return NULL;
    }
    
    if (inputName==NULL || outputName==NULL)
    {
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

static PyMethodDef warpCopyMethods[] = 
{
    {"warpCopy",warpCopy, METH_VARARGS, "reproject an image into another existing image"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initwarpCopy(void)
{
    (void) Py_InitModule("warpCopy", warpCopyMethods);
}