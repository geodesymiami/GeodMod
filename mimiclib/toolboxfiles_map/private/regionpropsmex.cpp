// Copyright 1993-2007 The MathWorks, Inc.
  
//////////////////////////////////////////////////////////////////////////////
//  Helper MEX-file for REGIONPROPS.
//  
//  Inputs:
//  prhs[0] - mxArray - label matrix
//  prhs[1] - int     - number of objects in image.
//
//  Output:
//  cell array containing a pixel index list for each region.
//////////////////////////////////////////////////////////////////////////////

#include "mex.h"
#include "iptutil_cpp.h"
#include "regionpropsmex.h"

//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
PixelIndexList::PixelIndexList()
{
    fLabel = NULL;
    fNumObj = -1;
    fIdxArray = NULL;
}

//////////////////////////////////////////////////////////////////////////////
// computeIndexLists
// Compute pixel Index Lists and store it in an output cell array.
//////////////////////////////////////////////////////////////////////////////
mxArray *PixelIndexList::computeIndexLists()
{

    mxAssert(fNumObj != -1, ERR_STRING("PixelIndexList::computeIndexLists()",
                                       "fNumObj is -1."));
    mxAssert(fLabel != NULL, ERR_STRING("PixelIndexList::computeIndexLists()",
                                       "fLabel is NULL."));

    const mxClassID labelClass = mxGetClassID(fLabel);
    void *In = mxGetData(fLabel);

    switch (labelClass)
    {
    case mxDOUBLE_CLASS:
        return(compute((double *) In));
        break;
        
    case mxSINGLE_CLASS:
        return(compute((float *) In));
        break;
        
    case mxUINT8_CLASS:
        return(compute((uint8_T *) In));
        break;
        
    case mxUINT16_CLASS:
        return(compute((uint16_T *) In));
        break;
        
    case mxUINT32_CLASS:
        return(compute((uint32_T *) In));
        break;
        
    case mxINT8_CLASS:
        return(compute((int8_T *) In));
        break;
        
    case mxINT16_CLASS:
        return(compute((int16_T *) In));
        break;
        
    case mxINT32_CLASS:
        return(compute((int32_T *) In));
        break;
        
    default:
        mexErrMsgIdAndTxt("Images:regionpropsmex:invalidType",
                          "%s", "fLabel must be numeric.");
        return(NULL);
    }
}

//////////////////////////////////////////////////////////////////////////////
// setLabelMatrix
// Set the label matrix.
//////////////////////////////////////////////////////////////////////////////
void PixelIndexList::setLabelMatrix(const mxArray *lab)
{
    fLabel = lab;
}

//////////////////////////////////////////////////////////////////////////////
// setNumberOfObjects
// Set the number of objects in the label matrix.
//////////////////////////////////////////////////////////////////////////////
void PixelIndexList::setNumberObjects(const mxArray *num)
{
    fNumObj = (int)mxGetScalar(num);
}

//////////////////////////////////////////////////////////////////////////////
// initializeIdxArray
// Initialize output cell array.
//////////////////////////////////////////////////////////////////////////////
void PixelIndexList::initializeIdxArray()
{
    mxAssert(fNumObj != -1, ERR_STRING("PixelIndexList::initializeIdxArray()",
                                       "fNumObj is -1."));
    fIdxArray = mxCreateCellMatrix(fNumObj,1);
    if (!fIdxArray)
    { 
        mexErrMsgIdAndTxt("Images:regionpropsmex:outOfMemory",
                              "%s", "Out of memory.");
    }

}


//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
PixelIndexList pixelIdxList;

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    (void) nlhs;  // unused parameter

    if (nrhs != 2)
    {
        mexErrMsgIdAndTxt("Images:regionpropsmex:invalidNumInputs",
                          "REGIONPROPSMEX requires 2 input arguments.");
    }

    if (mxIsSparse(prhs[0]) || (mxIsComplex(prhs[0])))
    {
        mexErrMsgIdAndTxt("Images:regionpropsmex:invalidInput",
                          "First input must be real and nonsparse.");
    }

    if (!(mxIsNumeric(prhs[0])))
    {
        mexErrMsgIdAndTxt("Images:regionpropsmex:invalidInput",
                          "REGIONPROPSMEX requires numeric input.");
    }        
    pixelIdxList.setLabelMatrix(prhs[0]);
    pixelIdxList.setNumberObjects(prhs[1]);
    pixelIdxList.initializeIdxArray();
    plhs[0] = pixelIdxList.computeIndexLists();
}
