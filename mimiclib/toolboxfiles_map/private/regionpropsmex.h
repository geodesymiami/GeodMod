// $Revision: 1.1.8.2 $
// Copyright 1993-2007 The MathWorks, Inc.

#ifndef _REGIONPROPSMEX_H
#define _REGIONPROPSMEX_H

#include "mex.h"
#include "sequencemex.h"
#include "iptutil_cpp.h"

class PixelIndexList
{
 private:
    const mxArray *fLabel;
    mwSignedIndex  fNumObj;
    mxArray       *fIdxArray;

//////////////////////////////////////////////////////////////////////////////
// computes a list of pixel indices for each object and stores it in an output
// cell array.
//////////////////////////////////////////////////////////////////////////////
template< typename _T>
mxArray *compute(_T *lMatrixPtr)
{
    mxAssert(fNumObj != -1, ERR_STRING("PixelIndexList::compute()",
                                       "fNumObj is -1."));
    mxAssert(fLabel != NULL, ERR_STRING("PixelIndexList::compute()",
                                       "fLabel is NULL."));

    Sequence<double>  *seq;               //array of sequences
    mwSignedIndex     label;              //pixel label
    bool              *isInitialized;     //flag if sequence is initialized
    mwSize            hint = 50;          //initial sequence length
    mwSize            oneSeqLength = 0;   //length of one sequence
    mxArray           *oneIdxArray;       //double matrix w/ pixel idx list.
 

    // Set isInitialized flag for each sequence to false.
    isInitialized = (bool *) mxCalloc(fNumObj,sizeof(bool));
    for (mwSignedIndex j = 0; j < fNumObj; j++)
    {
        isInitialized[j] = false;
    }

    // Create sequences.
    seq = (Sequence<double> *) mxCalloc(fNumObj,sizeof(Sequence<double>));
    mwSize numElements = mxGetNumberOfElements(fLabel);
    double val;
    for (mwSize p = 0; p < numElements; p++)
    {
        label = (mwSignedIndex) lMatrixPtr[p] - 1;  // region with label 1 in M 
                                                    // corresponds to the 0th sequence

        if (label != -1)
        {
            if (!isInitialized[label])
            {
                seq[label].initialize(hint);
                isInitialized[label] = true;
            }

            val = (double)p + 1;          // 0 idx in C++ corresponds to 1 
                                          // in M
            seq[label].addHigh(val);
        }
    }

    // Free isInitialized because we have finished using it.
    mxFree(isInitialized);

    // Populate cell array (fIdxArray) by pointing to data pointed to by
    // oneIdxArray.
    for (mwSignedIndex m = 0; m < fNumObj; m++)
    {
        // Get length of mth sequence.
        oneSeqLength = seq[m].getSequenceLength();

        // Create oneIdxArray that will contain the sequence's data.
        oneIdxArray = mxCreateDoubleMatrix(oneSeqLength,1,mxREAL);    
        if (!oneIdxArray)
        { 
            mexErrMsgIdAndTxt("Images:regionpropsmex:outOfMemory",
                          "%s", "Out of memory.");
        }
        double *a = (double *)mxGetData(oneIdxArray);

        // Copy each element in the sequence to oneIdxArray.
        seq[m].copyToBuf(a);
        seq[m].freeSequence();

        // Have the mth cell point to the data pointed to by oneIdxArray.
        mxSetCell(fIdxArray,m,oneIdxArray);
    }

    //clean up.
    mxFree(seq);

    return(fIdxArray);
}

 public:
    PixelIndexList();
    mxArray *computeIndexLists(void);
    void setLabelMatrix(const mxArray *labelMatrix);
    void setNumberObjects(const mxArray *numObjects);
    void initializeIdxArray(void);
};

#endif
