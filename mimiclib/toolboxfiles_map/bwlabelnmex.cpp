/*
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.1.6.4 $
 */

#include "unionfind.h"
#include "neighborhood.h"
#include "iptutil.h"
#include "mex.h"

#define DEFAULT_ALLOCATED_LENGTH 16
#define REALLOC_FACTOR 2

mwSize connected_components(mxLogical *BW, double *L, mwSize num_elements,
                          NeighborhoodWalker_T walker)
{
    mwSize p;
    mwSize q;
    mwSize next_label = 0;
    mwSize num_sets;
    ufType *uf;
    bool found;
    
    uf = uf_init(1000);
    
    for (p = 0; p < num_elements; p++)
    {
        if (BW[p] != 0)
        {
            /*
             * Find the first trailing nonzero neighbor.
             */
            found = false;
            nhSetWalkerLocation(walker, p);
            while (nhGetNextInboundsNeighbor(walker, &q, NULL))
            {
                if (BW[q] != 0)
                {
                    /*
                     * Assign L[p] the same label as L[q]
                     */
                    found = true;
                    L[p] = L[q];
                    break;
                }
            }

            if (! found)
            {
                /*
                 * There were no trailing nonzero neighbors.  Assign the next label.
                 */

                L[p] = (double) next_label;

                /*
                 * Increment next_label and make space for it in the unionfind structure.
                 */
                next_label++;
                uf_new_node(uf);
            }
            else
            {
                /*
                 * Look for additional trailing nonzero neighbors.  If any
                 * exist that have a different label than what has just
                 * been assigned to L[p], record the equivalence by calling
                 * uf_new_pair().
                 */
                while (nhGetNextInboundsNeighbor(walker, &q, NULL))
                {
                    if ((BW[q] != 0) && (L[q] != L[p]))
                    {
                        uf_new_pair(uf, (mwSize) L[p], (mwSize) L[q]);
                    }
                }
            }
        }
    }
    
    num_sets = uf_renumber(uf, 1);
    
    for (p = 0; p < num_elements; p++)
    {
        if (BW[p] != 0)
        {
            L[p] = (double) uf_query_set(uf, (mwSize) L[p]);
        }
    }
    
    uf_destroy(uf);

    return num_sets;
}


extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize num_elements;
    mwSize num_dims;
    mwSize num_sets;
    const mwSize *size;
    double *L;
    mxLogical *BW;
    Neighborhood_T nhood;
    NeighborhoodWalker_T walker;

    if (nrhs < 1)
    {
        mexErrMsgIdAndTxt("Images:bwlabelnmex:tooFewInputs",
                          "BWLABELNMEX needs at least one input.");
    }
    if (nrhs > 2)
    {
        mexErrMsgIdAndTxt("Images:bwlabelnmex:tooManyInputs",
                          "BWLABELNMEX cannot take more than two inputs.");
    }

    num_elements = mxGetNumberOfElements(prhs[0]);
    num_dims = mxGetNumberOfDimensions(prhs[0]);
    size = mxGetDimensions(prhs[0]);

    plhs[0] = mxCreateNumericArray(num_dims, size, mxDOUBLE_CLASS, mxREAL);
    L = (double *) mxGetData(plhs[0]);

    BW = mxGetLogicals(prhs[0]);
    
    if (nrhs < 2)
    {
        nhood = nhMakeDefaultConnectivityNeighborhood(num_dims);
    }
    else
    {
        nhood = nhMakeNeighborhood(prhs[1],NH_CENTER_MIDDLE_ROUNDDOWN);
    }
    
    walker = nhMakeNeighborhoodWalker(nhood, size, num_dims,
                                      NH_SKIP_CENTER | NH_SKIP_LEADING);
    nhDestroyNeighborhood(nhood);
    
    num_sets = connected_components(BW, L, num_elements, walker);
    nhDestroyNeighborhoodWalker(walker);

    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        *(mxGetPr(plhs[1])) = (double) num_sets;
    }
}
