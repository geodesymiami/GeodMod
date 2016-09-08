/*
 * Copyright 1993-2007 The MathWorks, Inc. 
 * $Revision: 1.2.4.2 $  $Date: 2007/06/04 21:10:12 $
 *
 * This file contains a function body for 2-D order-statistic
 * filtering.
 */

(TYPE *pA, TYPE *pB,
 mwSize startRow, mwSize startCol, 
 mwSize Mb, mwSize Nb, mwSize Ma, mwSize order, 
 mwSignedIndex *offsets, mwSize numOffsets
#ifdef ADD_OFFSET
 , double *add
#endif /* ADD_OFFSET */
)
{
    TYPE *vector;
    TYPE *p;
    mwSize col;
    mwSize row;
    mwSize k;

    vector = (TYPE *) mxCalloc(numOffsets, sizeof(*vector));
    
    for (col = 0; col < Nb; col++)
    {
        p = pA + (startCol+col)*Ma + startRow;
        for (row = 0; row < Mb; row++)
        {
            for (k = 0; k < numOffsets; k++)
            {
#ifdef ADD_OFFSET
		vector[k] = *(p + offsets[k]) + add[k];
#else
                vector[k] = *(p + offsets[k]);
#endif

            }
            *pB++ = 
		SELECT
		(vector, numOffsets, order);
            p++;
        }
    }

    mxFree((void *) vector);
}

#undef TYPE
#undef SELECT
#undef ADD_OFFSET

