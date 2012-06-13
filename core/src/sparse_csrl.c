/*! \file sparse_csrl.c
 *  \brief Functions for sparse matrices in CSRL format
 *
 *  \note For details of CSRL format, refer to 
 *        Optimizing sparse matrix vector product computations using unroll and jam
 *        by John Mellor-Crummey and John Garvin, Tech Report Rice Univ, Aug 2002.
 *
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dCSRLmat * fasp_dcsrl_create (INT num_rows, INT num_cols, INT num_nonzeros) 
 *
 * \brief Create a dCSRLmat object
 *
 * \param num_rows      Number of rows
 * \param num_cols      Number of cols
 * \param num_nonzeros  Number of nonzero entries
 *
 * \author Zhou Zhiyang
 * \date   2011/01/07
 */
dCSRLmat * fasp_dcsrl_create (INT num_rows, 
                              INT num_cols, 
                              INT num_nonzeros)
{
    dCSRLmat *A   = (dCSRLmat *)fasp_mem_calloc(1, sizeof(dCSRLmat));
    
    A -> row      = num_rows;
    A -> col      = num_cols;
    A -> nnz      = num_nonzeros;
    A -> nz_diff  = NULL;
    A -> index    = NULL;
    A -> start    = NULL;
    A -> ja       = NULL;
    A -> val      = NULL;
    
    return A;
}

/**
 * \fn void fasp_dcsrl_free ( dCSRLmat *A )
 *
 * \brief Destroy a dCSRLmat object
 *
 * \param A   Pointer to the dCSRLmat type matrix
 *
 * \author Zhou Zhiyang
 * \date   2011/01/07
 */
void fasp_dcsrl_free (dCSRLmat *A)
{
    if (A) {  
        if (A -> nz_diff) free(A -> nz_diff);
        if (A -> index)   free(A -> index);
        if (A -> start)   free(A -> start);
        if (A -> ja)      free(A -> ja);
        if (A -> val)     free(A -> val);
        free(A);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
