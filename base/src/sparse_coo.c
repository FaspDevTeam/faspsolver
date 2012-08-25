/*! \file sparse_coo.c
 *  \brief Functions for COO sparse matrices. 
 */

#include <math.h>
#include <omp.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dCOOmat fasp_dcoo_create (INT m, INT n, INT nnz)
 *
 * \brief Create IJ sparse matrix data memory space
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   The new dCOOmat matrix
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
 */
dCOOmat fasp_dcoo_create (INT m, 
                          INT n, 
                          INT nnz)
{    
    dCOOmat A;
    
    A.I   = (INT *)fasp_mem_calloc(nnz, sizeof(INT));     
    A.J   = (INT *)fasp_mem_calloc(nnz, sizeof(INT));     
    A.val = (REAL *)fasp_mem_calloc(nnz, sizeof(REAL)); 
    
    A.row=m; A.col=n; A.nnz=nnz;
    
    return A;
}

/**
 * \fn void fasp_dcoo_free (dCOOmat *A)
 *
 * \brief Free IJ sparse matrix data memeory space
 *
 * \param A   Pointer to the dCOOmat matrix
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
 */
void fasp_dcoo_free (dCOOmat *A)
{    
    if (A==NULL) return;
    
    fasp_mem_free(A->I);   A->I   = NULL;
    fasp_mem_free(A->J);   A->J   = NULL;
    fasp_mem_free(A->val); A->val = NULL;
}

/**
 * \fn void fasp_dcoo_shift (dCOOmat *A, INT offset)
 *
 * \brief Reindex a REAL matrix in IJ format to make the index starting from 0 or 1.
 *
 * \param A       Pointer to IJ matrix
 * \param offset  Size of offset (1 or -1)
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
 *
 * Modified by Chunsheng Feng, Zheng Li
 * \date  08/25/2012
 */
void fasp_dcoo_shift (dCOOmat *A,
                      INT offset)
{
    const INT nnz=A->nnz;
    INT i, *ai=A->I, *aj=A->J;
    
    // Variables for Openmp
    INT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;

#ifdef _OPENMP
    if (nnz > OPENMP_HOLDS) {
        use_openmp = TRUE;
	nthreads = FASP_GET_NUM_THREADS();
    }
#endif

    if (offset == 0) offset = ISTART;
    
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend)
#endif
        for (myid=0; myid<nthreads; myid++) {
	    FASP_GET_START_END(myid, nthreads, nnz, &mybegin, &myend);	
            for (i=mybegin; i<myend; ++i) {    
                ai[i]+=offset; aj[i]+=offset;
            }
        }
    }
    else {
        for (i=0;i<nnz;++i) {    
            ai[i]+=offset; aj[i]+=offset;
        }
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
