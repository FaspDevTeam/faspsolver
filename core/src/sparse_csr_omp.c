/*! \file sparse_csr_omp.c
 *  \brief Functions for CSR sparse matrices. 
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void void fasp_dcsr_getdiag_omp (INT n, dCSRmat *A, dvector *diag, INT nthreads, INT openmp_holds) 
 *
 * \brief Get first n diagonal entries of a CSR matrix A
 *
 * \param n     an interger (if n=0, then get all diagonal entries)
 * \param A     pointer to dCSRmat CSR matrix
 * \param diag  pointer to the diagonal as a dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012 Modified by FENG Chunsheng
 */
void fasp_dcsr_getdiag_omp (INT n, 
                            dCSRmat *A, 
                            dvector *diag, 
                            INT nthreads, 
                            INT openmp_holds) 
{
#if FASP_USE_OPENMP
    INT i,k,j,ibegin,iend;    
    
    if (n==0) n=MIN(A->row,A->col);
    
    fasp_dvec_alloc(n,diag);
    
    if (n > openmp_holds) {
        INT mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, k, j) 
        for (myid = 0; myid < nthreads; myid++ )
        {
    FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
    for (i=mybegin; i<myend; i++)
    {
    ibegin=A->IA[i]; iend=A->IA[i+1];
    for (k=ibegin;k<iend;++k) {
    j=N2C(A->JA[N2C(k)]);
    if ((j-i)==0) {
    diag->val[i] = A->val[N2C(k)]; break;
    } // end if
    } // end for k
    } // end for i
    }
    }
    else {
    for (i=0;i<n;++i) {
    ibegin=A->IA[i]; iend=A->IA[i+1];
    for (k=ibegin;k<iend;++k) {
    j=N2C(A->JA[N2C(k)]);
    if ((j-i)==0) {
    diag->val[i] = A->val[N2C(k)]; break;
    } // end if
    } // end for k
    } // end for i
    }
#endif
}

/**
 * \fn void fasp_dcsr_cp_omp (dCSRmat *A, dCSRmat *B, int nthreads, int openmp_holds)
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   pointer to the dCSRmat matrix
 * \param B   pointer to the dCSRmat matrix
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dcsr_cp_omp (dCSRmat *A, 
                       dCSRmat *B, 
                       int nthreads, 
                       int openmp_holds)
{
#if FASP_USE_OPENMP
    B->row=A->row;
    B->col=A->col;
    B->nnz=A->nnz;
    
    fasp_iarray_cp_omp(A->row+1, A->IA, B->IA, nthreads, openmp_holds);
    fasp_iarray_cp_omp(A->nnz, A->JA, B->JA, nthreads, openmp_holds);
    fasp_array_cp_omp(A->nnz, A->val, B->val, nthreads, openmp_holds);
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
