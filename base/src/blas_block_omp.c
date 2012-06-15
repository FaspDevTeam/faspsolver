/*! \file blas_block_omp.c
 *  \brief BLAS operations for sparse matrices in block CSR format.
 *
 */

#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_bdbsr_aAxpy_omp (double alpha, block_BSR *A, double *x, double *y,
 *                                     int nthreads, int openmp_holds)
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 * \param alpha real number
 * \param A pointer to block_BSR matrix
 * \param x pointer to dvector
 * \param y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_bdbsr_aAxpy_omp (double alpha, 
                                block_BSR *A, 
                                double *x, 
                                double *y, 
                                int nthreads, 
                                int openmp_holds)
{
#if FASP_USE_OPENMP
    register dBSRmat *Arr = &(A->ResRes);
    register dCSRmat *Arw = &(A->ResWel);
    register dCSRmat *Awr = &(A->WelRes);
    register dCSRmat *Aww = &(A->WelWel);
    
    unsigned int Nr = Arw->row;
    
    register double *xr = x;
    register double *xw = &(x[Nr]);
    register double *yr = y;
    register double *yw = &(y[Nr]);
    
    // yr = alpha*Arr*xr + alpha*Arw*xw + yr
    fasp_blas_dbsr_aAxpy_omp(alpha, Arr, xr, yr, nthreads, openmp_holds);
    fasp_blas_dcsr_aAxpy_omp(alpha, Arw, xw, yr, nthreads, openmp_holds); 
    
    // yw = alpha*Awr*xr + alpha*Aww*xw + yw
    fasp_blas_dcsr_aAxpy_omp(alpha, Awr, xr, yw, nthreads, openmp_holds); 
    fasp_blas_dcsr_aAxpy_omp(alpha, Aww, xw, yw, nthreads, openmp_holds); 
#endif    
}

/**
 * \fn void fasp_blas_bdbsr_mxv_omp (block_BSR *A, double *x, double *y, int nthreads, int openmp_holds)
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A pointer to block_BSR matrix
 * \param x pointer to dvector
 * \param y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_bdbsr_mxv_omp (block_BSR *A, 
                              double *x, 
                              double *y, 
                              int nthreads, 
                              int openmp_holds)
{
#if FASP_USE_OPENMP
    register dBSRmat *Arr = &(A->ResRes);
    register dCSRmat *Arw = &(A->ResWel);
    register dCSRmat *Awr = &(A->WelRes);
    register dCSRmat *Aww = &(A->WelWel);
    
    unsigned int Nr = Arw->row;
    
    register double *xr = x;
    register double *xw = &(x[Nr]);
    register double *yr = y;
    register double *yw = &(y[Nr]);
    
    // yr = alpha*Arr*xr + alpha*Arw*xw + yr
    fasp_blas_dbsr_mxv_omp(Arr, xr, yr, nthreads, openmp_holds);
    fasp_blas_dcsr_aAxpy_omp(1.0, Arw, xw, yr, nthreads, openmp_holds);
    
    // yw = alpha*Awr*xr + alpha*Aww*xw + yw
    fasp_blas_dcsr_mxv_omp(Awr, xr, yw, nthreads, openmp_holds);
    fasp_blas_dcsr_aAxpy_omp(1.0, Aww, xw, yw, nthreads, openmp_holds); 
#endif    
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
