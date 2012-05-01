/*! \file blas_bsr_omp.c
 *  \brief BLAS operations for sparse matrices in BSR format.
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*-----------------------------------omp--------------------------------------*/
/* Feng Chunsheng Yue Xiaoqiang  Expand the inner loop  /mar/14/2011/  */
/*!
 * \fn void fasp_blas_dbsr_aAxpy_omp ( double alpha, dBSRmat *A, double *x, 
 *                                     double *y, int nthreads, int openmp_holds )
 * \brief Compute y := alpha*A*x + y
 * \param alpha a real number
 * \param A pointer to the matrix
 * \param x pointer to the vector x
 * \param y pointer to the vector y
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/14/2011
 */
void fasp_blas_dbsr_aAxpy_omp (double alpha, 
     dBSRmat *A, 
     double *x, 
     double *y, 
     int nthreads, 
     int openmp_holds )
{
#if FASP_USE_OPENMP
    /* members of A */
    int     ROW = A->ROW;
    int     nb  = A->nb;
    int    *IA  = A->IA;
    int    *JA  = A->JA;
    double *val = A->val;
    
    /* local variables */
    int     size = ROW*nb;
    int     jump = nb*nb;
    int     i,j,k, iend;
    double  temp = 0.0;
    double *pA   = NULL;
    double *px0  = NULL;
    double *py0  = NULL;
    double *py   = NULL;
    
    //----------------------------------------------
    //   Treat (alpha == 0.0) computation 
    //----------------------------------------------
    
    if (alpha == 0.0)
    {
    return; // Nothting to compute
    }
    
    //-------------------------------------------------
    //   y = (1.0/alpha)*y
    //-------------------------------------------------
    
    if (alpha != 1.0)
    {
    temp = 1.0 / alpha;
    fasp_blas_array_scale_omp(size, temp, y, nthreads, openmp_holds);
    }
    
    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    
    switch (nb)
    {
    case 2: 
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend) //num_threads(nthreads)
    for (myid =0; myid < nthreads; myid++)
    {
    //    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*2];
    iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*4; // &val[k*jump];
    px0 = x+j*2; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc2( pA, px0, py );
    }
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*2];
                    iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*4; // &val[k*jump];
    px0 = x+j*2; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc2( pA, px0, py );
    }
    }
    }
    }
    break;

    case 3: 
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend ) 
    for (myid =0; myid < nthreads; myid++)
    {
    //    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*3];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );
    }
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*3];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );
    }
    }
    }
    }
    break;
    
    case 5:
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend ) //num_threads(nthreads)
    for (myid =0; myid < nthreads; myid++)
    {
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*5];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );
    }  
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*5];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );
    }
    }
    }
    }
    break;
    
    case 7:
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend ) //num_threads(nthreads)
    for (myid =0; myid < nthreads; myid++)
    {
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*7];

                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );
    }  
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*7];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );
    }
    }
    }
    }
    break;
    
    default: 
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend) //num_threads(nthreads)
    for (myid =0; myid < nthreads; myid++)
    {
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*nb];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );
    }  
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*nb];
                        iend = IA[i+1];
    for (k = IA[i]; k < iend; ++k)
    {
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );
    }
    }
    }
    }
    break;
    }
    
    //------------------------------------------
    //   y = alpha*y
    //------------------------------------------
    
    if (alpha != 1.0)
    {
    fasp_blas_array_scale_omp(size, alpha, y, nthreads, openmp_holds);
    } 
#endif    
    return;
}

/*!
 * \fn void fasp_blas_dbsr_mxv_omp ( dBSRmat *A, double *x, double *y, int nthreads, int openmp_holds )
 * \brief Compute y := A*x
 * \param A pointer to the matrix
 * \param x pointer to the vector x
 * \param y pointer to the vector y
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/14/2011
 */
void fasp_blas_dbsr_mxv_omp (dBSRmat *A, 
     double *x, 
     double *y, 
     int nthreads, 
     int openmp_holds )
{
#if FASP_USE_OPENMP
    /* members of A */
    int     ROW = A->ROW;
    int     nb  = A->nb;
    int    *IA  = A->IA;
    int    *JA  = A->JA;
    double *val = A->val;
    
    /* local variables */
    int     size = ROW*nb;
    int     jump = nb*nb;
    int     i,j,k, num_nnz_row;
    double *pA  = NULL;
    double *px0 = NULL;
    double *py0 = NULL;
    double *py  = NULL;
    
    //-----------------------------------------------------------------
    //  zero out 'y' 
    //-----------------------------------------------------------------
    fasp_array_set_omp(size, y, 0.0, nthreads,openmp_holds);
    
    //-----------------------------------------------------------------
    //   y = A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    
    switch (nb)
    {
    case 3: 
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py) //num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*3];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*9;
    px0 = x+j*3;
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );
    }
    break;
    }
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*3];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*9; // &val[k*jump];
    px0 = x+j*3; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc3( pA, px0, py );
    }
    break;
    }
    }
    }
    }
    break;
    
    case 5:
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py) //num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*5];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );
    }
    break;
    }
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*5];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*25; // &val[k*jump];
    px0 = x+j*5; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc5( pA, px0, py );
    }
    break;
    }
    }
    }
    }
    break;
    
    case 7:
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py) //num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*7];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );
    }
    break;
    }
    }
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*7];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    k ++;
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*49; // &val[k*jump];
    px0 = x+j*7; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx_nc7( pA, px0, py );
    }
    break;
    }
    }
    }
    }
    break;
    
    default: 
    {
    if (ROW > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py) //num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
    for (i=mybegin; i < myend; ++i)
    {
    py0 = &y[i*nb];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );
    }
    break;
    }
    }  
    }
    }
    else {
    for (i = 0; i < ROW; ++i)
    {
    py0 = &y[i*nb];
    num_nnz_row = IA[i+1] - IA[i];
    switch(num_nnz_row)
    {
    case 3:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 4:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 5:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 6:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    case 7:
    k = IA[i];
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    k ++;
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );

    break;
    default:
    for (k = IA[i]; k < IA[i+1]; ++k)
    {
    j = JA[k];
    pA = val+k*jump; // &val[k*jump];
    px0 = x+j*nb; // &x[j*nb];
    py = py0;
    fasp_blas_smat_ypAx( pA, px0, py, nb );
    }
    break;
    }
    }
    }
    }
    break;
    }
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
