/*! \file blas_vec_omp.c
 *  \brief BLAS operations for vectors
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_dvec_axpy_omp (const double a, dvector *x, dvector *y, int nthreads, int openmp_holds)
 * \brief y = a*x + y
 * \param a real number
 * \param x pointer to dvector
 * \param y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_dvec_axpy_omp (const double a, 
     dvector *x, 
     dvector *y, 
     int nthreads, 
     int openmp_holds)
{
#if FASP_USE_OPENMP
    unsigned int i, m=x->row;
    double *xpt=x->val, *ypt=y->val;
    
    if ((y->row-m)!=0) {
    printf("Error: two vectors have different length!\n");
    exit(ERROR_DATA_STRUCTURE);
    }
    
    if (m > openmp_holds) {
    int myid, mybegin, myend;
#pragma omp parallel private(myid,mybegin,myend,i) ////num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
    for (i=mybegin; i<myend; ++i) ypt[i] += a*xpt[i];
    }
    }
    else {
    for (i=0; i<m; ++i) ypt[i] += a*xpt[i];
    }
#endif
}

/**
 * \fn void fasp_blas_dvec_axpyz_omp (const double a, dvector *x, dvector *y, dvector *z, int nthreads, int openmp_holds)
 * \brief z = a*x + y, z is a third vector (z is cleared)
 * \param a real number
 * \param x pointer to dvector
 * \param y pointer to dvector
 * \param z pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_dvec_axpyz_omp (const double a, 
    dvector *x, 
    dvector *y, 
    dvector *z, 
    int nthreads, 
    int openmp_holds)
{
#if FASP_USE_OPENMP
    const int m=x->row;
    double *xpt=x->val, *ypt=y->val, *zpt=z->val;
    
    if ((y->row-m)!=0) {
    printf("Error: two vectors have different length!\n");
    exit(ERROR_DATA_STRUCTURE);
    }
    
    z->row = m;
    fasp_array_cp_omp(m, ypt, zpt, nthreads,openmp_holds);
    fasp_blas_array_axpy_omp(m, a, xpt, zpt, nthreads,openmp_holds);
#endif
}

/**
 * \fn double fasp_blas_dvec_dotprod_omp (dvector *x, dvector *y, int nthreads, int openmp_holds)
 * \brief Inner product of two vectors (x,y)
 * \param x pointer to dvector
 * \param y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return Inner product
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_dvec_dotprod_omp (dvector *x, 
    dvector *y, 
    int nthreads, 
    int openmp_holds)
{
    double value=0;
#if FASP_USE_OPENMP
    int i;
    const int length=x->row;
    double *xpt=x->val, *ypt=y->val;    
    
    if (length > openmp_holds) {
#pragma omp parallel for reduction(+:value) private(i) ////num_threads(nthreads)
    for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
    }
    else {
    for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
    }
#endif
    return value;
}

/**
 * \fn double fasp_blas_dvec_relerr_omp (dvector *x, dvector *y, int nthreads, int openmp_holds)
 * \brief Relative error of two dvector x and y
 * \param x pointer to dvector
 * \param y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return relative error ||x-y||/||x||
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_dvec_relerr_omp (dvector *x, 
     dvector *y, 
     int nthreads, 
     int openmp_holds)
{
    double diff=0, temp=0;
#if FASP_USE_OPENMP
    int i;
    const int length=x->row;
    double *xpt=x->val, *ypt=y->val;
    
    if (length!=y->row) {
    printf("Error: The lengths of vectors do not match! \n");
    exit(ERROR_DUMMY_VAR);    
    }
    
    if (length > openmp_holds) {
#pragma omp parallel for reduction(+:temp,diff) private(i) ////num_threads(nthreads)
    for (i=0;i<length;++i) {
    temp += xpt[i]*xpt[i];
    diff += pow(xpt[i]-ypt[i],2);
    }
    }
    else {
    for (i=0;i<length;++i) {
    temp += xpt[i]*xpt[i];
    diff += pow(xpt[i]-ypt[i],2);
    }
    }
#endif    
    return sqrt(diff/temp);
}

/**
 * \fn double fasp_blas_dvec_norm1_omp (dvector *x, int nthreads, int openmp_holds)
 * \brief L1 norm of dvector x
 * \param x pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return L1 norm of x
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_dvec_norm1_omp (dvector *x, 
    int nthreads, 
    int openmp_holds)
{
    double onenorm=0;
#if FASP_USE_OPENMP
    int i;
    const int length=x->row;
    double *xpt=x->val;
    if (length > openmp_holds) {
#pragma omp parallel for reduction(+:onenorm) private(i) ////num_threads(nthreads)
    for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
    }
    else {
    for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
    }
#endif
    return onenorm;
}

/**
 * \fn double fasp_blas_dvec_norm2_omp (dvector *x, int nthreads, int openmp_holds)
 * \brief L2 norm of dvector x
 * \param x pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return L2 norm of x
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Jan/11/2011
 */
double fasp_blas_dvec_norm2_omp (dvector *x, 
    int nthreads, 
    int openmp_holds)
{
    double twonorm=0;
#if FASP_USE_OPENMP
    int i;
    const int length=x->row;
    double *xpt=x->val;
    if (length > openmp_holds) {
#pragma omp parallel for reduction(+:twonorm) private(i) 
    for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
    }
    else {
    for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
    }
#endif
    return sqrt(twonorm);
}
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
