/*! \file array_omp.c
 *  \brief Array operations
 * 
 *  Some simple array operations ...
 *
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_array_set_omp (int n, double *x, double val, int nthreads, int openmp_holds)
 * \brief Set initial value for an array to be x=val
 * \param n number of variables
 * \param x pointer to the vector
 * \param val initial value for the double array
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_array_set_omp (int n, 
                         double *x, 
                         double val, 
                         int nthreads, 
                         int openmp_holds)
{
#if FASP_USE_OPENMP
    unsigned int i;
    if (val == 0.0) {
        if (n > openmp_holds) {
            int mybegin,myend,myid;
#pragma omp parallel for private(myid,mybegin,myend)
            for (myid = 0; myid < nthreads; myid ++)
                {
                    FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                    memset(&x[mybegin], 0x0, sizeof(double)*(myend-mybegin));
                }
        }
        else 
            memset(x, 0x0, sizeof(double)*n);
    }
    else {
        
        if (n > openmp_holds) {
            int mybegin,myend,myid;
#pragma omp parallel for private(myid,mybegin,myend,i)
            for (myid = 0; myid < nthreads; myid ++)
                {
                    FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                    for (i=mybegin; i<myend; ++i) x[i]=val;
                }
        }
        else {
            for (i=0; i<n; ++i) x[i]=val;
        }
    }
#endif
}

/**
 * \fn void fasp_iarray_set_omp (int n, int *x, int val, int nthreads, int openmp_holds)
 * \brief Set initial value for an array to be x=val
 * \param n number of variables
 * \param x pointer to the vector
 * \param val initial value for the double array
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_iarray_set_omp (int n, int *x, int val, int nthreads, int openmp_holds)
{
#if FASP_USE_OPENMP
    unsigned int i;
    if (val == 0) {
        if (n > openmp_holds) {
            int mybegin,myend,myid;
#pragma omp parallel for private(myid,mybegin,myend)
            for (myid = 0; myid < nthreads; myid ++)
                {
                    FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                    memset(&x[mybegin], 0, sizeof(int)*(myend-mybegin));
                }
        }
        else {
            memset(x, 0, sizeof(int)*n);
        }
    }
    else {
        if (n > openmp_holds) {
            int mybegin,myend,myid;
#pragma omp parallel for private(myid,mybegin,myend,i)
            for (myid = 0; myid < nthreads; myid ++)
                {
                    FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                    for (i=mybegin; i<myend; ++i) x[i]=val;
                }
        }
        else {
            for (i=0; i<n; ++i) x[i]=val;
        }
    }
#endif
}

/**
 * \fn void fasp_array_cp_omp (int n, double *x, double *y, int nthreads, int openmp_holds) 
 * \brief Copy an array to the other y=x
 * \param n number of variables
 * \param x pointer to the original vector
 * \param y pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012  Modified by FENG Chunsheng
 */
void fasp_array_cp_omp (int n, 
                        double *x, 
                        double *y, 
                        int nthreads, 
                        int openmp_holds)
{
#if FASP_USE_OPENMP
#if 1 
    if (n > openmp_holds) {
        int mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin,myend) 
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                memcpy(&y[mybegin],&x[mybegin], sizeof(double)*(myend-mybegin));
            }
    }
    else
#endif    
        {
            memcpy(y,x,n*sizeof(double));
        }
#endif
}

/**
 * \fn void fasp_iarray_cp_omp (int n, int *x, int *y, int nthreads, int openmp_holds) 
 * \brief Copy an array to the other y=x
 * \param n number of variables
 * \param x pointer to the original vector
 * \param y pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012  Modified by FENG Chunsheng
 */
void fasp_iarray_cp_omp (int n, 
                         int *x, 
                         int *y, 
                         int nthreads, 
                         int openmp_holds) 
{
#if FASP_USE_OPENMP
#if 1 
    if (n > openmp_holds) {
        int mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin,myend) 
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                memcpy(&y[mybegin],&x[mybegin], sizeof(int)*(myend-mybegin));
            }
    }
    else
#endif    
        {
            memcpy(y,x,n*sizeof(int));
        }
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

