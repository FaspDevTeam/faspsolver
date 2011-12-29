/*! \file vec.c
 *  \brief Simple operations for vectors (INT and REAL). 
 *
 *  \note 
 *  Every structures should be initialized before usage.
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_dvec_set_omp (int n, dvector *x, double val, int nthreads, int openmp_holds)
 * \brief Initialize dvector x=val
 *
 * \param n number of variables
 * \param *x pointer to dvector
 * \param val initial value for the dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dvec_set_omp (int n, 
                        dvector *x, 
                        double val, 
                        int nthreads, 
                        int openmp_holds)
{
#if FASP_USE_OPENMP
	int i;
	double *xpt=x->val;
	
	if (n>0) x->row=n;
	else n=x->row;
	
	if (val == 0.0) {
		if (n > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
				FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
				memset(&xpt[mybegin], 0x0, sizeof(double)*(myend-mybegin));
			}
		}
		else {
			memset(xpt, 0x0, sizeof(double)*n);
		}
	}
	else {
		if (n > openmp_holds) {
#pragma omp parallel for private(i) schedule(static)
			for (i=0; i<n; ++i) xpt[i]=val;
		}
		else
		{
			for (i=0; i<n; ++i) xpt[i]=val;
		}
	}
#endif
}

/**
 * \fn void fasp_dvec_cp_omp (dvector *x, dvector *y, int nthreads, int openmp_holds) 
 * \brief Copy dvector x to dvector y
 *
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dvec_cp_omp (dvector *x, 
                       dvector *y, 
                       int nthreads, 
                       int openmp_holds) 
{
#if FASP_USE_OPENMP
	int row=x->row;
	double *x_data=x->val,*y_data=y->val;
	
	y->row=row;
	if (row > openmp_holds) {
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			memcpy(&y_data[mybegin], &x_data[mybegin], sizeof(double)*(myend-mybegin));
		}
	}
	else {
		memcpy(y_data,x_data,row*sizeof(double));
	}
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
