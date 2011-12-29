/*! \file array_omp.c
 *  \brief Array operations
 * 
 *  Some simple array operations ...
 *
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_array_set_omp (INT n, REAL *x, REAL val, INT nthreads, INT openmp_holds)
 * \brief Set initial value for an array to be x=val
 * \param n number of variables
 * \param *x pointer to the vector
 * \param val initial value for the REAL array
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_array_set_omp (INT n, 
                         REAL *x, 
                         REAL val, 
                         INT nthreads, 
                         INT openmp_holds)
{
#if FASP_USE_OPENMP
	INT i;
	if (val == 0.0) {
		if (n > openmp_holds) {
			INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
				FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
				memset(&x[mybegin], 0x0, sizeof(REAL)*(myend-mybegin));
			}
		}
		else {
			memset(x, 0x0, sizeof(REAL)*n);
		}
	}
	else {
		if (n > openmp_holds) {
#pragma omp parallel for private(i) schedule(static)
			for (i=0; i<n; ++i) x[i]=val;
		}
		else {
			for (i=0; i<n; ++i) x[i]=val;
		}
	}
#endif
}

/**
 * \fn void fasp_iarray_set_omp (INT n, REAL *x, REAL val, INT nthreads, INT openmp_holds)
 * \brief Set initial value for an array to be x=val
 * \param n number of variables
 * \param *x pointer to the vector
 * \param val initial value for the REAL array
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_iarray_set_omp (INT n, INT *x, INT val, INT nthreads, INT openmp_holds)
{
#if FASP_USE_OPENMP
	INT i;
	if (val == 0) {
		if (n > openmp_holds) {
			INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
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
#pragma omp parallel for private(i) schedule(static) ////num_threads(nthreads)
			for (i=0; i<n; ++i) x[i]=val;
		}
		else {
			for (i=0; i<n; ++i) x[i]=val;
		}
	}
#endif
}

/**
 * \fn void fasp_array_cp_omp (INT n, REAL *x, REAL *y, INT nthreads, INT openmp_holds) 
 * \brief Copy an array to the other y=x
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_array_cp_omp (INT n, 
                        REAL *x, 
                        REAL *y, 
                        INT nthreads, 
                        INT openmp_holds)
{
#if FASP_USE_OPENMP
	if (n > openmp_holds) {
		INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
			memcpy(&y[mybegin], &x[mybegin], sizeof(REAL)*(myend-mybegin));
		}
	}
	else {
		memcpy(y, x, n*sizeof(REAL));
	}
#endif
}

/**
 * \fn void fasp_iarray_cp_omp (INT n, INT *x, INT *y, INT nthreads, INT openmp_holds) 
 * \brief Copy an array to the other y=x
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_iarray_cp_omp (INT n, 
                         INT *x, 
                         INT *y, 
                         INT nthreads, 
                         INT openmp_holds) 
{
#if FASP_USE_OPENMP
	if (n > openmp_holds) {
		INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
			memcpy(&y[mybegin], &x[mybegin], sizeof(int)*(myend-mybegin));
		}
	}
	else {
		memcpy(y, x, n*sizeof(int));
	}
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

