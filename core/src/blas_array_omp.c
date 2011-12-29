/*! \file blas_array_omp.c
 *  \brief BLAS operations for arrays.
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
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_blas_array_scale_omp(int n, double a, double *x, 
 *                                    int nthreads, int openmp_holds) 
 * \brief x = a*x
 * \param n number of variables
 * \param a  a real number
 * \param *x pointer to the vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_array_scale_omp (int n, 
																const double a, 
																double *x, 
																int nthreads, 
																int openmp_holds)
{
#if FASP_USE_OPENMP
	int i;
	if (a == 1.0) {
	}
	else {
		if (n > openmp_holds) {
#pragma omp parallel for private(i) schedule(static)  
			for (i=0; i<n; ++i) x[i] *= a;
		}
		else {
			for (i=0; i<n; ++i) x[i] *= a;
		}
	}
#endif
}

/**
 * \fn void fasp_blas_array_axpy_omp (int n, const double a, double *x, 
 *                                    double *y, int nthreads, int openmp_holds)
 * \brief y = a*x + y
 * \param n number of variables
 * \param a a real number
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_array_axpy_omp (int n, 
															 const double a, 
															 double *x, 
															 double *y, 
															 int nthreads, 
															 int openmp_holds)
{
#if FASP_USE_OPENMP
	unsigned int i;
	
	if (a==1.0) {
		if (n > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
				FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
				for (i=mybegin; i<myend; ++i) y[i] += x[i];
			}
		}
		else {
			for (i=0; i<n; ++i) y[i] += x[i];
		}
	}
	else if (a==-1.0) {
		if (n > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
				FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
				for (i=mybegin; i<myend; ++i) y[i] -= x[i];
			}
		}
		else {
			for (i=0; i<n; ++i) y[i] -= x[i];
		}
	}
	else {
		if (n > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
				FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
				for (i=mybegin; i<myend; ++i) y[i] += a*x[i];
			}
		}
		else {
			for (i=0; i<n; ++i) y[i] += a*x[i];
		}
	}
#endif
}

/**
 * \fn void fasp_blas_array_axpyz_omp (int n, const double a, double *x, 
 *                         double *y, double *z, int nthreads, int openmp_holds) 
 * \brief z = a*x + y, z is the third vector
 * \param n number of variables
 * \param a a real number
 * \param *x pointer to the original vector 1
 * \param *y pointer to the original vector 2
 * \param *z pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_array_axpyz_omp (int n, 
																const double a, 
																double *x, 
																double *y, 
																double *z, 
																int nthreads, 
																int openmp_holds)
{
#if FASP_USE_OPENMP
	unsigned int i;
	if (n > openmp_holds) {
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
			for (i=mybegin; i<myend; ++i) z[i] = a*x[i]+y[i];
		}
	}
	else {
		for (i=0; i<n; ++i) z[i] = a*x[i]+y[i];
	}
#endif
}

/**
 * \fn void fasp_blas_array_axpby_omp (int n, const double a, double *x, 
 *                    const double b, double *y, int nthreads, int openmp_holds) 
 * \brief y = a*x + b*y
 * \param n number of variables
 * \param a real number
 * \param b real number
 * \param *x pointer to the origianl vector
 * \param *y pointer to the destination vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_array_axpby_omp (int n, 
																const double a, 
																double *x, 
																const double b, 
																double *y, 
																int nthreads, 
																int openmp_holds) 
{
#if FASP_USE_OPENMP
	unsigned int i;
	if (n > openmp_holds) {
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
			for (i=mybegin; i<myend; ++i) y[i] = a*x[i]+b*y[i];
		}
	}
	else {
		for (i=0; i<n; ++i) y[i] = a*x[i]+b*y[i];
	}
#endif
}

/**
 * \fn double fasp_blas_array_dotprod_omp (int n, double *x, double *y, 
 *                                         int nthreads, int openmp_holds) 
 * \brief Inner product of two arraies (x,y)
 * \param n number of variables
 * \param *x pointer to vector 1
 * \param *y pointer to vector 2
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return inner product
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_array_dotprod_omp (int n, 
																		double *x, 
																		double *y, 
																		int nthreads, 
																		int openmp_holds)
{
	double value=0.0;
#if FASP_USE_OPENMP
	int i;
	
	if (n > openmp_holds) {
#pragma omp parallel for reduction(+:value) private(i) ////num_threads(nthreads)
		for (i=0;i<n;++i) value+=x[i]*y[i];
	}
	else {
		for (i=0;i<n;++i) value+=x[i]*y[i];
	}
#endif
	return value;
}

/**
 * \fn double fasp_blas_array_norm1_omp (int n, double *x, int nthreads, 
 *                                       int openmp_holds) 
 * \brief L1 norm of array x
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return L1 norm of x
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_array_norm1_omp (int n, 
																	double *x, 
																	int nthreads, 
																	int openmp_holds)
{
	double onenorm=0;
#if FASP_USE_OPENMP
	int i;
	if (n > openmp_holds) {
#pragma omp parallel for reduction(+:onenorm) private(i) ////num_threads(nthreads)
		for (i=0;i<n;++i) onenorm+=ABS(x[i]);
	}
	else {
		for (i=0;i<n;++i) onenorm+=ABS(x[i]);
	}
#endif
	return onenorm;
}

/**
 * \fn double fasp_blas_array_norm2_omp (int n, double *x, int nthreads, 
 *                                       int openmp_holds)
 * \brief L2 norm of array x
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return L2 norm of x
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_array_norm2_omp (int n, 
																	double *x, 
																	int nthreads, 
																	int openmp_holds)
{
	double twonorm=0;
#if FASP_USE_OPENMP
	int i;
	if (n > openmp_holds) {
#pragma omp parallel for reduction(+:twonorm) private(i) ////num_threads(nthreads)
		for (i=0;i<n;++i) twonorm+=x[i]*x[i];
	}
	else {
		for (i=0;i<n;++i) twonorm+=x[i]*x[i];
	}
#endif
	return sqrt(twonorm);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
