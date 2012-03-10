/*! \file precond_bsr_omp.c
 *  \brief Preconditioners for sparse matrices in BSR format.
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_dbsr_diag_omp (double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief Diagonal preconditioner z=inv(D)*r
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \note: works for general nb
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Nov/28/2011
 */
void fasp_precond_dbsr_diag_omp (double *r, 
                                 double *z, 
                                 void *data, 
                                 int nthreads, 
                                 int openmp_holds)
{
#if FASP_USE_OPENMP
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	const int nb = diag->nb; 
	
	switch (nb)
	{
		case 2:
			fasp_precond_dbsr_diag_nc2_omp( r, z, diag, nthreads, openmp_holds );
			break;
		case 3:
			fasp_precond_dbsr_diag_nc3_omp( r, z, diag, nthreads, openmp_holds );
			break;
			
		case 5:
			fasp_precond_dbsr_diag_nc5_omp( r, z, diag, nthreads, openmp_holds );
			break;
			
		case 7:
			fasp_precond_dbsr_diag_nc7_omp( r, z, diag, nthreads, openmp_holds );
			break;
			
		default:
		{
			double *diagptr = diag->diag.val;
			const int nb2 = nb*nb;
			const int m = diag->diag.row/nb2;	
			
			unsigned int i;
			if (m > openmp_holds) {
				int myid;
				int mybegin;
				int myend;
				int stride_i = m/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i) num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1)  myend = mybegin+stride_i;
					else myend = m;
					for (i=mybegin; i < myend; ++i)
					{
						fasp_blas_smat_mxv(&(diagptr[i*nb2]),&(r[i*nb]),&(z[i*nb]),nb);
					}
				}
			}
			else {
				for (i = 0; i < m; ++i) 
				{
					fasp_blas_smat_mxv(&(diagptr[i*nb2]),&(r[i*nb]),&(z[i*nb]),nb);
				}
			}
		}
			break;
	}
#endif
}

/**
 * \fn void fasp_precond_dbsr_diag_nc2_omp (double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \note: works for 2-component 
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Nov/28/2011
 */
void fasp_precond_dbsr_diag_nc2_omp (double *r, 
                                     double *z, 
                                     void *data, 
                                     int nthreads, 
                                     int openmp_holds )
{
#if FASP_USE_OPENMP
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/4;	
	
	unsigned int i;
	if (m > openmp_holds) {
		int myid;
		int mybegin;
		int myend;
		int stride_i = m/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i) ///num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1)  myend = mybegin+stride_i;
			else myend = m;
			for (i=mybegin; i < myend; ++i)
			{
				fasp_blas_smat_mxv_nc2(&(diagptr[i*4]),&(r[i*2]),&(z[i*2]));
			}
		}
	}
	else {
		for (i = 0; i < m; ++i) 
		{
			fasp_blas_smat_mxv_nc2(&(diagptr[i*4]),&(r[i*2]),&(z[i*2]));
		}
	}
#endif
}

/**
 * \fn void fasp_precond_dbsr_diag_nc3_omp (double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \note: works for 3-component 
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Nov/28/2011
 */
void fasp_precond_dbsr_diag_nc3_omp (double *r, 
                                     double *z, 
                                     void *data, 
                                     int nthreads, 
                                     int openmp_holds )
{
#if FASP_USE_OPENMP
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/9;	
	
	unsigned int i;
	if (m > openmp_holds) {
		int myid;
		int mybegin;
		int myend;
		int stride_i = m/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1)  myend = mybegin+stride_i;
			else myend = m;
			for (i=mybegin; i < myend; ++i)
			{
				fasp_blas_smat_mxv_nc3(&(diagptr[i*9]),&(r[i*3]),&(z[i*3]));
			}
		}
	}
	else {
		for (i = 0; i < m; ++i) 
		{
			fasp_blas_smat_mxv_nc3(&(diagptr[i*9]),&(r[i*3]),&(z[i*3]));
		}
	}
#endif
}

/**
 * \fn void fasp_precond_dbsr_diag_nc5_omp (double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note: works for 5-component
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Nov/28/2011
 */
void fasp_precond_dbsr_diag_nc5_omp (double *r, 
                                     double *z, 
                                     void *data, 
                                     int nthreads, 
                                     int openmp_holds)
{
#if FASP_USE_OPENMP
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/25;
	
	unsigned int i;
	if (m > openmp_holds) {
		int myid;
		int mybegin;
		int myend;
		int stride_i = m/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1)  myend = mybegin+stride_i;
			else myend = m;
			for (i=mybegin; i < myend; ++i)
			{
				fasp_blas_smat_mxv_nc5(&(diagptr[i*25]),&(r[i*5]),&(z[i*5]));
			}
		}
	}
	else {
		for (i = 0; i < m; ++i) 
		{
			fasp_blas_smat_mxv_nc5(&(diagptr[i*25]),&(r[i*5]),&(z[i*5]));
		}
	}
#endif
}

/**
 * \fn void fasp_precond_dbsr_diag_nc7_omp ( double *r, double *z, void *data, int nthreads, int openmp_holds )
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note: works for 7-component
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Nov/28/2011
 */
void fasp_precond_dbsr_diag_nc7_omp (double *r, 
                                     double *z, 
                                     void *data, 
                                     int nthreads, 
                                     int openmp_holds )
{
#if FASP_USE_OPENMP
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/49;	
	
	unsigned int i;
	if (m > openmp_holds) {
		int myid;
		int mybegin;
		int myend;
		int stride_i = m/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1)  myend = mybegin+stride_i;
			else myend = m;
			for (i=mybegin; i < myend; ++i)
			{
				fasp_blas_smat_mxv_nc7(&(diagptr[i*49]),&(r[i*7]),&(z[i*7]));
			}
		}
	}
	else {
		for (i = 0; i < m; ++i) 
		{
			fasp_blas_smat_mxv_nc7(&(diagptr[i*49]),&(r[i*7]),&(z[i*7]));
		}
	}
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
