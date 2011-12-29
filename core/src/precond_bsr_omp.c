/*! \file precond_bsr_omp.c
 *  \brief Preconditioners for sparse matrices in BSR format.
 */

#include "fasp.h"
#include "fasp_functs.h"
/*---------------------------------omp----------------------------------------*/

/**
 * \fn void fasp_precond_dbsr_diag_omp (double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief Diagonal preconditioner z=inv(D)*r
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \note: works for general nb (Xiaozhe)
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
				int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend,i) ////num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
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
 * \fn void fasp_precond_dbsr_diag_nc3_omp (double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \note: works for 3-component (Xiaozhe)
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
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
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
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note: works for 5-component (Xiaozhe)
 * \date 01/06/2011
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
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
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
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note: works for 7-component (Xiaozhe)
 * \date 01/06/2011
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
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
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
