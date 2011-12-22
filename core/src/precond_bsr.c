/*! \file precond_bsr.c
 *  \brief Preconditioners for sparse matrices in BSR format.
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_dbsr_diag (double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Zhou Zhiyang, Xiaozhe Hu
 *
 * \note: works for general nb (Xiaozhe)
 * \date 10/26/2010
 */
void fasp_precond_dbsr_diag (double *r, 
														 double *z, 
														 void *data)
{
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	const int nb = diag->nb; 
	
	switch (nb)
	{
		case 2:
			fasp_precond_dbsr_diag_nc2( r, z, diag );
			break;
			
		case 3:
			fasp_precond_dbsr_diag_nc3( r, z, diag );
			break;
			
		case 5:
			fasp_precond_dbsr_diag_nc5( r, z, diag );
			break;
			
		case 7:
			fasp_precond_dbsr_diag_nc7( r, z, diag );
			break;
			
		default:
		{
			double *diagptr = diag->diag.val;
			const int nb2 = nb*nb;
			const int m = diag->diag.row/nb2;	
			
			unsigned int i;
			for (i = 0; i < m; ++i) 
			{
				fasp_blas_smat_mxv(&(diagptr[i*nb2]),&(r[i*nb]),&(z[i*nb]),nb);
			}
		}
			break;
	}
	
}

/**
 * \fn void fasp_precond_dbsr_diag_nc2 (double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Zhou Zhiyang, Xiaozhe Hu
 *
 * \note: works for 2-component (Xiaozhe)
 * \date 11/18/2011
 */
void fasp_precond_dbsr_diag_nc2 (double *r, 
								 double *z, 
								 void *data )
{
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/4;	
	
	unsigned int i;
	for (i = 0; i < m; ++i) 
	{
		fasp_blas_smat_mxv_nc2(&(diagptr[i*4]),&(r[i*2]),&(z[i*2]));
	}	
}

/**
 * \fn void fasp_precond_dbsr_diag_nc3 (double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Zhou Zhiyang, Xiaozhe Hu
 *
 * \note: works for 3-component (Xiaozhe)
 * \date 01/06/2011
 */
void fasp_precond_dbsr_diag_nc3 (double *r, 
																 double *z, 
																 void *data )
{
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/9;	
	
	unsigned int i;
	for (i = 0; i < m; ++i) 
	{
		fasp_blas_smat_mxv_nc3(&(diagptr[i*9]),&(r[i*3]),&(z[i*3]));
	}	
}

/**
 * \fn void fasp_precond_dbsr_diag_nc5 (double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Zhou Zhiyang, Xiaozhe Hu
 *
 * \note: works for 5-component (Xiaozhe)
 * \date 01/06/2011
 */
void fasp_precond_dbsr_diag_nc5 (double *r, 
																 double *z, 
																 void *data )
{
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/25;	
	
	unsigned int i;
	for (i = 0; i < m; ++i) 
	{
		fasp_blas_smat_mxv_nc5(&(diagptr[i*25]),&(r[i*5]),&(z[i*5]));
	}	
}

/**
 * \fn void fasp_precond_dbsr_diag_nc7 (double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r.
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Zhou Zhiyang, Xiaozhe Hu
 *
 * \note: works for 7-component (Xiaozhe)
 * \date 01/06/2011
 */
void fasp_precond_dbsr_diag_nc7 (double *r, 
																 double *z, 
																 void *data )
{
	precond_diagbsr *diag   = (precond_diagbsr *)data;
	double          *diagptr = diag->diag.val;
	
	const int m = diag->diag.row/49;	
	
	unsigned int i;
	for (i = 0; i < m; ++i) 
	{
		fasp_blas_smat_mxv_nc7(&(diagptr[i*49]),&(r[i*7]),&(z[i*7]));
	}	
}

/**
 * \fn void fasp_precond_dbsr_ilu (double *r, double *z, void *data)
 * \brief preconditioning using ILU decomposition
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Shiquan Zhang
 *
 * \note: works for general nb (Xiaozhe)
 * \date 11/09/2010
 */
void fasp_precond_dbsr_ilu (double *r, 
														double *z, 
														void *data)
{
	ILU_data    *iludata=(ILU_data *)data;
	const int    m=iludata->row, mm1=m-1, mm2=m-2, memneed=2*m;
	const int    nb=iludata->nb, nb2=nb*nb, size=m*nb;
	int         *ijlu=iludata->ijlu;
	double       *lu=iludata->luval;
	
	int         ib, ibstart,ibstart1;
	int         i, j, jj, begin_row, end_row;
	double      *zz, *zr, *mult;	   
	
	if (iludata->nwork<memneed) {
		printf("Error: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
		exit(ERROR_ALLOC_MEM);
	}
	
	zz   = iludata->work; 
	zr   = zz + size;
	mult = zr + size;
	
	memcpy(zr, r, size*sizeof(double));
	
	switch (nb) {
			
		case 1:
			
			// forward sweep: solve unit lower matrix equation L*zz=zr
			zz[0]=zr[0];
			for (i=1;i<=mm1;++i) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				for (j=begin_row;j<=end_row;++j) {
					jj=ijlu[j];
					if (jj<i) zr[i]-=lu[j]*zz[jj];
					else break;
				}
				zz[i]=zr[i];
			}
			
			// backward sweep: solve upper matrix equation U*z=zz
			z[mm1]=zz[mm1]*lu[mm1];
			for (i=mm2;i>=0;i--) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				for (j=end_row;j>=begin_row;j--) {
					jj=ijlu[j];
					if (jj>i) zz[i]-=lu[j]*z[jj];
					else break;
				} 
				z[i]=zz[i]*lu[i];
			}
			
			break; //end (if nb==1) 
			
		case 3:
			
			// forward sweep: solve unit lower matrix equation L*zz=zr
			zz[0] = zr[0];
			zz[1] = zr[1];
			zz[2] = zr[2];
			
			for (i=1;i<=mm1;++i) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb;
				for (j=begin_row;j<=end_row;++j) {
					jj=ijlu[j];
					if (jj<i)
					{
						fasp_blas_smat_mxv_nc3(&(lu[j*nb2]),&(zz[jj*nb]),mult);
						for (ib=0;ib<nb;++ib) zr[ibstart+ib]-=mult[ib];
					}
					else break;
				}
				
				zz[ibstart]   = zr[ibstart];
				zz[ibstart+1] = zr[ibstart+1];
				zz[ibstart+2] = zr[ibstart+2];
			}
			
			// backward sweep: solve upper matrix equation U*z=zz
			ibstart=mm1*nb2;
			ibstart1=mm1*nb;
			fasp_blas_smat_mxv_nc3(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]));
			
			for (i=mm2;i>=0;i--) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb2;
				ibstart1=i*nb;
				for (j=end_row;j>=begin_row;j--) {
					jj=ijlu[j];
					if (jj>i) 
						
					{
						fasp_blas_smat_mxv_nc3(&(lu[j*nb2]),&(z[jj*nb]),mult);
						for (ib=0;ib<nb;++ib) zz[ibstart1+ib]-=mult[ib];
					}
					
					else break;
				} 
				
				fasp_blas_smat_mxv_nc3(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]));
				
			}
			
			break; // end (if nb=3)
			
		case 5:
			
			// forward sweep: solve unit lower matrix equation L*zz=zr
			fasp_array_cp(nb,&(zr[0]),&(zz[0]));
			
			for (i=1;i<=mm1;++i) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb;
				for (j=begin_row;j<=end_row;++j) {
					jj=ijlu[j];
					if (jj<i)
					{
						fasp_blas_smat_mxv_nc5(&(lu[j*nb2]),&(zz[jj*nb]),mult);
						for (ib=0;ib<nb;++ib) zr[ibstart+ib]-=mult[ib];
					}
					else break;
				}
				
				fasp_array_cp(nb,&(zr[ibstart]),&(zz[ibstart]));
			}
			
			// backward sweep: solve upper matrix equation U*z=zz
			ibstart=mm1*nb2;
			ibstart1=mm1*nb;
			fasp_blas_smat_mxv_nc5(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]));
			
			for (i=mm2;i>=0;i--) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb2;
				ibstart1=i*nb;
				for (j=end_row;j>=begin_row;j--) {
					jj=ijlu[j];
					if (jj>i) 
						
					{
						fasp_blas_smat_mxv_nc5(&(lu[j*nb2]),&(z[jj*nb]),mult);
						for (ib=0;ib<nb;++ib) zz[ibstart1+ib]-=mult[ib];
					}
					
					else break;
				} 
				
				fasp_blas_smat_mxv_nc5(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]));
				
			}
			
			break; //end (if nb==5)
			
		case 7:
			
			// forward sweep: solve unit lower matrix equation L*zz=zr
			fasp_array_cp(nb,&(zr[0]),&(zz[0]));
			
			for (i=1;i<=mm1;++i) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb;
				for (j=begin_row;j<=end_row;++j) {
					jj=ijlu[j];
					if (jj<i)
					{
						fasp_blas_smat_mxv_nc7(&(lu[j*nb2]),&(zz[jj*nb]),mult);
						for (ib=0;ib<nb;++ib) zr[ibstart+ib]-=mult[ib];
					}
					else break;
				}
				
				fasp_array_cp(nb,&(zr[ibstart]),&(zz[ibstart]));
			}
			
			// backward sweep: solve upper matrix equation U*z=zz
			ibstart=mm1*nb2;
			ibstart1=mm1*nb;
			fasp_blas_smat_mxv_nc7(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]));
			
			for (i=mm2;i>=0;i--) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb2;
				ibstart1=i*nb;
				for (j=end_row;j>=begin_row;j--) {
					jj=ijlu[j];
					if (jj>i) 
						
					{
						fasp_blas_smat_mxv_nc7(&(lu[j*nb2]),&(z[jj*nb]),mult);
						for (ib=0;ib<nb;++ib) zz[ibstart1+ib]-=mult[ib];
					}
					
					else break;
				} 
				
				fasp_blas_smat_mxv_nc7(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]));
				
			}
			
			break; //end (if nb==7)
			
		default:
			
			// forward sweep: solve unit lower matrix equation L*zz=zr
			fasp_array_cp(nb,&(zr[0]),&(zz[0]));
			
			for (i=1;i<=mm1;++i) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb;
				for (j=begin_row;j<=end_row;++j) {
					jj=ijlu[j];
					if (jj<i)
					{
						fasp_blas_smat_mxv(&(lu[j*nb2]),&(zz[jj*nb]),mult,nb);
						for (ib=0;ib<nb;++ib) zr[ibstart+ib]-=mult[ib];
					}
					else break;
				}
				
				fasp_array_cp(nb,&(zr[ibstart]),&(zz[ibstart]));
			}
			
			// backward sweep: solve upper matrix equation U*z=zz
			ibstart=mm1*nb2;
			ibstart1=mm1*nb;
			fasp_blas_smat_mxv(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]),nb);
			
			for (i=mm2;i>=0;i--) {
				begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
				ibstart=i*nb2;
				ibstart1=i*nb;
				for (j=end_row;j>=begin_row;j--) {
					jj=ijlu[j];
					if (jj>i) 
						
					{
						fasp_blas_smat_mxv(&(lu[j*nb2]),&(z[jj*nb]),mult,nb);
						for (ib=0;ib<nb;++ib) zz[ibstart1+ib]-=mult[ib];
					}
					
					else break;
				} 
				
				fasp_blas_smat_mxv(&(lu[ibstart]),&(zz[ibstart1]),&(z[ibstart1]),nb);
				
			}
			
			break; // end everything else 
	}
	
	return;
}

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
