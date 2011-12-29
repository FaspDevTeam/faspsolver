/*! \file smoother.c
 *  \brief Smoothers for sparse matrix in CSR format
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_smoother_gs_cf_omp( dvector *u, dCSRmat *A, dvector *b, int L, int *mark, int order, int nthreads, int openmp_holds )
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \param u     initial guess and the new approximation to the solution obtained after L GS steps
 * \param *A    pointer to stiffness matrix
 * \param *b    pointer to right hand side
 * \param L     number of iterations
 * \param *mark C/F marker array
 * \param order C/F ordering: -1: F-first; 1: C-first
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return      void
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_smoother_dcsr_gs_cf_omp (dvector *u, 
																	 dCSRmat *A, 
																	 dvector *b, 
																	 int L, 
																	 int *mark, 
																	 int order, 
																	 int nthreads, 
																	 int openmp_holds)
{
#if FASP_USE_OPENMP
	const int *ia=A->IA,*ja=A->JA;
	
	int i,j,k,begin_row,end_row;
	int size = b->row;
	
	double *aj=A->val,*bval=b->val,*uval=u->val;
	double t,d=0.0;
	int myid, mybegin, myend;
	
	if (order == -1) // F-point first
	{	
		while (L--) 
		{
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] != 1)
						{
							t = bval[i];
							begin_row = ia[i], end_row = ia[i+1];
							for (k = begin_row; k < end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] != 1)
					{
						t = bval[i];
						begin_row = ia[i], end_row = ia[i+1];
						for (k = begin_row; k < end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
			
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] == 1)
						{
							t = bval[i];
							begin_row = ia[i], end_row = ia[i+1];
							for (k = begin_row; k < end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];					
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] == 1)
					{
						t = bval[i];
						begin_row = ia[i], end_row = ia[i+1];
						for (k = begin_row; k < end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
		} // end while		
	}
	else
	{
		while (L--) 
		{
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] == 1)
						{
							t = bval[i];
							begin_row = ia[i],end_row = ia[i+1]-1;
							for (k = begin_row; k <= end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] == 1)
					{
						t = bval[i];
						begin_row = ia[i],end_row = ia[i+1]-1;
						for (k = begin_row; k <= end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
			
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] != 1)
						{
							t = bval[i];
							begin_row = ia[i],end_row = ia[i+1]-1;
							for (k = begin_row; k <= end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] != 1)
					{
						t = bval[i];
						begin_row = ia[i],end_row = ia[i+1]-1;
						for (k = begin_row; k <= end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
		} // end while	
	}		
#endif
	return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
