/*! \file sparse_csr_omp.c
 *  \brief Functions for CSR sparse matrices. 
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_dcsr_getdiag_omp (int n, dCSRmat *A, dvector *diag, int nthreads, int openmp_holds) 
 * \brief Get first n diagonal entries of a CSR matrix A
 *
 * \param  n     an interger (if n=0, then get all diagonal entries)
 * \param *A     pointer to dCSRmat CSR matrix
 * \param *diag  pointer to the diagonal as a dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dcsr_getdiag_omp (int n, 
														dCSRmat *A, 
														dvector *diag, 
														int nthreads, 
														int openmp_holds) 
{
#if FASP_USE_OPENMP
	int i,k,j,ibegin,iend;	
	
	if (n==0) n=MIN(A->row,A->col);
	
	fasp_dvec_alloc(n,diag);
	
	if (n > openmp_holds) {
		int myid;
		int mybegin;
		int myend;
		int stride_i = n/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i, ibegin, iend, k, j) //num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1) myend=mybegin+stride_i;
			else myend=n;
			for (i=mybegin; i<myend; i++)
			{
				ibegin=A->IA[i]; iend=A->IA[i+1];
				for (k=ibegin;k<iend;++k) {
					j=N2C(A->JA[N2C(k)]);
					if ((j-i)==0) {
						diag->val[i] = A->val[N2C(k)]; break;
					} // end if
				} // end for k
			} // end for i
		}
	}
	else {
		for (i=0;i<n;++i) {
			ibegin=A->IA[i]; iend=A->IA[i+1];
			for (k=ibegin;k<iend;++k) {
				j=N2C(A->JA[N2C(k)]);
				if ((j-i)==0) {
					diag->val[i] = A->val[N2C(k)]; break;
				} // end if
			} // end for k
		} // end for i
	}
#endif
}

/**
 * \fn void fasp_dcsr_cp_omp (dCSRmat *A, dCSRmat *B, int nthreads, int openmp_holds)
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param *A   pointer to the dCSRmat matrix
 * \param *B   pointer to the dCSRmat matrix
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dcsr_cp_omp (dCSRmat *A, 
											 dCSRmat *B, 
											 int nthreads, 
											 int openmp_holds)
{
#if FASP_USE_OPENMP
	B->row=A->row;
	B->col=A->col;
	B->nnz=A->nnz;
	
	fasp_iarray_cp_omp(A->row+1, A->IA, B->IA, nthreads, openmp_holds);
	fasp_iarray_cp_omp(A->nnz, A->JA, B->JA, nthreads, openmp_holds);
	fasp_array_cp_omp(A->nnz, A->val, B->val, nthreads, openmp_holds);
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
