/*! \file sparse_block.c
 *  \brief Functions and operation for block sparse matrices. 
 *
 */

#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_dcsr_getblk (dCSRmat *A, int *Is, int *Js, int m, int n, dCSRmat *B)
 * \brief get a sub CSR matrix of A with specified rows and colums
 *
 * \param *A   pointer to dCSRmat CSR matrix
 * \param *B   pointer to dCSRmat CSR matrix
 * \param *Is  pointer to selected rows
 * \param *Js  pointer to selected colums
 * \param m    the number of selected rows
 * \param n    the number of selected colums
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \data 12/25/2010
 */
int fasp_dcsr_getblk (dCSRmat *A, 
											int *Is, 
											int *Js, 
											int m, 
											int n, 
											dCSRmat *B)
{
	int i,j,k,nnz=0;
	int status = SUCCESS;
	int *col_flag;
	
	// create colum flags
	col_flag=(int*)fasp_mem_calloc(A->col,sizeof(int)); 
	
	B->row=m; B->col=n;
	
	B->IA=(int*)fasp_mem_calloc(m+1,sizeof(int));
	
	for (i=0;i<n;++i) col_flag[N2C(Js[i])]=i+1;
	
	// first pass: count nonzeros for sub matrix
	B->IA[0]=0;
	for (i=0;i<m;++i){		
		for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
			j=A->JA[N2C(k)];
			if (col_flag[N2C(j)]>0) nnz++;
		} /* end for k */
		B->IA[i+1]=nnz;
	} /* end for i */
	B->nnz=nnz;
	
	// allocate 
	B->JA=(int*)fasp_mem_calloc(nnz,sizeof(int)); 
	
	B->val=(double*)fasp_mem_calloc(nnz,sizeof(double));
	
	// second pass: copy data to B
	// no need to do the following loop, need to be modified!!  Xiaozhe 
	nnz = 0;
	for (i=0;i<m;++i){		
		for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
			j=A->JA[N2C(k)];
			if (col_flag[N2C(j)]>0) {
				B->JA[nnz]=col_flag[j]-1;
				B->val[nnz]=A->val[N2C(k)];
				nnz++;
			}
		} /* end for k */
	} /* end for i */
	
	fasp_mem_free(col_flag);   
	
	return(status);
}

/**
 * \fn void fasp_dbsr_getblk (dBSRmat *A, int *Is, int *Js, int m, int n, dBSRmat *B)
 * \brief get a sub BSR matrix of A with specified rows and columns. 
 *
 * \param *A   pointer to dBSRmat BSR matrix
 * \param *B   pointer to dBSRmat BSR matrix
 * \param *Is  pointer to selected rows
 * \param *Js  pointer to selected colums
 * \param m    the number of selected rows
 * \param n    the number of selected colums
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \data 12/25/2010
 */
int fasp_dbsr_getblk (dBSRmat *A, 
											int *Is, 
											int *Js, 
											int m, 
											int n, 
											dBSRmat *B)
{
	int i,j,k,nnz=0;
	int status = SUCCESS;
	int *col_flag;
	
	const int nb = A->nb;
	const int nb2=nb*nb;
	
	// create colum flags
	col_flag=(int*)fasp_mem_calloc(A->COL,sizeof(int)); 
	
	B->ROW=m; B->COL=n; B->nb=nb; B->storage_manner=A->storage_manner;  
	
	B->IA=(int*)fasp_mem_calloc(m+1,sizeof(int));
	
	for (i=0;i<n;++i) col_flag[N2C(Js[i])]=i+1;
	
	// first pass: count nonzeros for sub matrix
	B->IA[0]=0;
	for (i=0;i<m;++i){		
		for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
			j=A->JA[N2C(k)];
			if (col_flag[N2C(j)]>0) nnz++;
		} /* end for k */
		B->IA[i+1]=nnz;
	} /* end for i */
	B->NNZ=nnz;
	
	// allocate 
	B->JA=(int*)fasp_mem_calloc(nnz,sizeof(int)); 
	
	B->val=(double*)fasp_mem_calloc(nnz*nb2,sizeof(double));
	
	// second pass: copy data to B
	// no need to do the following loop, need to be modified!!  Xiaozhe 
	nnz = 0;
	for (i=0;i<m;++i){		
		for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
			j=A->JA[N2C(k)];
			if (col_flag[N2C(j)]>0) {
				B->JA[nnz]=col_flag[j]-1;
				memcpy(B->val+nnz*nb2, A->val+N2C(k)*nb2, nb2*sizeof(double));
				nnz++;
			}
		} /* end for k */
	} /* end for i */
	
	fasp_mem_free(col_flag);   
	
	return(status);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
