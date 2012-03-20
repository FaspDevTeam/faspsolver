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
 * \fn SHORT fasp_dcsr_getblk (dCSRmat *A, INT *Is, INT *Js, INT m, INT n, dCSRmat *B)
 *
 * \brief Get a sub CSR matrix of A with specified rows and colums
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param B     Pointer to dCSRmat CSR matrix
 * \param Is    Pointer to selected rows
 * \param Js    Pointer to selected colums
 * \param m     Number of selected rows
 * \param n     Number of selected colums
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   12/25/2010
 */
SHORT fasp_dcsr_getblk (dCSRmat *A, 
                        INT *Is, 
                        INT *Js, 
                        INT m, 
                        INT n, 
                        dCSRmat *B)
{
    
	INT    i,j,k,nnz=0;
	SHORT  status = SUCCESS;
	
	// create colum flags
	INT *col_flag=(int*)fasp_mem_calloc(A->col,sizeof(INT)); 
	
	B->row=m; B->col=n;
	
	B->IA=(int*)fasp_mem_calloc(m+1,sizeof(INT));
	
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
	B->JA=(int*)fasp_mem_calloc(nnz,sizeof(INT)); 
	
	B->val=(REAL*)fasp_mem_calloc(nnz,sizeof(REAL));
	
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
 * \fn SHORT fasp_dbsr_getblk (dBSRmat *A, INT *Is, INT *Js, INT m, INT n, dBSRmat *B)
 *
 * \brief Get a sub BSR matrix of A with specified rows and columns. 
 *
 * \param A     Pointer to dBSRmat BSR matrix
 * \param B     Pointer to dBSRmat BSR matrix
 * \param Is    Pointer to selected rows
 * \param Js    Pointer to selected colums
 * \param m     Number of selected rows
 * \param n     Number of selected colums
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   12/25/2010
 */
SHORT fasp_dbsr_getblk (dBSRmat *A, 
                        INT *Is, 
                        INT *Js, 
                        INT m, 
                        INT n, 
                        dBSRmat *B)
{
	const INT nb = A->nb;
	const INT nb2=nb*nb;
	
	INT i,j,k,nnz=0;
	SHORT status = SUCCESS;
	
	// create colum flags
	INT *col_flag=(int*)fasp_mem_calloc(A->COL,sizeof(INT)); 
	
	B->ROW=m; B->COL=n; B->nb=nb; B->storage_manner=A->storage_manner;  
	
	B->IA=(int*)fasp_mem_calloc(m+1,sizeof(INT));
	
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
	B->JA=(int*)fasp_mem_calloc(nnz,sizeof(INT)); 
	
	B->val=(REAL*)fasp_mem_calloc(nnz*nb2,sizeof(REAL));
	
	// second pass: copy data to B
	// no need to do the following loop, need to be modified!!  Xiaozhe 
	nnz = 0;
	for (i=0;i<m;++i){		
		for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
			j=A->JA[N2C(k)];
			if (col_flag[N2C(j)]>0) {
				B->JA[nnz]=col_flag[j]-1;
				memcpy(B->val+nnz*nb2, A->val+N2C(k)*nb2, nb2*sizeof(REAL));
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
