/*! \file blas_csr.c
 *  \brief BLAS operations for sparse matrices in CSR format.
 *
 *  \note Sparse functions usually contain three runs. The three runs are all the 
 *        same but thy serve different purpose. 
 * 
 *  Example: If you do c=a+b: 
 *			- first do a dry run to find the number of non-zeroes in the result and form ic; 
 *			- allocate space (memory) for jc and form this one; If you only care about a "boolean" result of the addition, you stop here. 
 *			- you call another routine, which uses ic and jc to perform the addition. 
 *
 *		We need to redo these routines later in this way!
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_dcsr_add (dCSRmat *A, const REAL alpha, dCSRmat *B,  
 *                              const REAL beta, dCSRmat *C)
 *
 * \brief compute C = alpha*A + beta*B in CSR format 
 *
 * \param *A      pointer to CSR matrix
 * \param  alpha  real number 
 * \param *B      pointer to CSR matrix
 * \param  beta   real number 
 * \param *C      pointer to CSR matrix
 *
 * \return        SUCCESS if succees, RUN_FAIL if not
 *
 * \author Xiaozhe Hu
 * \date 11/07/2009
 */
INT fasp_blas_dcsr_add (dCSRmat *A, 
                        const REAL alpha, 
                        dCSRmat *B, 
                        const REAL beta, 
                        dCSRmat *C)
{
	INT i,j,k,l;
	INT count=0, added, countrow;
	INT status = SUCCESS;
	
	if (A->row != B->row || A->col != B->col){
		printf("fasp_blas_dcsr_add: The dimension of two matrices does not match!\n");
		status = ERROR_DATA_STRUCTURE;
		goto FINISHED;
	}
	
	if (A == NULL && B == NULL) {
		C->row=0; C->col=0; C->nnz=0;
		status=SUCCESS; goto FINISHED;
	}
	
	if (A->nnz == 0 && B->nnz == 0) {
		C->row=A->row; C->col=A->col; C->nnz=A->nnz;
		status=SUCCESS; goto FINISHED;
	}
	
	// empty matrix A
	if (A->nnz == 0 || A == NULL) {
		fasp_dcsr_alloc(B->row,B->col,B->nnz,C);
		memcpy(C->IA,B->IA,(B->row+1)*sizeof(int));
		memcpy(C->JA,B->JA,(B->nnz)*sizeof(int));
		for (i=0;i<A->nnz;++i) C->val[i]=B->val[i]*beta;
		status=SUCCESS; goto FINISHED;
	}
	
	// empty matrix B
	if (B->nnz == 0 || B == NULL) {
		fasp_dcsr_alloc(A->row,A->col,A->nnz,C);
		memcpy(C->IA,A->IA,(A->row+1)*sizeof(int));
		memcpy(C->JA,A->JA,(A->nnz)*sizeof(int));
		for (i=0;i<A->nnz;++i) C->val[i]=A->val[i]*alpha;
		status=SUCCESS; goto FINISHED;
	}
	
	C->row=A->row; C->col=A->col;
	
	C->IA=(int*)fasp_mem_calloc(C->row+1,sizeof(int));
	
	// allocate work space for C->JA and C->val
	C->JA=(INT *)fasp_mem_calloc(A->nnz+B->nnz,sizeof(int));
	
	C->val=(REAL *)fasp_mem_calloc(A->nnz+B->nnz,sizeof(REAL));
	
	for (i=0; i<A->row; ++i){
		
		countrow = 0;
		for (j=A->IA[i]; j<A->IA[i+1]; ++j){
			C->val[count] = alpha * A->val[N2C(j)];
			C->JA[count] = A->JA[N2C(j)];
			C->IA[i+1]++;
			count++;
			countrow++;
		} // end for j
		
		for (k=B->IA[i]; k<B->IA[i+1]; ++k){
			added = 0;
			
			for (l=C->IA[i]; l<C->IA[i]+countrow+1; l++){
				if (B->JA[N2C(k)] == C->JA[N2C(l)]){
					C->val[N2C(l)] = C->val[N2C(l)] + beta * B->val[N2C(k)];
					added = 1;
					break;
				}
			} // end for l
			
			if (added == 0){
				C->val[count] = beta * B->val[N2C(k)];
				C->JA[count] = B->JA[N2C(k)];
				C->IA[i+1]++;
				count++;
			}	
			
		} // end for k
		
		C->IA[i+1] += C->IA[i];
		
	}
	
	C->nnz = count;
	C->JA  = (INT *)fasp_mem_realloc(C->JA, (count)*sizeof(int));
	C->val = (REAL *)fasp_mem_realloc(C->val, (count)*sizeof(REAL));
	
FINISHED:
	return status;
	
}

/**
 * \fn void fasp_blas_dcsr_axm (dCSRmat *A, const REAL alpha)
 *
 * \brief Multiply a sparse matrix A in CSR format by a scalar alpha.
 *
 * \param *A      pointer to CSR matrix
 * \param  alpha  a real number
 *
 * \return        1 if succees, 0 if fail
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_dcsr_axm (dCSRmat *A, 
                         const REAL alpha)
{
	const INT nnz=A->nnz;
	INT i;
	
	for (i=0; i<nnz; ++i) A->val[i] = A->val[i] * alpha;
}

/**
 * \fn void fasp_blas_dcsr_mxv (dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_dcsr_mxv (dCSRmat *A, 
                         REAL *x, 
                         REAL *y)
{
	const INT   m  = A->row;
	const INT  *ia = A->IA, *ja = A->JA;
	const REAL *aj = A->val;
	INT i, k, begin_row, end_row;		
	register REAL temp;
	
	for (i=0;i<m;++i) {
		temp=0.0; 
		begin_row=ia[i]; end_row=ia[i+1]; 
		for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
		y[i]=temp;
	}
}

/**
 * \fn void fasp_blas_dcsr_mxv_agg (dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = A*x, where the entries of A are all ones.
 *
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 *
 * \author: Xiaozhe Hu
 * \date: 02/22/2011
 */
void fasp_blas_dcsr_mxv_agg (dCSRmat *A, 
                             REAL *x, 
                             REAL *y)
{
	const INT  m  = A->row;
	const INT *ia = A->IA, *ja = A->JA;
	INT i, k, begin_row, end_row;		
	register REAL temp;
	
	for (i=0;i<m;++i) {
		temp=0.0; 
		begin_row=ia[i]; end_row=ia[i+1]; 
		for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
		y[i]=temp;
	}
}

/**
 * \fn void fasp_blas_dcsr_aAxpy (const REAL alpha, dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  real number
 * \param *A     pointer to dCSRmat CSR matrix
 * \param *x     pointer to dvector
 * \param *y     pointer to dvector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_dcsr_aAxpy (const REAL alpha, 
                           dCSRmat *A, 
                           REAL *x,
                           REAL *y)
{
	const INT  m  = A->row;
	const INT *ia = A->IA, *ja = A->JA;
	const REAL *aj = A->val;
	INT i, k, begin_row, end_row;		
	register REAL temp;
	
	if ( alpha == 1.0 ) {
		for (i=0;i<m;++i) {
			temp=0.0; 
			begin_row=ia[i]; end_row=ia[i+1]; 
			for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
			y[i]+=temp;	
		}
	}
	
	else if ( alpha == -1.0 ) {
		for (i=0;i<m;++i) {
			temp=0.0; 
			begin_row=ia[i]; end_row=ia[i+1]; 
			for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
			y[i]-=temp;
		}	
	}
	
	else {
		for (i=0;i<m;++i) {
			temp=0.0; 
			begin_row=ia[i]; end_row=ia[i+1]; 
			for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
			y[i]+=temp*alpha;
		}
	}
	
}

/**
 * \fn void fasp_blas_dcsr_aAxpy_agg (const REAL alpha, dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where the entries of A are all ones
 *
 * \param alpha  a real number
 * \param *A     pointer to dCSRmat CSR matrix
 * \param *x     pointer to dvector
 * \param *y     pointer to dvector
 *
 * \author: Xiaozhe Hu
 * \date: 02/22/2011
 */
void fasp_blas_dcsr_aAxpy_agg (const REAL alpha, 
                               dCSRmat *A, 
                               REAL *x, 
                               REAL *y)
{
	const INT  m  = A->row;
	const INT *ia = A->IA, *ja = A->JA;
	INT i, k, begin_row, end_row;		
	register REAL temp;
	
	if ( alpha == 1.0 ) {
		for (i=0;i<m;++i) {
			temp=0.0; 
			begin_row=ia[i]; end_row=ia[i+1]; 
			for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
			y[i]+=temp;	
		}
	}
	
	else if ( alpha == -1.0 ) {
		for (i=0;i<m;++i) {
			temp=0.0; 
			begin_row=ia[i]; end_row=ia[i+1]; 
			for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
			y[i]-=temp;
		}	
	}
	
	else {
		for (i=0;i<m;++i) {
			temp=0.0; 
			begin_row=ia[i]; end_row=ia[i+1]; 
			for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
			y[i]+=temp*alpha;
		}
	}
	
}

/**
 * \fn REAL fasp_blas_dcsr_vmv (dCSRmat *A, REAL *x, REAL *y) 
 * \brief vector-Matrix-vector multiplication alpha = y'*A*x
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param alpha real number
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
REAL fasp_blas_dcsr_vmv (dCSRmat *A, 
                         REAL *x, 
                         REAL *y)
{
	const INT m=A->row;
	INT i, k, begin_row, end_row;		
	INT *ia=A->IA, *ja=A->JA;
	REAL *aj=A->val;
	register REAL temp;
	register REAL value=0.0;
	
	for (i=0;i<m;++i) {
		temp=0.0; 
		begin_row=ia[i]; end_row=ia[i+1];
		for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
		value+=y[i]*temp;
	}
	
	return value;
}

/**
 * \fn void fasp_blas_dcsr_mxm (dCSRmat *A, dCSRmat *B, dCSRmat *C)
 * \brief Sparse matrix multiplication C=A*B
 *
 * \param *A   pointer to the dCSRmat matrix
 * \param *B   pointer to the dCSRmat matrix
 * \param *C   pointer to dCSRmat matrix equal to A*B
 * \return void
 *
 * \author Xiaozhe Hu
 * \date 11/07/2009
 */
void fasp_blas_dcsr_mxm (dCSRmat *A, 
                         dCSRmat *B, 
                         dCSRmat *C)
{
	// this needs to be changed!!!
	// this fct will be replaced!
	
	INT i,j,k,l,count;
	INT *JD = (INT *)fasp_mem_calloc(B->col,sizeof(int));
	
	C->row=A->row;
	C->col=B->col;
	
	C->val = NULL;
	C->JA  = NULL;	
	C->IA  = (int*)fasp_mem_calloc(C->row+1,sizeof(int));
	
	for (i=0;i<B->col;++i) JD[i]=-1;
	
	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;++i)
	{
		count=0;
		
		for(k=A->IA[i];k<A->IA[i+1];++k)
		{
			for(j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
				for(l=0;l<count;l++) {
					if(JD[l]==B->JA[j]) break;
				}
				
				if(l==count) {
					JD[count]=B->JA[j];
					count++;
				}
			}
		}
		C->IA[i+1]=count;
		for(j=0;j<count;++j) {
			JD[j]=-1;
		}		
	}	
	
	for(i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
	
	// step 2: Find the structure JA of C
	INT countJD;
	
	C->JA=(int*)fasp_mem_calloc(C->IA[C->row],sizeof(int));
	
	for(i=0;i<C->row;++i)
	{
		countJD=0;
		count=C->IA[i];
		for(k=A->IA[i];k<A->IA[i+1];++k)
		{
			for(j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j)
			{
				for(l=0;l<countJD;l++) {
					if(JD[l]==B->JA[j]) break;
				}
				
				if(l==countJD) {
					C->JA[count]=B->JA[j];
					JD[countJD]=B->JA[j];
					count++;
					countJD++;
				}
			}
		}
		
		for(j=0;j<countJD;++j) JD[j]=-1;
	}	
	
	fasp_mem_free(JD);
	
	// step 3: Find the structure A of C
	C->val=(REAL*)fasp_mem_calloc(C->IA[C->row],sizeof(REAL));
	
	for(i=0;i<C->row;++i)
	{
		for(j=C->IA[i];j<C->IA[i+1];++j)
		{
			C->val[j]=0;
			for(k=A->IA[i];k<A->IA[i+1];++k)
			{
				for(l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++)
				{
					if(B->JA[l]==C->JA[j]) {
						C->val[j]+=A->val[k]*B->val[l];
					} // end if
				} // end for l
			} // end for k
		} // end for j
	}	// end for i
	
	C->nnz = C->IA[C->row]-C->IA[0];
	
}

/**
 * \fn void fasp_blas_dcsr_rap (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *B   pointer to dCSRmat matrix equal to R*A*P
 * \return     void
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date 05/10/2010
 */
void fasp_blas_dcsr_rap (dCSRmat *R, 
                         dCSRmat *A, 
                         dCSRmat *P, 
                         dCSRmat *B)
{
	const INT row=R->row, col=P->col;
	unsigned INT nB=A->nnz;
	INT i,i1,j,jj,k,length;	
	INT begin_row,end_row,begin_rowA,end_rowA,begin_rowR,end_rowR;
	INT istart,iistart,count;
	
	REAL *rj=R->val, *aj=A->val, *pj=P->val, *acj;
	INT *ir=R->IA, *ia=A->IA, *ip=P->IA, *iac;
	INT *jr=R->JA, *ja=A->JA, *jp=P->JA, *jac;
	
	INT *index=(INT *)fasp_mem_calloc(A->col,sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (A->col)*sizeof(int);
#endif
	
	INT *iindex=(INT *)fasp_mem_calloc(col,sizeof(int));	
	
#if CHMEM_MODE
	total_alloc_mem += (col)*sizeof(int);
#endif
	
	for (i=0; i<A->col; ++i)  index[i] = -2;
	
	memcpy(iindex,index,col*sizeof(int));
	
	jac=(int*)fasp_mem_calloc(nB,sizeof(int));	
	
#if CHMEM_MODE
	total_alloc_mem += (nB)*sizeof(int);
#endif
	
	iac=(int*)fasp_mem_calloc(row+1,sizeof(int));	
#if CHMEM_MODE
	total_alloc_mem += (row+1)*sizeof(int);
#endif
	
	REAL *temp=(REAL*)fasp_mem_calloc(A->col,sizeof(REAL));
	
#if CHMEM_MODE
	total_alloc_mem += (A->col)*sizeof(REAL);
#endif
	
	iac[0] = 0;
    
	// First loop: form sparsity partern of R*A*P
	for (i=0; i < row; ++i) {		
		// reset istart and length at the begining of each loop
		istart = -1; length = 0; i1 = i+1;
		
		// go across the rows in R
		begin_rowR=ir[i]; end_rowR=ir[i1];
		for (jj=begin_rowR; jj<end_rowR; ++jj){
			j = jr[N2C(jj)];
			// for each column in A
			begin_rowA=ia[j]; end_rowA=ia[j+1];			
			for (k=begin_rowA; k<end_rowA; ++k){
				if (index[N2C(ja[N2C(k)])] == -2){
					index[N2C(ja[N2C(k)])] = istart;
					istart = ja[N2C(k)];
					++length;
				}
			}
		}    
		
		// book-keeping [reseting length and setting iistart]
		count = length; iistart = -1; length = 0;
		
		// use each column that would have resulted from R*A
		for (j=0; j < count; ++j){
			jj = istart;
			istart = index[istart];
			index[N2C(jj)] = -2;
			
			// go across the row of P
			begin_row=ip[jj]; end_row=ip[jj+1];
			for (k=begin_row; k<end_row; ++k){
				// pull out the appropriate columns of P
				if (iindex[N2C(jp[N2C(k)])] == -2){
					iindex[N2C(jp[N2C(k)])] = iistart;
					iistart = jp[N2C(k)];
					++length;
				}
			} // end for k
		} // end for j
		
		// set B->IA
		iac[i1]=iac[i]+length;
		
		if (iac[i1]>nB) {
			nB=nB*2;
			jac=(int*)fasp_mem_realloc(jac, nB*sizeof(int));
		}
		
		// put the correct columns of p into the column list of the products
		begin_row=iac[i]; end_row=iac[i1];
		for (j=begin_row; j<end_row; ++j){
			// put the value in B->JA
			jac[N2C(j)] = iistart;
			// set istart to the next value
			iistart = iindex[N2C(iistart)];
			// set the iindex spot to 0
			iindex[N2C(jac[j])] = -2;
		} // end j
		
	} // end i: First loop
	
	jac=(int*)fasp_mem_realloc(jac,(iac[row])*sizeof(int));
	
	acj=(REAL*)fasp_mem_calloc(iac[row],sizeof(REAL));
	
#if CHMEM_MODE
	total_alloc_mem += (iac[row])*sizeof(REAL);
#endif
	
	INT *BTindex=(int*)fasp_mem_calloc(col,sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (col)*sizeof(int);
#endif
	
	// Second loop: compute entries of R*A*P
	for (i=0; i<row; ++i) {
		i1 = i+1;
		
		// each col of B
		begin_row=iac[i]; end_row=iac[i1];		
		for (j=begin_row; j<end_row; ++j) {
			BTindex[N2C(jac[N2C(j)])]=j;
		}
		
		// reset istart and length at the begining of each loop
		istart = -1; length = 0;
		
		// go across the rows in R
		begin_rowR=ir[i]; end_rowR=ir[i1];		
		for ( jj=begin_rowR; jj<end_rowR; ++jj ){
			j = jr[N2C(jj)];
			
			// for each column in A
			begin_rowA=ia[j]; end_rowA=ia[j+1];					
			for (k=begin_rowA; k<end_rowA; ++k){
				if (index[N2C(ja[N2C(k)])] == -2){
					index[N2C(ja[N2C(k)])] = istart;
					istart = ja[N2C(k)];
					++length;
				}
				temp[N2C(ja[N2C(k)])]+=rj[N2C(jj)]*aj[N2C(k)];
			}
		} 
		
		// book-keeping [reseting length and setting iistart]
		// use each column that would have resulted from R*A
		for (j=0; j<length; ++j) {
			jj = N2C(istart);
			istart = index[N2C(istart)];
			index[N2C(jj)] = -2;
			
			// go across the row of P
			begin_row=ip[jj]; end_row=ip[jj+1];		
			for (k=begin_row; k<end_row; ++k) {
				// pull out the appropriate columns of P
				acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
			}
			temp[jj]=0.0;
		}
		
	} // end for i: Second loop
	
	// setup coarse matrix B
	B->row=row; B->col=col;
	B->IA=iac; B->JA=jac; B->val=acj;	
	B->nnz=B->IA[B->row]-B->IA[0];
	
	fasp_mem_free(temp);
	fasp_mem_free(index);
	fasp_mem_free(iindex);
	fasp_mem_free(BTindex);
}

/**
 * \fn void fasp_blas_dcsr_rap_agg (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 * \brief Triple sparse matrix multiplication B=R*A*P, where the entries of R and P are all ones.
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *B   pointer to dCSRmat matrix equal to R*A*P
 * \return     void
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Xiaozhe Hu
 * \date 02/21/2011
 */
void fasp_blas_dcsr_rap_agg (dCSRmat *R, 
                             dCSRmat *A, 
                             dCSRmat *P, 
                             dCSRmat *B)
{
	const INT row=R->row, col=P->col;
	unsigned INT nB=A->nnz;
	INT i,i1,j,jj,k,length;	
	INT begin_row,end_row,begin_rowA,end_rowA,begin_rowR,end_rowR;
	INT istart,iistart,count;
	
	REAL *aj=A->val, *acj;
	INT *ir=R->IA, *ia=A->IA, *ip=P->IA, *iac;
	INT *jr=R->JA, *ja=A->JA, *jp=P->JA, *jac;
	
	INT *index=(INT *)fasp_mem_calloc(A->col,sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (A->col)*sizeof(int);
#endif
	
	INT *iindex=(INT *)fasp_mem_calloc(col,sizeof(int));	
	
#if CHMEM_MODE
	total_alloc_mem += (col)*sizeof(int);
#endif
	
	for (i=0; i<A->col; ++i)  index[i] = -2;
	
	memcpy(iindex,index,col*sizeof(int));
	
	jac=(int*)fasp_mem_calloc(nB,sizeof(int));	
	
#if CHMEM_MODE
	total_alloc_mem += (nB)*sizeof(int);
#endif
	
	iac=(int*)fasp_mem_calloc(row+1,sizeof(int));	
#if CHMEM_MODE
	total_alloc_mem += (row+1)*sizeof(int);
#endif
	
	REAL *temp=(REAL*)fasp_mem_calloc(A->col,sizeof(REAL));
	
#if CHMEM_MODE
	total_alloc_mem += (A->col)*sizeof(REAL);
#endif
	
	iac[0] = 0;
	
	// First loop: form sparsity partern of R*A*P
	for (i=0; i < row; ++i) {		
		// reset istart and length at the begining of each loop
		istart = -1; length = 0; i1 = i+1;
		
		// go across the rows in R
		begin_rowR=ir[i]; end_rowR=ir[i1];
		for (jj=begin_rowR; jj<end_rowR; ++jj){
			j = jr[N2C(jj)];
			// for each column in A
			begin_rowA=ia[j]; end_rowA=ia[j+1];			
			for (k=begin_rowA; k<end_rowA; ++k){
				if (index[N2C(ja[N2C(k)])] == -2){
					index[N2C(ja[N2C(k)])] = istart;
					istart = ja[N2C(k)];
					++length;
				}
			}
		}    
		
		// book-keeping [reseting length and setting iistart]
		count = length; iistart = -1; length = 0;
		
		// use each column that would have resulted from R*A
		for (j=0; j < count; ++j){
			jj = istart;
			istart = index[istart];
			index[N2C(jj)] = -2;
			
			// go across the row of P
			begin_row=ip[jj]; end_row=ip[jj+1];
			for (k=begin_row; k<end_row; ++k){
				// pull out the appropriate columns of P
				if (iindex[N2C(jp[N2C(k)])] == -2){
					iindex[N2C(jp[N2C(k)])] = iistart;
					iistart = jp[N2C(k)];
					++length;
				}
			} // end for k
		} // end for j
		
		// set B->IA
		iac[i1]=iac[i]+length;
		
		if (iac[i1]>nB) {
			nB=nB*2;
			jac=(int*)fasp_mem_realloc(jac, nB*sizeof(int));
		}
		
		// put the correct columns of p into the column list of the products
		begin_row=iac[i]; end_row=iac[i1];
		for (j=begin_row; j<end_row; ++j){
			// put the value in B->JA
			jac[N2C(j)] = iistart;
			// set istart to the next value
			iistart = iindex[N2C(iistart)];
			// set the iindex spot to 0
			iindex[N2C(jac[j])] = -2;
		} // end j
		
	} // end i: First loop
	
	jac=(int*)fasp_mem_realloc(jac,(iac[row])*sizeof(int));
	
	acj=(REAL*)fasp_mem_calloc(iac[row],sizeof(REAL));
	
#if CHMEM_MODE
	total_alloc_mem += (iac[row])*sizeof(REAL);
#endif
	
	INT *BTindex=(int*)fasp_mem_calloc(col,sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (col)*sizeof(int);
#endif
	
	// Second loop: compute entries of R*A*P
	for (i=0; i<row; ++i) {
		i1 = i+1;
		
		// each col of B
		begin_row=iac[i]; end_row=iac[i1];		
		for (j=begin_row; j<end_row; ++j) {
			BTindex[N2C(jac[N2C(j)])]=j;
		}
		
		// reset istart and length at the begining of each loop
		istart = -1; length = 0;
		
		// go across the rows in R
		begin_rowR=ir[i]; end_rowR=ir[i1];		
		for ( jj=begin_rowR; jj<end_rowR; ++jj ){
			j = jr[N2C(jj)];
			
			// for each column in A
			begin_rowA=ia[j]; end_rowA=ia[j+1];					
			for (k=begin_rowA; k<end_rowA; ++k){
				if (index[N2C(ja[N2C(k)])] == -2){
					index[N2C(ja[N2C(k)])] = istart;
					istart = ja[N2C(k)];
					++length;
				}
				temp[N2C(ja[N2C(k)])]+=aj[N2C(k)];
			}
		} 
		
		// book-keeping [reseting length and setting iistart]
		// use each column that would have resulted from R*A
		for (j=0; j<length; ++j) {
			jj = N2C(istart);
			istart = index[N2C(istart)];
			index[N2C(jj)] = -2;
			
			// go across the row of P
			begin_row=ip[jj]; end_row=ip[jj+1];		
			for (k=begin_row; k<end_row; ++k) {
				// pull out the appropriate columns of P
				acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj];
			}
			temp[jj]=0.0;
		}
		
	} // end for i: Second loop
	
	// setup coarse matrix B
	B->row=row; B->col=col;
	B->IA=iac; B->JA=jac; B->val=acj;	
	B->nnz=B->IA[B->row]-B->IA[0];
	
	fasp_mem_free(temp);
	fasp_mem_free(index);
	fasp_mem_free(iindex);
	fasp_mem_free(BTindex);
}

/**
 * \fn dCSRmat fasp_blas_dcsr_ptap (dCSRmat *A, dCSRmat *P)
 * \brief Triple sparse matrix multiplication B=P'*A*P
 *
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \return *B  pointer to dCSRmat matrix equal to P'*A*P
 *
 * Driver to compute triple matrix product P'*A*P using ltz CSR format. 
 * In ltx format: ia[0]=1, ja[0] and a[0] are used as usual. When called 
 * from Fortran, ia[0], ja[0] and a[0] will be just ia(1),ja(1),a(1).
 * For the indices, 
 *			ia_ltz[k] = ia_usual[k]+1, 
 *			ja_ltz[k] = ja_usual[k]+1,
 *			 a_ltz[k] =  a_usual[k].
 *
 * \author Ludmil Zikatanov, Chensong Zhang
 * \date 05/10/2010
 */
void fasp_blas_dcsr_ptap (dCSRmat *Pt,
                          dCSRmat *A, 
                          dCSRmat *P, 
                          dCSRmat *Ac)
{
	INT nc=Pt->row, n=Pt->col, nnzP=P->nnz, nnzA=A->nnz, maxrpout;
	INT i;
	
	// shift A from usual to ltz format
	for (i=0;i<=n;++i)   { A->IA[i]++; P->IA[i]++; }
	for (i=0;i<nnzA;++i) { A->JA[i]++; }
	for (i=0;i<=nc;++i)  { Pt->IA[i]++; }
	for (i=0;i<nnzP;++i) { P->JA[i]++; Pt->JA[i]++; }
	
	// compute P' A P
	dCSRmat PtAP = fasp_blas_dcsr_rap2(Pt->IA,Pt->JA,Pt->val,A->IA,A->JA,A->val, \
                                       Pt->IA,Pt->JA,Pt->val,n,nc,&maxrpout,     \
                                       P->IA,P->JA);
	
	Ac->row = PtAP.row;
	Ac->col = PtAP.col;
	Ac->nnz = PtAP.nnz;
	Ac->IA  = PtAP.IA;
	Ac->JA  = PtAP.JA;
	Ac->val = PtAP.val;
	
	// shift A back from ltz format
	for (i=0;i<=Ac->row;++i) Ac->IA[i]--;
	for (i=0;i< Ac->nnz;++i) Ac->JA[i]--;
	for (i=0;i<=n;++i)   { A->IA[i]--; P->IA[i]--; }
	for (i=0;i<nnzA;++i) { A->JA[i]--; }
	for (i=0;i<=nc;++i)  { Pt->IA[i]--; }
	for (i=0;i<nnzP;++i) { P->JA[i]--; Pt->JA[i]--; }
	
	return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
