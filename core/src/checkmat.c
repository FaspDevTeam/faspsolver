/*! \file checkmat.c
 *  \brief Check matrix properties
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_check_diagpos (dCSRmat *A)
 *
 * \brief Check positivity of diagonal entries of a CSR sparse matrix.
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \return Number of negative diagonal entries 
 *
 * \author Shuo Zhang
 * \date   03/29/2009
 */
INT fasp_check_diagpos (dCSRmat *A)
{
	const INT      m=A->row;
	unsigned INT   i, num_neg;
	dvector        diag; fasp_dcsr_getdiag(m,A,&diag);
    
#if DEBUG_MODE
	printf("### DEBUG: nr = %d, nc = %d, nnz = %d\n", A->row, A->col, A->nnz);
#endif
    
	for (num_neg=i=0;i<m;++i) if (diag.val[i]<0) num_neg++;
	
	printf("Number of negative diagonal entries = %d\n", num_neg);
	
	fasp_dvec_free(&diag);
    
	return num_neg;
}

/**
 * \fn SHORT fasp_check_diagzero (dCSRmat *A)
 *
 * \brief Check wether a CSR sparse matrix has diagonal entries that are very close to zero.
 *
 * \param A pointr to the dCSRmat matrix
 * 
 * \return SUCCESS if no diagonal entry is clase to zero, else ERROR (negative value)
 *
 * \author Shuo Zhang
 * \date   03/29/2009
 */
SHORT fasp_check_diagzero (dCSRmat *A)
{
	const INT    m = A->row;
	const INT   *ia=A->IA, *ja=A->JA;
	const REAL  *aj=A->val;
	INT          i,j,k,begin_row,end_row;
	SHORT        status;
	
	for (i=0;i<m;++i) {
		begin_row=ia[i],end_row=ia[i+1];
		for (k=begin_row;k<end_row;++k) {
			j=ja[k];
			if (i==j) {
				if (ABS(aj[k]) < SMALLREAL) {
					printf("### ERROR: Diagonal entry (%d,%e) close to zero!\n", i, aj[k]);
					status = ERROR_DATA_ZERODIAG; 
					goto FINISHED;
				}
			}				
		} // end for k
	} // end for i
	
	status = SUCCESS;
	
FINISHED:	
	return status;
}

/**
 * INT fasp_check_diagdom (dCSRmat *A)
 *
 * \brief Check whether a matrix is diagonal dominant.
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \return Number of the rows which are diagonal dominant
 *
 * \note The routine chechs whether the sparse matrix is diagonal dominant on every row.
 *	     It will print out the percentage of the rows which are diagonal dominant and 
 *       which are not; the routine will return the number of the rows which are diagonal 
 *       dominant.
 *
 * \author Shuo Zhang
 * \date   03/29/2009
 */
INT fasp_check_diagdom (dCSRmat *A)
{
	const INT   nn  = A->row;
	const INT   nnz = A->IA[nn]-A->IA[0];
	INT         i,j,k=0;
	REAL        sum;
	
    INT *rowp = (INT *)fasp_mem_calloc(nnz,sizeof(INT));
	
	for (i=0;i<nn;++i) {
		for (j=A->IA[i];j<A->IA[i+1];++j) rowp[j]=i;
	}
	
	for (i=0;i<nn;++i) {
		sum=0.0;
		for (j=A->IA[i];j<A->IA[i+1];++j) {
			if (A->JA[j]==i) sum=sum+A->val[j];
			if (A->JA[j]!=i) sum=sum-fabs(A->val[j]);
		}
		if (sum<-SMALLREAL) ++k;
	}
	
	printf("Percentage of the diagonal-dominant rows is %3.2lf%s\n", 
           100.0*(REAL)(nn-k)/(REAL)nn,"%");
	
	fasp_mem_free(rowp);
	
	return k;
}

/**
 * \fn INT fasp_check_symm (dCSRmat *A)
 *
 * \brief Check symmetry of a sparse matrix of CSR format.
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \return 1 and 2 if the structure of the matrix is not symmetric;
 *         0 if the structure of the matrix is symmetric,
 * 
 * \note Print the maximal relative difference between matrix and its transpose.
 *
 * \author Shuo Zhang
 * \date   03/29/2009
 */
INT fasp_check_symm (dCSRmat *A)
{
	const REAL symmetry_tol=1.0e-12;
	
	INT i,j,mdi,mdj,nnz,nn;
	INT *rowp,*rows[2],*cols[2];
	INT nns[2],tnizs[2];
	INT type=0;
	
	REAL maxdif,dif;
	REAL *vals[2];
	
	nn=A->row;
	nnz=A->IA[nn]-A->IA[0];
	
	if (nnz!=A->nnz) {
		printf("### ERROR: nnz=%d, ia[n]-ia[0]=%d does NOT match!!!\n",A->nnz,nnz);
		exit(ERROR_WRONG_FILE);
	}
	
	rowp=(INT *)fasp_mem_calloc(nnz,sizeof(INT));
	
	for (i=0;i<nn;++i) {
		for (j=A->IA[i];j<A->IA[i+1];++j) rowp[N2C(j)]=C2N(i);
	}
	
	rows[0]=(INT *)fasp_mem_calloc(nnz,sizeof(INT));
	cols[0]=(INT *)fasp_mem_calloc(nnz,sizeof(INT));
	vals[0]=(REAL *)fasp_mem_calloc(nnz,sizeof(REAL));
	
	memcpy(rows[0],rowp,nnz*sizeof(INT));
	memcpy(cols[0],A->JA,nnz*sizeof(INT));
	memcpy(vals[0],A->val,nnz*sizeof(REAL));
	
	nns[0]=nn;
	nns[1]=A->col;
	tnizs[0]=nnz;	
	
	rows[1]=(INT *)fasp_mem_calloc(nnz,sizeof(INT));	
	cols[1]=(INT *)fasp_mem_calloc(nnz,sizeof(INT));	
	vals[1]=(REAL *)fasp_mem_calloc(nnz,sizeof(REAL));
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	memcpy(rows[0],rows[1],nnz*sizeof(INT));
	memcpy(cols[0],cols[1],nnz*sizeof(INT));
	memcpy(vals[0],vals[1],nnz*sizeof(REAL));
	nns[0]=A->col;
	nns[1]=nn;
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	maxdif=0.;
	mdi=0;
	mdj=0;
	for (i=0;i<nnz;++i)
	{
		rows[0][i]=rows[1][i]-rows[0][i];
		if (rows[0][i]!=0) {
			type=-1;
			mdi=rows[1][i];
			break;
		}
		
		cols[0][i]=cols[1][i]-cols[0][i];
		if (cols[0][i]!=0) {
			type=-2;
			mdj=cols[1][i];
			break;
		}
		
		if (fabs(vals[0][i])>SMALLREAL||fabs(vals[1][i])>SMALLREAL) {
			dif=fabs(vals[1][i]-vals[0][i])/(fabs(vals[0][i])+fabs(vals[1][i]));
			if (dif>maxdif) {
				maxdif=dif;
				mdi=rows[0][i];
				mdj=cols[0][i];
			}
		}
	}
	
	if (maxdif>symmetry_tol) type=-3;
	
    switch (type) {
        case 0:
            printf("Matrix is symmetric with max relative difference is %1.3le\n",maxdif);
            break;
        case 3:
            printf("Matrix is nonsymmetric with max relative difference is %1.3le\n",maxdif);
            break;
        case -1:
            printf("Matrix has nonsymmetric pattern, check the %d-th, %d-th and %d-th rows and cols\n",
                   mdi-1,mdi,mdi+1);
            break;
        case -2:
            printf("Matrix has nonsymmetric pattern, check the %d-th, %d-th and %d-th cols and rows\n",
                   mdj-1,mdj,mdj+1);
            break;
        default:
            break;
    }

	fasp_mem_free(rowp);
	fasp_mem_free(rows[1]);
	fasp_mem_free(cols[1]);
	fasp_mem_free(vals[1]);	
	fasp_mem_free(rows[0]);
	fasp_mem_free(cols[0]);
	fasp_mem_free(vals[0]);
	
	return type;
}

/**
 * \fn SHORT fasp_check_dCSRmat (dCSRmat *A)
 *
 * \brief Check whether an dCSRmat matrix is valid or not
 *
 * \param A   Pointer to the matrix in dCSRmat format
 *
 * \author Shuo Zhang
 * \date   03/29/2009
 */
SHORT fasp_check_dCSRmat (dCSRmat *A)
{	
	INT i;	
	
	if (A->row != A->col) {
		printf("### ERROR: Non-square CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);		
	}
	
	if ((A->nnz==0)|(A->row==0)|(A->col==0)) {
		printf("### ERROR: Empty CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	for (i=0;i<A->nnz;++i) {
		if ((N2C(A->JA[i])<0)|(N2C(A->JA[i])-A->col>=0)) {
			printf("### ERROR: Wrong CSR matrix format!\n");
			exit(ERROR_DATA_STRUCTURE);
		}
	}
	
	return SUCCESS;
}

/**
 * \fn SHORT fasp_check_iCSRmat (iCSRmat *A)
 *
 * \brief Check whether an iCSRmat matrix is valid or not
 *
 * \param A   Pointer to the matrix in iCSRmat format
 *
 * \author Shuo Zhang
 * \date   03/29/2009
 */
SHORT fasp_check_iCSRmat (iCSRmat *A)
{	
	INT i;	
	
	if (A->row != A->col) {
		printf("### ERROR: Non-square CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);		
	}
	
	if ((A->nnz==0)|(A->row==0)|(A->col==0)) {
		printf("### ERROR: Empty CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	for (i=0;i<A->nnz;++i) {
		if ((N2C(A->JA[i])<0)|(N2C(A->JA[i])-A->col>=0)) {
			printf("### ERROR: Wrong CSR matrix!\n");
			exit(ERROR_DATA_STRUCTURE);
		}
	}
	
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
