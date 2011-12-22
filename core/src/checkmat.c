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
 * \fn int fasp_check_diagpos (dCSRmat *A)
 * \brief Check positivity of diagonal entries of a CSR sparse matrix.
 * \param *A pointer to the dCSRmat matrix
 * \return number of negative entries 
 *
 * \author Shuo Zhang
 * \date 03/29/2009
 */
int fasp_check_diagpos (dCSRmat *A)
{
	const int m=A->row, n=A->col;
	dvector diag; fasp_dcsr_getdiag(m,A,&diag);
	
	unsigned int i, num_neg;
	for (num_neg=i=0;i<m;++i) if (diag.val[i]<0) num_neg++;
	
	printf("check_diagpos: nr = %d, nc = %d, nnz = %d\n", m,n,A->nnz);		
	printf("check_diagpos: number of negative diagonal entries = %d\n", num_neg);
	
	fasp_dvec_free(&diag);
	return num_neg;
}

/**
 * \fn int fasp_check_diagzero(dCSRmat *A)
 * \brief Check wether a CSR sparse matrix has diagonal entries that are very close to zero.
 * \param *A pointr to the dCSRmat matrix
 * 
 * \return SUCCESS (0) if no diagonal entry is clase to zero, else ERROR (negative value)
 *
 * \author Shuo Zhang
 * \date 03/29/2009
 */
int fasp_check_diagzero (dCSRmat *A)
{
	const int m = A->row;
	const int *ia=A->IA, *ja=A->JA;
	const double *aj=A->val;
	int i,j,k,begin_row,end_row;
	int status;
	
	for (i=0;i<m;++i) {
		begin_row=ia[i],end_row=ia[i+1];
		for (k=begin_row;k<end_row;++k) {
			j=ja[k];
			if (i==j) {
				if (ABS(aj[k]) < SMALLREAL) {
					printf("Error: diagonal entry (%d,%e) is close to zero!\n",i,aj[k]);
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
 * int fasp_check_diagdom(dCSRmat *A)
 * \brief Check whether a matrix is diagonal dominant.
 * \param *A pointer to the dCSRmat matrix
 * \return print the percentage of the rows which are diagonal dominant and not
 * the number of the rows which are diagonal dominant
 *
 * The routine chechs whether the sparse matrix is diagonal dominant on every row.
 *	 It will print out the percentage of the rows which are diagonal dominant and 
 * which are not; the routine will return the number of the rows which are diagonal 
 * dominant.
 *
 * \author Shuo Zhang
 * \date 03/29/2009
 */
int fasp_check_diagdom (dCSRmat *A)
{
	const int nn  = A->row;
	const int nnz = A->IA[nn]-A->IA[0];
	int i,j,k=0;
	double sum;
	int *rowp;
	
	rowp=(int *)fasp_mem_calloc(nnz,sizeof(int));
	
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
	
	printf("check_diagdom: percentage of the diagonal-dominant rows is %3.2lf%s\n", 
				 100.0*(double)(nn-k)/(double)nn,"%");
	
	fasp_mem_free(rowp);
	
	return k;
}

/**
 * \fn int fasp_check_symm(dCSRmat *A)
 * \brief Check symmetry of a sparse matrix of CSR format.
 * \param *A pointer to the dCSRmat matrix
 * \return 1 and 2 if the structure of the matrix is not symmetric;
 * \return 0 if the structure of the matrix is symmetric,
 * 
 * \note Print the maximal relative difference between matrix and its transpose.
 *
 * \author Shuo Zhang
 * \date 03/29/2009
 */
int fasp_check_symm (dCSRmat *A)
{
	const double symmetry_tol=1.0e-12;
	
	int i,j,mdi,mdj,nnz,nn;
	int *rowp,*rows[2],*cols[2];
	int nns[2],tnizs[2];
	int type=0;
	
	double maxdif,dif;
	double *vals[2];
	
	nn=A->row;
	nnz=A->IA[nn]-A->IA[0];
	
	if (nnz!=A->nnz) {
		printf("check_symm: nnz of the matrix is wrong!!!\n");
		printf("nnz=%d, ia[n]-ia[0]=%d\n",A->nnz,nnz);
		exit(ERROR_WRONG_FILE);
	}
	
	rowp=(int *)fasp_mem_calloc(nnz,sizeof(int));
	
	for (i=0;i<nn;++i) {
		for (j=A->IA[i];j<A->IA[i+1];++j) rowp[N2C(j)]=C2N(i);
	}
	
	rows[0]=(int *)fasp_mem_calloc(nnz,sizeof(int));
	cols[0]=(int *)fasp_mem_calloc(nnz,sizeof(int));
	vals[0]=(double *)fasp_mem_calloc(nnz,sizeof(double));
	
	memcpy(rows[0],rowp,nnz*sizeof(int));
	memcpy(cols[0],A->JA,nnz*sizeof(int));
	memcpy(vals[0],A->val,nnz*sizeof(double));
	
	nns[0]=nn;
	nns[1]=A->col;
	tnizs[0]=nnz;	
	
	rows[1]=(int *)fasp_mem_calloc(nnz,sizeof(int));	
	cols[1]=(int *)fasp_mem_calloc(nnz,sizeof(int));	
	vals[1]=(double *)fasp_mem_calloc(nnz,sizeof(double));
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	memcpy(rows[0],rows[1],nnz*sizeof(int));
	memcpy(cols[0],cols[1],nnz*sizeof(int));
	memcpy(vals[0],vals[1],nnz*sizeof(double));
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
	
	if (type==0) printf("check_symm: matrix is symmetric with max relative difference is %1.3le\n",maxdif);
	if (type==-3) printf("check_symm: matrix is nonsymmetric with max relative difference is %1.3le\n",maxdif);
	if (type==-1) printf("check_symm: matrix has nonsymmetric sp pattern, check the %d-th, %d-th and %d-th rows and cols\n",mdi-1,mdi,mdi+1);
	if (type==-2) printf("check_symm: matrix has nonsymmetric sp pattern, check the %d-th, %d-th and %d-th cols and rows\n",mdj-1,mdj,mdj+1);
	
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
 * \fn int fasp_check_dCSRmat(dCSRmat *A)
 * \brief check whether a dCSRmat is valid or not
 * \param *A pointer to the dCSRmat matrix
 *
 * \author Shuo Zhang
 * \date 03/29/2009
 */
int fasp_check_dCSRmat (dCSRmat *A)
{	
	int i;	
	
	if (A->row != A->col) {
		printf("Error: non-square CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);		
	}
	
	if ((A->nnz==0)|(A->row==0)|(A->col==0)) {
		printf("Error: empty CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	for (i=0;i<A->nnz;++i) {
		if ((N2C(A->JA[i])<0)|(N2C(A->JA[i])-A->col>=0)) {
			printf("Error: wrong CSR matrix format!\n");
			exit(ERROR_DATA_STRUCTURE);
		}
	}
	
	return SUCCESS;
}

/**
 * \fn int fasp_check_iCSRmat(iCSRmat *A)
 * \brief check whether an iCSRmat is valid or not
 * \param *A pointer to the iCSRmat matrix
 *
 * \author Shuo Zhang
 * \date 03/29/2009
 */
int fasp_check_iCSRmat (iCSRmat *A)
{	
	int i;	
	
	if (A->row != A->col) {
		printf("Error: non-square CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);		
	}
	
	if ((A->nnz==0)|(A->row==0)|(A->col==0)) {
		printf("Error: empty CSR matrix!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	for (i=0;i<A->nnz;++i) {
		if ((N2C(A->JA[i])<0)|(N2C(A->JA[i])-A->col>=0)) {
			printf("Error: wrong CSR matrix!\n");
			exit(ERROR_DATA_STRUCTURE);
		}
	}
	
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
