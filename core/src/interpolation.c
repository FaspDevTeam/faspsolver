/*
 *  interpolation.c
 *  Form interpolation matrix for AMG
 *
 *------------------------------------------------------
 *  Created by Xuehai Huang on 1/31/09.
 *  Modifed by Chensong Zhang on 04/04/2010.
 *------------------------------------------------------
 *
 */

/*! \file interpolation.c
 *  \brief Interpolation operators for AMG
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

static void interp_RS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param);
static void interp_EM(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param);
static int invden(int nn, double *mat, double *invmat);
static int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask);
static int gentisquare_nomass(dCSRmat *A, int mm, int *Ii, double *ima, int *mask);
static int getinonefull(int **mat, double **matval, int *lengths, int mm, int *Ii, double *ima);
static int orderone(int **mat, double **matval, int *lengths);
static int genintval(dCSRmat *A, int **itmat, double **itmatval, int ittniz, int *isol, int numiso, int nf, int nc);
static int getiteval(dCSRmat *A, dCSRmat *it);

#if FASP_USE_OPENMP
static void interp_RS_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads, int openmp_holds);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_amg_interp (dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
 * \brief Generate interpolation P 
 *
 * \param *A          pointer to the stiffness matrix
 * \param *vertices   pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *P          pointer to the dCSRmat matrix of resulted interpolation
 * \param *param      pointer to AMG parameters
 * \return            SUCCESS or error message
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date 04/04/2010 
 */
int fasp_amg_interp (dCSRmat *A, 
										 ivector *vertices, 
										 dCSRmat *P, 
										 AMG_param *param)
{
	const int interp_type=param->interpolation_type;
	int status = SUCCESS;
	
#if DEBUG_MODE
	printf("fasp_amg_interp ...... [Start]\n");
#endif
	
	/*-- Standard interpolation operator --*/
	switch (interp_type) {
		case 3: // Energy min interpolation in C
			interp_EM(A, vertices, P, param);
			break;
		default: // R-S interpolation
			interp_RS(A, vertices, P, param); 
	}
	
#if DEBUG_MODE
	printf("fasp_amg_interp ...... [Finish]\n");
#endif
	
	if (status<0) exit(status);
	else return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void interp_RS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param)
 * \brief Direct interpolation 
 *
 * \param *A         pointer to the stiffness matrix
 * \param *vertices  pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *Ptr       pointer to the dCSRmat matrix of resulted interpolation
 * \param *param     pointer to AMG parameters
 *
 * Refer to P479, U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. 
 *	          Academic Press Inc., San Diego, CA, 2001. 
 *          With contributions by A. Brandt, P. Oswald and K. St¨uben.
 *
 * \author Xuehai Huang
 * \date 01/31/2009  
 */
static void interp_RS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param)
{
	double epsilon_tr = param->truncation_threshold;
	double amN, amP, apN, apP;
	double alpha, beta, aii=0;
	int *vec = vertices->val;
	int countPplus, diagindex;
	
	unsigned int i,j,k,l,index=0;
	int begin_row, end_row;
	
	/** Generate interpolation P */
	dCSRmat P=fasp_dcsr_create(Ptr->row,Ptr->col,Ptr->nnz);
	
	/** step 1: Find first the structure IA of P */
	memcpy(P.IA,Ptr->IA,(P.row+1)*sizeof(int)); 	//for(i=0;i<=P.row;++i) P.IA[i]=Ptr->IA[i];
	
	/** step 2: Find the structure JA of P */
	memcpy(P.JA,Ptr->JA,P.nnz*sizeof(int)); 	//for(i=0;i<P.nnz;++i) P.JA[i]=Ptr->JA[i];
	
	/** step 3: Fill the data of P */
	for(i=0;i<A->row;++i)
	{
		begin_row=A->IA[i]; end_row=A->IA[i+1]-1;	
		
		for(diagindex=begin_row;diagindex<=end_row;diagindex++) {
			if (A->JA[diagindex]==i) {
				aii=A->val[diagindex];
				break;
			}
		}
		
		if(vec[i]==0)  // if node i is on fine grid 
		{
			amN=0, amP=0, apN=0, apP=0,  countPplus=0;
			
			for(j=begin_row;j<=end_row;++j)
			{
				if(j==diagindex) continue;
				
				for(k=Ptr->IA[i];k<Ptr->IA[i+1];++k) {
					if(Ptr->JA[k]==A->JA[j]) break;
				}
				
				if(A->val[j]>0) {
					apN+=A->val[j];
					if(k<Ptr->IA[i+1]) {
						apP+=A->val[j];
						countPplus++;
					}
				}
				else
				{
					amN+=A->val[j];
					if(k<Ptr->IA[i+1]) {
						amP+=A->val[j];
					}
				}
			} // j
			
			alpha=amN/amP;
			if(countPplus>0) {
				beta=apN/apP;
			}
			else {
				beta=0;
				aii+=apN;
			}
			
			for(j=P.IA[i];j<P.IA[i+1];++j)
			{
				k=P.JA[j];
				for(l=A->IA[i];l<A->IA[i+1];l++)
				{
					if(A->JA[l]==k) break;
				}
				
				if(A->val[l]>0)
				{
					P.val[j]=-beta*A->val[l]/aii;
				}
				else
				{
					P.val[j]=-alpha*A->val[l]/aii;
				}
			}
		}
		else if(vec[i]==2)  // if node i is a special fine node
		{
		}
		else // if node i is on coarse grid 
		{
			P.val[P.IA[i]]=1;
		}
	}	
	
	fasp_mem_free(Ptr->IA);
	fasp_mem_free(Ptr->JA);
	fasp_mem_free(Ptr->val);	

#if 0	// Changed vertices->row to A->row, Edited by Feng Chunsheng 2011/04/11	
	int *CoarseIndex=(int*)fasp_mem_calloc(vertices->row, sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (vertices->row)*sizeof(int);
#endif
	
	index=0;
        for(i=0;i<vertices->row;++i)
#else	
	int *CoarseIndex=(int*)fasp_mem_calloc(A->row, sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (A->row)*sizeof(int);
#endif
	
	index=0;
	for(i=0;i< A->row;++i)    
#endif  // Changed vertices->row to A->row, Edited by Feng Chunsheng 2011/04/11
	{
		if(vec[i]==1)
		{
			CoarseIndex[i]=index;
			index++;
		}
	}
	for(i=0;i<P.IA[P.row];++i)
	{
		j=P.JA[i];
		P.JA[i]=CoarseIndex[j];
	}
	fasp_mem_free(CoarseIndex);
	
	/** Truncation of interpolation */
	double mMin, pMax;
	double mSum, pSum;
	double mTruncedSum, pTruncedSum;
	int mTruncCount, pTruncCount;
	int num_lost=0;
	
	Ptr->val=(double*)fasp_mem_calloc(P.IA[Ptr->row],sizeof(double));
#if CHMEM_MODE
	total_alloc_mem += (P.IA[Ptr->row])*sizeof(double);
#endif
	Ptr->JA=(int*)fasp_mem_calloc(P.IA[Ptr->row],sizeof(int));	
#if CHMEM_MODE
	total_alloc_mem += (P.IA[Ptr->row])*sizeof(int);
#endif
	Ptr->IA=(int*)fasp_mem_calloc(Ptr->row+1, sizeof(int));
#if CHMEM_MODE
	total_alloc_mem += (Ptr->row+1)*sizeof(int);
#endif
	
	int index1=0, index2=0;
	for(i=0;i<P.row;++i)
	{
		mMin=0;
		pMax=0;
		mSum=0;
		pSum=0;
		mTruncedSum=0;
		pTruncedSum=0;
		mTruncCount=0;
		pTruncCount=0;
		
		Ptr->IA[i]-=num_lost;
		
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				mSum+=P.val[j];
				if(P.val[j]<mMin)
				{
					mMin=P.val[j];
				}
			}
			
			if(P.val[j]>0)
			{
				pSum+=P.val[j];
				if(P.val[j]>pMax)
				{
					pMax=P.val[j];
				}
			}
		}
		
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				if(P.val[j]>mMin*epsilon_tr)
				{
					mTruncCount++;
				}
				else
				{
					num_lost--;
				}
			}
			
			if(P.val[j]>0)
			{
				if(P.val[j]<pMax*epsilon_tr)
				{
					pTruncCount++;
				}
				else
				{
					num_lost--;
				}
			}
		}
		
		// step 2: Find the structure JA and fill the data A of Ptr
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				if(!(P.val[j]>mMin*epsilon_tr))
				{
					Ptr->JA[index1]=P.JA[j];
					mTruncedSum+=P.val[j];
					index1++;
				}
			}
			
			if(P.val[j]>0)
			{
				if(!(P.val[j]<pMax*epsilon_tr))
				{
					Ptr->JA[index1]=P.JA[j];
					pTruncedSum+=P.val[j];
					index1++;
				}
			}
		}
		
		// step 3: Fill the data A of Ptr
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				if(!(P.val[j]>mMin*epsilon_tr))
				{
					Ptr->val[index2]=P.val[j]/mTruncedSum*mSum;
					index2++;
				}
			}
			
			if(P.val[j]>0)
			{
				if(!(P.val[j]<pMax*epsilon_tr))
				{
					Ptr->val[index2]=P.val[j]/pTruncedSum*pSum;
					index2++;
				}
			}
		}
	}
	Ptr->IA[P.row]-=num_lost;
	Ptr->nnz=Ptr->IA[Ptr->row];
	
	Ptr->JA=(int*)fasp_mem_realloc(Ptr->JA, Ptr->IA[Ptr->row]*sizeof(int));	
	Ptr->val=(double*)fasp_mem_realloc(Ptr->val, Ptr->IA[Ptr->row]*sizeof(double));
	
	fasp_mem_free(P.IA);
	fasp_mem_free(P.JA);
	fasp_mem_free(P.val);
}

/**
 * \fn static int invden(int nn, double *mat, double *invmat)
 * \brief the routine is to find the inverse of a dense matrix
 *
 * \param nn      scale of the matrix
 * \param mat     the double pointer to the full matrix
 * \param invmat  the double pointer to the full inverse matrix
 * \return        SUCCESS or error message
 * 
 * \note this routine works for symmetric matrix.
 *
 * \author Xuehai Huang
 * \date 04/04/2009 
 */
static int invden(int nn, double *mat, double *invmat)
{
	int indic,i,j;
	int status = SUCCESS;
	int *pivot;
	double *rhs, *sol;
	
	pivot=(int *)fasp_mem_calloc(nn,sizeof(int));	
	rhs=(double *)fasp_mem_calloc(nn,sizeof(double));	
	sol=(double *)fasp_mem_calloc(nn,sizeof(double));
	
	indic=fasp_smat_lu_decomp(mat,pivot,nn);
	
	for (i=0;i<nn;++i) {
		for (j=0;j<nn;++j) rhs[j]=0.;
		rhs[i]=1.;
		fasp_smat_lu_solve(mat,rhs,pivot,sol,nn);
		for (j=0;j<nn;++j) invmat[i*nn+j]=sol[j];
	}
	
	fasp_mem_free(pivot);
	fasp_mem_free(rhs);
	fasp_mem_free(sol);
	
	return status;
}

/**
 * \fn static int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask)
 * \brief Get a local block from a CSR sparse matrix
 *
 * \param *A    pointer to a sparse matrix
 * \param m     number of rows of the local block matrix
 * \param n     number of columns of the local block matrix
 * \param rows  indices to local rows
 * \param cols  indices to local columns
 * \param Aloc  local dense matrix saved as an array
 * \param *mask working array, which should be a negative number initially
 * \return      SUCCESS or error message
 *
 * \author Xuehai Huang
 * \date 04/04/2009 
 */
static int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask)
{
	unsigned int i, j, k, iloc;
	
	for ( i=0; i<m; ++i ) {
		for ( j=0; j<n; ++j ) {
			Aloc[i*n+j] = 0.0; // initialize Aloc
		}			
	}
	
	for ( j=0; j<n; ++j ) {
		mask[N2C(cols[j])] = j; // initialize mask, mask stores C indices 0,1,... 
	}		
	
	for ( i=0; i<m; ++i ) {
		iloc=N2C(rows[i]);
		for ( k=A->IA[iloc]; k<A->IA[iloc+1]; ++k ) {
			j = N2C(A->JA[N2C(k)]);
			if (mask[j]>=0) Aloc[i*n+mask[j]]=A->val[k];
		} /* end for k */
	} /* enf for i */
	
	for ( j=0; j<n; ++j ) mask[N2C(cols[j])] = -1; // re-initialize mask
	
	return SUCCESS;
}

/**
 * \fn static int gentisquare_nomass(dCSRmat *A, int mm, int *Ii, double *ima, int *mask)
 * \brief given the row indices and col indices, to find a block submatrix and get its inverse
 *
 * \param *A     pointer to the whole matrix
 * \param mm     integer of the scale of the submatrix
 * \param Ii     integer to an integer array, to record the indices of row (also col)
 * \param *ima   pointer to the inverse of the full submatrix, the storage is row by row
 * \param *mask  working interger array
 * \return       SUCCESS or error message
 *
 * \author Xuehai Huang
 * \date 04/04/2010 
 */
static int gentisquare_nomass(dCSRmat *A, int mm, int *Ii, double *ima, int *mask)
{
	int indic;
	int status = SUCCESS;
	
	double *ms=(double *)fasp_mem_calloc(mm*mm,sizeof(double));	
	
	indic=get_block(A,mm,mm,Ii,Ii,ms,mask);
	indic=invden(mm,ms,ima);
	
	fasp_mem_free(ms);
	return status;
}

/**
 * \fn static int getinonefull(int **mat, double **matval, int *lengths, int mm, int *Ii, double *ima)
 * \brief to add a small submatrix to a big matrix with respect to its row and cols in the big matrix
 *
 * \param mat      a double pointer pointing to the structure of the matrix
 * \param matval   a double pointer pointing to the values according to the structure
 * \param lengths  a 2d array, the second entry is the lengths of matval
 * \param mm       the number of the rows (also the columns) 
 * \param Ii       an integer pointer to the array to store the relative position of the rows and cols
 * \param ima      the pointer to the full submatrix, the sequence is row by row
 * \return         SUCCESS or error message
 *
 * \author Xuehai Huang
 * \date 04/04/2009 
 */
static int getinonefull(int **mat, double **matval, int *lengths, int mm, int *Ii, double *ima)
{
	int tniz,i,j;
	
	tniz=lengths[1];
	for (i=0;i<mm;++i) {
		for (j=0;j<mm;++j) {		
			mat[0][tniz+i*mm+j]=Ii[i];
			mat[1][tniz+i*mm+j]=Ii[j];
			matval[0][tniz+i*mm+j]=ima[i*mm+j];
		}
	}
	lengths[1]=tniz+mm*mm;
	
	return 1;
}

/**
 * \fn static int orderone(int **mat, double **matval, int *lengths)
 * \brief Order a cluster of entries in a sequence
 *
 * \param **mat     a double pointer to the relative position of the entries
 * \param **matval  a double pointer to the values corresponding to the position
 * \param lengths   int array, to record the number of rows, number of cols and number of nonzerow
 * \return          SUCCESS or error message
 *
 * \author Xuehai Huang
 * \date 04/04/2009 
 */
static int orderone(int **mat, double **matval, int *lengths)
//	lengths[0] for the number of rows
//	lengths[1] for the number of cols
//	lengths[2] for the number of nonzeros
{
	int *rows[2],*cols[2],nns[2],tnizs[2];
	double *vals[2];
	int status = SUCCESS;
	int tniz,i;
	
	nns[0]=lengths[0];
	nns[1]=lengths[1];
	tnizs[0]=lengths[2];
	tniz=lengths[2];
	
	rows[0]=(int *)fasp_mem_calloc(tniz,sizeof(int));	
	cols[0]=(int *)fasp_mem_calloc(tniz,sizeof(int));	
	vals[0]=(double *)fasp_mem_calloc(tniz,sizeof(double));
	
	for (i=0;i<tniz;++i) 
	{
		rows[0][i]=mat[0][i];
		cols[0][i]=mat[1][i];
		vals[0][i]=matval[0][i];
	}
	
	rows[1]=(int *)fasp_mem_calloc(tniz,sizeof(int));	
	cols[1]=(int *)fasp_mem_calloc(tniz,sizeof(int));	
	vals[1]=(double *)fasp_mem_calloc(tniz,sizeof(double));
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	//	all the nonzeros with same col are gathering together
	
	for (i=0;i<tniz;++i)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	//	all the nozeros with same col and row are gathering togheter	
	
	for (i=0;i<tniz;++i)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	
	tniz=tnizs[0];
	for (i=0;i<tniz-1;++i) {
		if (rows[0][i]==rows[0][i+1]&&cols[0][i]==cols[0][i+1]) {
			vals[0][i+1]+=vals[0][i];
			rows[0][i]=nns[0];
			cols[0][i]=nns[1];
		}
	}
	nns[0]=nns[0]+1;
	nns[1]=nns[1]+1;
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	for (i=0;i<tniz;++i)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	for (i=0;i<tniz;++i)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	
	tniz=0;
	for (i=0;i<tnizs[0];++i)
		if (rows[0][i]<nns[0]-1) tniz++;
	
	for (i=0;i<tniz;++i)
	{
		mat[0][i]=rows[0][i];
		mat[1][i]=cols[0][i];
		matval[0][i]=vals[0][i];
	}
	nns[0]=nns[0]-1;
	nns[1]=nns[1]-1;
	lengths[0]=nns[0];
	lengths[1]=nns[1];
	lengths[2]=tniz;
	
	fasp_mem_free(rows[0]);
	fasp_mem_free(rows[1]);
	fasp_mem_free(cols[0]);
	fasp_mem_free(cols[1]);
	fasp_mem_free(vals[0]);
	fasp_mem_free(vals[1]);
	
	return(status);
}


/**
 * \fn static int genintval(dCSRmat *A, int **itmat, double **itmatval, int ittniz, int nf, int nc) 
 * \brief Given the structure of the interpolation, to get the evaluation of the interpolation
 *
 * \param *A         pointer to the dCSRmat matrix
 * \param **itmat    a double integer pointer pointing to the structure of the interpolation
 * \param **itmatval a double double pointer to the evaluation of the interpolation
 * \param ittniz     int, the length of interpolation
 * \param nf         int, the number of fine-level nodes
 * \param nc         int, the number of coarse-level nodes
 * \return           SUCCESS or error message
 *
 * \note 
 *  nf=number fine, nc= n coarse	
 *  Suppose that the structure of the interpolation is known. 
 *  It is N*m matrix, N>m, recorded in itmat.
 *  We record its row index and col index, note that the same col indices gather together
 *  the itma and itmatval have a special data structure
 *  to be exact, the same columns gather together
 *  itmat[0] record the column number, and itmat[1] record the row number.
 */
static int genintval(dCSRmat *A, int **itmat, double **itmatval, int ittniz, int *isol, int numiso, int nf, int nc) 
{
	int *Ii=NULL, *mask=NULL;
	double *ima=NULL, *pex=NULL, **imas=NULL;
	int **mat=NULL;
	double **matval;
	int lengths[3];
	dCSRmat T;
	int tniz;
	dvector sol, rhs;
	
	int mm,sum,i,j,k,status=SUCCESS;
	int *iz,*izs,*izt,*izts;
	
	itsolver_param itparam;
	itparam.print_level    = 0;
	itparam.itsolver_type  = 1;
	itparam.stop_type      = 1;
	itparam.tol            = 1e-3; 
	itparam.maxit          = 100;
	itparam.restart        = 100;
	
	mask=(int *)fasp_mem_calloc(nf,sizeof(int));		
	iz=(int *)fasp_mem_calloc(nc,sizeof(int));		
	izs=(int *)fasp_mem_calloc(nc,sizeof(int));		
	izt=(int *)fasp_mem_calloc(nf,sizeof(int));		
	izts=(int *)fasp_mem_calloc(nf,sizeof(int));
	
	for (i=0;i<nf;++i) mask[i]=-1;
	
	for (i=0;i<nc;++i) iz[i]=0;
	
	for (i=0;i<ittniz;++i) iz[itmat[0][i]]++;
	
	izs[0]=0;
	for (i=1;i<nc;++i) izs[i]=izs[i-1]+iz[i-1];
	
	for (sum=i=0;i<nc;++i) sum+=iz[i]*iz[i];
	
	imas=(double **)fasp_mem_calloc(nc,sizeof(double *));
	
	for (i=0;i<nc;++i) {
		imas[i]=(double *)fasp_mem_calloc(iz[i]*iz[i],sizeof(double));
	}
	
	mat=(int **)fasp_mem_calloc(2,sizeof(int *));
	mat[0]=(int *)fasp_mem_calloc((sum+numiso),sizeof(int));	
	mat[1]=(int *)fasp_mem_calloc((sum+numiso),sizeof(int));		
	matval=(double **)fasp_mem_calloc(1,sizeof(double *));		
	matval[0]=(double *)fasp_mem_calloc(sum+numiso,sizeof(double));
	
	lengths[1]=0;
	
	for (i=0;i<nc;++i) {	
		
		mm=iz[i]; 
		Ii=(int *)fasp_mem_realloc(Ii,mm*sizeof(int));
		
		for (j=0;j<mm;++j) Ii[j]=itmat[1][izs[i]+j];
		
		ima=(double *)fasp_mem_realloc(ima,mm*mm*sizeof(double));
		
		gentisquare_nomass(A,mm,Ii,ima,mask);
		
		getinonefull(mat,matval,lengths,mm,Ii,ima);
		
		for (j=0;j<mm*mm;++j) imas[i][j]=ima[j];
	}
	
	for (i=0;i<numiso;++i) {
		mat[0][sum+i]=isol[i];
		mat[1][sum+i]=isol[i];
		matval[0][sum+i]=1.0;
	}
	
	lengths[0]=nf;
	lengths[2]=lengths[1]+numiso;
	lengths[1]=nf;
	orderone(mat,matval,lengths);
	tniz=lengths[2];
	
	sol.row=nf;
	sol.val=(double*)fasp_mem_calloc(nf,sizeof(double));
	
	for (i=0;i<nf;++i) izt[i]=0;
	
	for (i=0;i<tniz;++i) izt[mat[0][i]]++;
	
	T.IA=(int*)fasp_mem_calloc((nf+1),sizeof(int));
	
	T.row=nf;
	T.col=nf;
	T.nnz=tniz;
	T.IA[0]=0;
	for (i=1;i<nf+1;++i) T.IA[i]=T.IA[i-1]+izt[i-1];
	
	T.JA=(int*)fasp_mem_calloc(tniz,sizeof(int));
	
	for (j=0;j<tniz;++j) T.JA[j]=mat[1][j];
	
	T.val=(double*)fasp_mem_calloc(tniz,sizeof(double));
	
	for (j=0;j<tniz;++j) T.val[j]=matval[0][j];
	
	rhs.val=(double*)fasp_mem_calloc(nf,sizeof(double));
	
	for (i=0;i<nf;++i) rhs.val[i]=1.0;
	rhs.row=nf;
	
	fasp_solver_dcsr_krylov_diag(&T,&rhs,&sol,&itparam);
	
	for (i=0;i<nc;++i)
	{
		mm=iz[i];
		
		ima=(double *)fasp_mem_realloc(ima,mm*mm*sizeof(double));
		
		pex=(double *)fasp_mem_realloc(pex,mm*sizeof(double));
		
		Ii=(int *)fasp_mem_realloc(Ii,mm*sizeof(int));
		
		for (j=0;j<mm;++j) Ii[j]=itmat[1][izs[i]+j];
		
		for (j=0;j<mm*mm;++j) ima[j]=imas[i][j];
		
		for (k=0;k<mm;++k)
		{
			for(pex[k]=j=0;j<mm;++j) pex[k]+=ima[k*mm+j]*sol.val[Ii[j]];
		}
		for (j=0;j<mm;++j) itmatval[0][izs[i]+j]=pex[j];
		
	}
	
	fasp_mem_free(ima);
	fasp_mem_free(pex);
	fasp_mem_free(Ii);
	fasp_mem_free(mask);
	fasp_mem_free(iz);
	fasp_mem_free(izs);
	fasp_mem_free(izt);
	fasp_mem_free(izts);
	fasp_mem_free(mat[0]);
	fasp_mem_free(mat[1]);
	fasp_mem_free(matval[0]);
	
	if (status<0) {
		printf("genintval: Cannot allocate memory!\n");
		exit(status);	
	}
	else {
		return status;
	}
	
}

/**
 * \fn static int getiteval(dCSRmat *A, dCSRmat *it)
 * \brief given a coarsening (in the form of an interpolation operator), inherit the structure,
 *        get new evaluation
 *
 * \param *A    ponter to the dCSRmat matrix
 * \param *it   pointer to the interpolation matrix
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date 10/29/2010  
 */
static int getiteval(dCSRmat *A, dCSRmat *it)
{
	int nf,nc,ittniz;
	int *itmat[2];
	double **itmatval;
	int *rows[2],*cols[2];
	double *vals[2];
	int nns[2],tnizs[2];
	int i,j,indic,numiso;
	int *isol;
	int status = SUCCESS;
	
	nf=A->row;
	nc=it->col;
	ittniz=it->IA[nf]; 
	
	itmat[0]=(int *)fasp_mem_calloc(ittniz,sizeof(int));	
	itmat[1]=(int *)fasp_mem_calloc(ittniz,sizeof(int));	
	itmatval=(double **)fasp_mem_calloc(1,sizeof(double *));	
	itmatval[0]=(double *)fasp_mem_calloc(ittniz,sizeof(double));	
	isol=(int *)fasp_mem_calloc(nf,sizeof(int));
	
	numiso=0;
	for (i=0;i<nf;++i) {
		if (it->IA[i]==it->IA[i+1]) {
			isol[numiso]=i;
			numiso++;
		}
	}
	
	for (i=0;i<nf;++i) {
		for (j=it->IA[i];j<it->IA[i+1];++j) itmat[0][j]=i;
	}
	
	for (j=0;j<ittniz;++j) itmat[1][j]=it->JA[j];
	
	for (j=0;j<ittniz;++j) itmatval[0][j]=it->val[j];
	
	rows[0]=(int *)fasp_mem_calloc(ittniz,sizeof(int));	
	cols[0]=(int *)fasp_mem_calloc(ittniz,sizeof(int));	
	vals[0]=(double *)fasp_mem_calloc(ittniz,sizeof(double));
	
	for (i=0;i<ittniz;++i)
	{
		rows[0][i]=itmat[0][i];
		cols[0][i]=itmat[1][i];
		vals[0][i]=itmat[0][i];
	}
	
	nns[0]=nf;
	nns[1]=nc;
	tnizs[0]=ittniz;
	
	rows[1]=(int *)fasp_mem_calloc(ittniz,sizeof(int));	
	cols[1]=(int *)fasp_mem_calloc(ittniz,sizeof(int));	
	vals[1]=(double *)fasp_mem_calloc(ittniz,sizeof(double));
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	for (i=0;i<ittniz;++i)
	{
		itmat[0][i]=rows[1][i];
		itmat[1][i]=cols[1][i];
		itmatval[0][i]=vals[1][i];
	}
	indic=genintval(A,itmat,itmatval,ittniz,isol,numiso,nf,nc);
	
	for (i=0;i<ittniz;++i)
	{
		rows[0][i]=itmat[0][i];
		cols[0][i]=itmat[1][i];
		vals[0][i]=itmatval[0][i];
	}
	nns[0]=nc;
	nns[1]=nf;
	tnizs[0]=ittniz;
	
	fasp_dcsr_transpose(rows,cols,vals,nns,tnizs);
	
	for (i=0;i<ittniz;++i) it->val[i]=vals[1][i];
	
	fasp_mem_free(isol);
	fasp_mem_free(itmat[0]); 
	fasp_mem_free(itmat[1]);
	fasp_mem_free(itmatval[0]); 
	fasp_mem_free(itmatval);
	fasp_mem_free(rows[0]);
	fasp_mem_free(rows[1]);
	fasp_mem_free(cols[0]);
	fasp_mem_free(cols[1]);
	fasp_mem_free(vals[0]);
	fasp_mem_free(vals[1]);
	
	if (status<0) {
		printf("getiteval: Cannot allocate memory!\n");
		exit(status);	
	}
	else {
		return status;
	}
	
}

/**
 * \fn static void interp_EM(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
 * \brief Energy-min interpolation 
 *
 * \param *A          pointer to the stiffness matrix
 * \param *vertices   pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *P          pointer to the dCSRmat matrix of resulted interpolation
 * \param *param      pointer to AMG parameters 
 *
 * Refer to J. Xu and L. Zikatanov
 *          On An Energy Minimazing Basis in Algebraic Multigrid Methods
 *          Computing and visualization in sciences, 2003
 *
 * \author Shuo Zhang, Xuehai Huang
 * \date 04/04/2010 
 */
static void interp_EM(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
{
	int i, j, index;
	int *vec=vertices->val;	
	int *CoarseIndex=(int*)fasp_mem_calloc(vertices->row,sizeof(int));
	
	for (index=i=0;i<vertices->row;++i) {
		if (vec[i]==1) {
			CoarseIndex[i]=index;
			index++;
		}
	}
	
	for (i=0;i<P->nnz;++i) {
		j=P->JA[i];
		P->JA[i]=CoarseIndex[j];
	}
	
	// clean up memory
	fasp_mem_free(CoarseIndex);		
	
	// main part 
	getiteval(A, P); 
	
	return;
}


/*-----------------------------------omp--------------------------------------*/

/**
 * \fn int fasp_amg_interp_omp (dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param, int nthreads, int openmp_holds)
 * \brief Generate interpolation P 
 *
 * \param *A          pointer to the stiffness matrix
 * \param *vertices   pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *P          pointer to the dCSRmat matrix of resulted interpolation
 * \param *param      pointer to AMG parameters
 * \param nthreads    number of threads
 * \param openmp_holds threshold of parallelization
 * \return            SUCCESS or error message
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_interp_omp (dCSRmat *A, 
												 ivector *vertices, 
												 dCSRmat *P, 
												 AMG_param *param, 
												 int nthreads, 
												 int openmp_holds)
{
	int status = SUCCESS;
#if FASP_USE_OPENMP
	const int interp_type=param->interpolation_type;
	
#if DEBUG_MODE
	printf("fasp_amg_interp ...... [Start]\n");
#endif
	
	/*-- Standard interpolation operator --*/
	switch (interp_type) {
		case 3: // Energy min interpolation in C
			interp_EM(A, vertices, P, param);
			break;
		default: // R-S interpolation
			interp_RS_omp(A, vertices, P, param, nthreads,openmp_holds);
	}
	
#if DEBUG_MODE
	printf("fasp_amg_interp ...... [Finish]\n");
#endif
	
#endif
	if (status<0) exit(status);
	else return status;
}

#if FASP_USE_OPENMP
/**
 * \fn static void interp_RS_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads,, int openmp_holds)
 * \brief Direct interpolation 
 *
 * \param *A         pointer to the stiffness matrix
 * \param *vertices  pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *Ptr       pointer to the dCSRmat matrix of resulted interpolation
 * \param *param     pointer to AMG parameters
 * \param  nthreads  number of threads
 * \param openmp_holds threshold of parallelization
 *
 * Refer to P479, U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. 
 *	          Academic Press Inc., San Diego, CA, 2001. 
 *          With contributions by A. Brandt, P. Oswald and K. St¨uben.
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
static void interp_RS_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads, int openmp_holds)
{
	double epsilon_tr = param->truncation_threshold;
	double amN, amP, apN, apP;
	double alpha, beta, aii=0;
	int *vec = vertices->val;
	int countPplus, diagindex;
	
	unsigned int i,j,k,l,index=0;
	int begin_row, end_row;
	int myid, mybegin, myend;
	
	/** Generate interpolation P */
	dCSRmat P=fasp_dcsr_create(Ptr->row,Ptr->col,Ptr->nnz);
	
	/** step 1: Find first the structure IA of P */
	fasp_iarray_cp_omp(P.row+1, Ptr->IA, P.IA, nthreads,openmp_holds);
	
	/** step 2: Find the structure JA of P */
	fasp_iarray_cp_omp(P.nnz, Ptr->JA, P.JA, nthreads,openmp_holds);
	
	/** step 3: Fill the data of P */
	if (A->row > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i, begin_row, end_row, diagindex, aii, amN, amP, apN, apP, countPplus, j, k, alpha, beta, l) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, A->row, mybegin, myend);
			for (i=mybegin; i<myend; ++i)
			{
				begin_row=A->IA[i]; end_row=A->IA[i+1]-1;	
				
				for(diagindex=begin_row;diagindex<=end_row;diagindex++) {
					if (A->JA[diagindex]==i) {
						aii=A->val[diagindex];
						break;
					}
				}
				
				if(vec[i]==0)  // if node i is on fine grid 
				{
					amN=0, amP=0, apN=0, apP=0,  countPplus=0;
					
					for(j=begin_row;j<=end_row;++j)
					{
						if(j==diagindex) continue;
						
						for(k=Ptr->IA[i];k<Ptr->IA[i+1];++k) {
							if(Ptr->JA[k]==A->JA[j]) break;
						}
						
						if(A->val[j]>0) {
							apN+=A->val[j];
							if(k<Ptr->IA[i+1]) {
								apP+=A->val[j];
								countPplus++;
							}
						}
						else
						{
							amN+=A->val[j];
							if(k<Ptr->IA[i+1]) {
								amP+=A->val[j];
							}
						}
					} // j
					
					alpha=amN/amP;
					if(countPplus>0) {
						beta=apN/apP;
					}
					else {
						beta=0;
						aii+=apN;
					}
					
					for(j=P.IA[i];j<P.IA[i+1];++j)
					{
						k=P.JA[j];
						for(l=A->IA[i];l<A->IA[i+1];l++)
						{
							if(A->JA[l]==k) break;
						}
						
						if(A->val[l]>0)
						{
							P.val[j]=-beta*A->val[l]/aii;
						}
						else
						{
							P.val[j]=-alpha*A->val[l]/aii;
						}
					}
				}
				else if(vec[i]==2)  // if node i is a special fine node
				{
				}
				else // if node i is on coarse grid 
				{
					P.val[P.IA[i]]=1;
				}
			}
		}
	}
	else {
		for(i=0;i<A->row;++i)
		{
			begin_row=A->IA[i]; end_row=A->IA[i+1]-1;	
			
			for(diagindex=begin_row;diagindex<=end_row;diagindex++) {
				if (A->JA[diagindex]==i) {
					aii=A->val[diagindex];
					break;
				}
			}
			
			if(vec[i]==0)  // if node i is on fine grid 
			{
				amN=0, amP=0, apN=0, apP=0,  countPplus=0;
				
				for(j=begin_row;j<=end_row;++j)
				{
					if(j==diagindex) continue;
					
					for(k=Ptr->IA[i];k<Ptr->IA[i+1];++k) {
						if(Ptr->JA[k]==A->JA[j]) break;
					}
					
					if(A->val[j]>0) {
						apN+=A->val[j];
						if(k<Ptr->IA[i+1]) {
							apP+=A->val[j];
							countPplus++;
						}
					}
					else
					{
						amN+=A->val[j];
						if(k<Ptr->IA[i+1]) {
							amP+=A->val[j];
						}
					}
				} // j
				
				alpha=amN/amP;
				if(countPplus>0) {
					beta=apN/apP;
				}
				else {
					beta=0;
					aii+=apN;
				}
				
				for(j=P.IA[i];j<P.IA[i+1];++j)
				{
					k=P.JA[j];
					for(l=A->IA[i];l<A->IA[i+1];l++)
					{
						if(A->JA[l]==k) break;
					}
					
					if(A->val[l]>0)
					{
						P.val[j]=-beta*A->val[l]/aii;
					}
					else
					{
						P.val[j]=-alpha*A->val[l]/aii;
					}
				}
			}
			else if(vec[i]==2)  // if node i is a special fine node
			{
			}
			else // if node i is on coarse grid 
			{
				P.val[P.IA[i]]=1;
			}
		}
	}
	
	fasp_mem_free(Ptr->IA);
	fasp_mem_free(Ptr->JA);
	fasp_mem_free(Ptr->val);	
	

#if 0	// Changed vertices->row to A->row, Edited by Feng Chunsheng 2011/04/11	
	int *CoarseIndex=(int*)fasp_mem_calloc(vertices->row, sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (vertices->row)*sizeof(int);
#endif
	
	index=0;
        for(i=0;i<vertices->row;++i)
#else	
	int *CoarseIndex=(int*)fasp_mem_calloc(A->row, sizeof(int));
	
#if CHMEM_MODE
	total_alloc_mem += (A->row)*sizeof(int);
#endif
	
	index=0;
	for(i=0;i< A->row;++i)    
#endif  // Changed vertices->row to A->row, Edited by Feng Chunsheng 2011/04/11
	{
		if(vec[i]==1)
		{
			CoarseIndex[i]=index;
			index++;
		}
	}
	if (P.IA[P.row] > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i, j) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, P.IA[P.row], mybegin, myend);
			for (i=mybegin; i<myend; ++i)
			{
				j=P.JA[i];
				P.JA[i]=CoarseIndex[j];
			}
		}
	}
	else {
		for(i=0;i<P.IA[P.row];++i)
		{
			j=P.JA[i];
			P.JA[i]=CoarseIndex[j];
		}
	}
	fasp_mem_free(CoarseIndex);
	
	/** Truncation of interpolation */
	double mMin, pMax;
	double mSum, pSum;
	double mTruncedSum, pTruncedSum;
	int mTruncCount, pTruncCount;
	int num_lost=0;
	
	Ptr->val=(double*)fasp_mem_calloc(P.IA[Ptr->row],sizeof(double));
#if CHMEM_MODE
	total_alloc_mem += (P.IA[Ptr->row])*sizeof(double);
#endif
	Ptr->JA=(int*)fasp_mem_calloc(P.IA[Ptr->row],sizeof(int));	
#if CHMEM_MODE
	total_alloc_mem += (P.IA[Ptr->row])*sizeof(int);
#endif
	Ptr->IA=(int*)fasp_mem_calloc(Ptr->row+1, sizeof(int));
#if CHMEM_MODE
	total_alloc_mem += (Ptr->row+1)*sizeof(int);
#endif
	
	int index1=0, index2=0;
	for(i=0;i<P.row;++i)
	{
		mMin=0;
		pMax=0;
		mSum=0;
		pSum=0;
		mTruncedSum=0;
		pTruncedSum=0;
		mTruncCount=0;
		pTruncCount=0;
		
		Ptr->IA[i]-=num_lost;
		
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				mSum+=P.val[j];
				if(P.val[j]<mMin)
				{
					mMin=P.val[j];
				}
			}
			
			if(P.val[j]>0)
			{
				pSum+=P.val[j];
				if(P.val[j]>pMax)
				{
					pMax=P.val[j];
				}
			}
		}
		
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				if(P.val[j]>mMin*epsilon_tr)
				{
					mTruncCount++;
				}
				else
				{
					num_lost--;
				}
			}
			
			if(P.val[j]>0)
			{
				if(P.val[j]<pMax*epsilon_tr)
				{
					pTruncCount++;
				}
				else
				{
					num_lost--;
				}
			}
		}
		
		// step 2: Find the structure JA and fill the data A of Ptr
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				if(!(P.val[j]>mMin*epsilon_tr))
				{
					Ptr->JA[index1]=P.JA[j];
					mTruncedSum+=P.val[j];
					index1++;
				}
			}
			
			if(P.val[j]>0)
			{
				if(!(P.val[j]<pMax*epsilon_tr))
				{
					Ptr->JA[index1]=P.JA[j];
					pTruncedSum+=P.val[j];
					index1++;
				}
			}
		}
		
		// step 3: Fill the data A of Ptr
		for(j=P.IA[i];j<P.IA[i+1];++j)
		{
			if(P.val[j]<0)
			{
				if(!(P.val[j]>mMin*epsilon_tr))
				{
					Ptr->val[index2]=P.val[j]/mTruncedSum*mSum;
					index2++;
				}
			}
			
			if(P.val[j]>0)
			{
				if(!(P.val[j]<pMax*epsilon_tr))
				{
					Ptr->val[index2]=P.val[j]/pTruncedSum*pSum;
					index2++;
				}
			}
		}
	}
	Ptr->IA[P.row]-=num_lost;
	Ptr->nnz=Ptr->IA[Ptr->row];
	
	Ptr->JA=(int*)fasp_mem_realloc(Ptr->JA, Ptr->IA[Ptr->row]*sizeof(int));	
	Ptr->val=(double*)fasp_mem_realloc(Ptr->val, Ptr->IA[Ptr->row]*sizeof(double));
	
	fasp_mem_free(P.IA);
	fasp_mem_free(P.JA);
	fasp_mem_free(P.val);

}

#endif // endif With_OPEN_MP

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
