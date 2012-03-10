/*! \file interpolation_omp.c
 *  \brief Interpolation operators for AMG
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if FASP_USE_OPENMP
static void interp_RS_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads, int openmp_holds);
static void fasp_get_nbl_nbr_ysk_omp(dCSRmat *A, int *nbl_ptr, int *nbr_ptr, int nthreads, int openmp_holds);
static void fasp_mod_coarse_index_omp( int nrows, int *CoarseIndex, int nthreads, int openmp_holds );
int fasp_BinarySearch(int *list, int value, int list_length);
static void fasp_get_icor_ysk_omp(int nrows, int ncols, int *CoarseIndex, int nbl_ysk, int nbr_ysk, int *CF_marker, int *icor_ysk, int nthreads, int openmp_holds);
static void interp_RS1_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int *icor_ysk, int nthreads, int openmp_holds);
static void interp_RS2_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads, int openmp_holds);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn int fasp_amg_interp_omp (dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param, int nthreads, int openmp_holds)
 * \brief Generate interpolation P 
 *
 * \param A          pointer to the stiffness matrix
 * \param vertices   pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param P          pointer to the dCSRmat matrix of resulted interpolation
 * \param param      pointer to AMG parameters
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

/**
 * \fn int fasp_amg_interp1_omp (dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param, int nthreads, int openmp_holds)
 * \brief Generate interpolation P 
 *
 * \param A          pointer to the stiffness matrix
 * \param vertices   pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param P          pointer to the dCSRmat matrix of resulted interpolation
 * \param param      pointer to AMG parameters
 * \param nthreads    number of threads
 * \param openmp_holds threshold of parallelization
 * \return            SUCCESS or error message
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_interp1_omp (dCSRmat *A, 
												 ivector *vertices, 
												 dCSRmat *P, 
												 AMG_param *param, 
												 int *icor_ysk,
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
			interp_RS1_omp(A, vertices, P, param, icor_ysk, nthreads, openmp_holds);
	}
	
#if DEBUG_MODE
	printf("fasp_amg_interp ...... [Finish]\n");
#endif
	
#endif
	if (status<0) exit(status);
	else return status;
}

/**
 * \fn int fasp_amg_interp2_omp (dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param, int nthreads, int openmp_holds)
 * \brief Generate interpolation P 
 *
 * \param A          pointer to the stiffness matrix
 * \param vertices   pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param P          pointer to the dCSRmat matrix of resulted interpolation
 * \param param      pointer to AMG parameters
 * \param nthreads    number of threads
 * \param openmp_holds threshold of parallelization
 * \return            SUCCESS or error message
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_interp2_omp (dCSRmat *A, 
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
			interp_RS2_omp(A, vertices, P, param, nthreads,openmp_holds);
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
 * \param A         pointer to the stiffness matrix
 * \param vertices  pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param Ptr       pointer to the dCSRmat matrix of resulted interpolation
 * \param param     pointer to AMG parameters
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
	int myid;
	int mybegin;
	int myend;
	int stride_i;
	
	/** Generate interpolation P */
	dCSRmat P=fasp_dcsr_create(Ptr->row,Ptr->col,Ptr->nnz);
	
	/** step 1: Find first the structure IA of P */
	fasp_iarray_cp_omp(P.row+1, Ptr->IA, P.IA, nthreads,openmp_holds);
	
	/** step 2: Find the structure JA of P */
	fasp_iarray_cp_omp(P.nnz, Ptr->JA, P.JA, nthreads,openmp_holds);
	
	/** step 3: Fill the data of P */
	if (A->row > openmp_holds) {
		stride_i = A->row/nthreads;
#pragma omp parallel private(myid,mybegin,myend,i,begin_row,end_row,diagindex,aii,amN,amP,apN,apP,countPplus,j,k,alpha,beta,l) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1) myend = mybegin+stride_i;
			else myend = A->row;
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
#if 1
	int *CoarseIndex=(int*)fasp_mem_calloc(A->row, sizeof(int));
#else
	int *CoarseIndex=(int*)fasp_mem_calloc(vertices->row, sizeof(int));
#endif
	
#if CHMEM_MODE
#if 1
	total_alloc_mem += (A->row)*sizeof(int);
#else
	total_alloc_mem += (vertices->row)*sizeof(int);
#endif
#endif
	
	index=0;
#if 1
	for(i=0;i<A->row;++i)
#else
	for(i=0;i<vertices->row;++i)
#endif
	{
		if(vec[i]==1)
		{
			CoarseIndex[i]=index;
			index++;
		}
	}
	if (P.IA[P.row] > openmp_holds) {
		stride_i = P.IA[P.row]/nthreads;
#pragma omp parallel private(myid,mybegin,myend,i,j) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1) myend = mybegin+stride_i;
			else myend = P.IA[P.row];
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

/**
 * \fn static void fasp_get_nbl_nbr_ysk(dCSRmat *A, int *nbl_ptr, int *nbr_ptr, int nthreads, int openmp_holds)
 * \brief Get the bandwidth of the matrix
 *
 * \param A         pointer to the stiffness matrix
 * \param nbl_ptr   the left bandwidth
 * \param nbr_ptr   the right bandwidth
 * \param  nthreads  number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
static void fasp_get_nbl_nbr_ysk_omp(dCSRmat *A, int *nbl_ptr, int *nbr_ptr, int nthreads, int openmp_holds)
{
	int *IA = A->IA;
	int *JA = A->JA;
	int myid, mybegin, myend;
	int max_l, max_r;
	int i, end_row_A, j;
	
	int *max_left_right = (int *)fasp_mem_calloc(2*nthreads, sizeof(int));
#pragma omp parallel for private(myid,mybegin,myend,max_l,max_r,i,end_row_A,j)
	for (myid = 0; myid < nthreads; myid ++)
	{
		FASP_GET_START_END(myid, nthreads, A->row, mybegin, myend);
		max_l = 0;
		max_r = 0;
		for (i = mybegin; i < myend; i ++) {
			end_row_A = IA[i+1];
			for (j = IA[i]; j < end_row_A; j ++) {
				max_l = MAX(i-JA[j], max_l);
				max_r = MAX(JA[j]-i, max_r);
			}
		}
		max_left_right[myid*2] = max_l;
		max_left_right[myid*2+1] = max_r;
	}
	max_l = max_left_right[0];
	max_r = max_left_right[1];
	for (i = 1; i < nthreads; i ++) {
		max_l = MAX(max_l, max_left_right[i*2]);
        max_r = MAX(max_r, max_left_right[i*2+1]);
    }
    fasp_mem_free(max_left_right);
   *nbl_ptr = max_l;
   *nbr_ptr = max_r;
}

/**
 * \fn static void fasp_mod_coarse_index_omp( int nrows, int *CoarseIndex, int nthreads, int openmp_holds )
 * \brief Modify CoarseIndex
 *
 * \param nrows         length of CoarseIndex
 * \param CoarseIndex  pointer to index of coarse points
 * \param  nthreads     number of threads
 * \param openmp_holds  threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
static void fasp_mod_coarse_index_omp( int nrows, int *CoarseIndex, int nthreads, int openmp_holds )
{
	int myid, mybegin, myend, i;
	
#pragma omp parallel for private(myid,mybegin,myend,i)
	for (myid = 0; myid < nthreads; myid ++)
	{
		FASP_GET_START_END(myid, nthreads, nrows, mybegin, myend);
		if (myid == 0)
		{
			mybegin ++;
		}
		for (i = mybegin; i < myend; i ++)
		{
			if (CoarseIndex[i] < CoarseIndex[i-1])
			{
				CoarseIndex[i] = CoarseIndex[i-1];
			}
		}
	}
}

/**
 * \fn int fasp_BinarySearch(int *list, int value, int list_length)
 * \brief Binary Search
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_BinarySearch(int *list, int value, int list_length)
{
	int not_found = 1;
	int low, high, m;
	
	low = 0;
	high = list_length - 1;
	while (not_found && low <= high)
	{
		m = (low + high) / 2;
		if (value < list[m])
		{
			high = m - 1;
		}
		else if (value > list[m])
		{
			low = m + 1;
		}
		else
		{
			not_found = 0;
			return m;
		}
	}
	
	return -1;
}

/**
 * \fn static void fasp_get_icor_ysk_omp
 * \brief Get Icor_ysk
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
static void fasp_get_icor_ysk_omp(int nrows, int ncols, int *CoarseIndex, int nbl_ysk, int nbr_ysk, int *CF_marker, int *icor_ysk, int nthreads, int openmp_holds)
{
	int myid, FiveMyid, mybegin, myend, min_A, max_A, i, first_f_node, min_P, max_P, myend_minus_one;
	int lengthAA = 0, lengthPP = 0;
	
#pragma omp parallel for private(myid,FiveMyid,mybegin,myend,min_A,max_A,i,first_f_node,min_P,max_P,myend_minus_one) reduction(+: lengthAA,lengthPP)
	for (myid = 0; myid < nthreads; myid ++)
	{
		FiveMyid = myid * 5;
		FASP_GET_START_END(myid, nthreads, ncols, mybegin, myend);
		icor_ysk[FiveMyid] = mybegin;
		if (mybegin == myend) {
			lengthAA = 0;
			lengthPP = 0;
			icor_ysk[FiveMyid+1] = 0;
			icor_ysk[FiveMyid+3] = 0;
		}
		else {
			first_f_node = fasp_BinarySearch(CoarseIndex, mybegin, nrows);
			for (i = first_f_node; i > -1; i --) {
				if (CoarseIndex[i] != mybegin) {
					break;
				}
			}
			min_A = i + 1;
			min_A = MAX(0, min_A-2*nbl_ysk);
			myend_minus_one = myend - 1;
			first_f_node = fasp_BinarySearch(CoarseIndex, myend_minus_one, nrows);
			for (i = first_f_node; i > -1; i --) {
				if (CoarseIndex[i] != myend_minus_one) {
					max_A = i;
					break;
				}
			}
			max_A = MIN(nrows, max_A+2*nbr_ysk+1);
			lengthAA = max_A - min_A + 2;
			icor_ysk[FiveMyid+1] = lengthAA;
			icor_ysk[FiveMyid+2] = min_A;
			for (i = min_A; i >= 0; i --) {
				if (CF_marker[i] == 0) {
					first_f_node = i;
					break;
				}
			}
			if (i == -1) {
				min_P = 0;
			}
			else {
				first_f_node -= nbl_ysk;
				if (first_f_node <= 0) {
					min_P = 0;
				}
				else {
					for (i = first_f_node; i >= 0; i --) {
						if (CF_marker[i] == 1) {
							min_P = CoarseIndex[i];
							break;
						}
					}
					if (i == -1) {
						min_P = 0;
					}
				}
			}
			for (i = max_A-1; i < nrows; i ++) {
				if (CF_marker[i] == 0) {
					first_f_node = i;
					break;
				}
			}
			if (i == nrows) {
				max_P = ncols;
			}
			else {
				first_f_node += nbr_ysk;
				if (first_f_node >= nrows) {
					max_P = ncols;
				}
				else {
					for (i = first_f_node; i < nrows; i ++) {
						if (CF_marker[i] == 1) {
							max_P = CoarseIndex[i] + 1;
							break;
						}
					}
					if (i == nrows) {
						max_P = ncols;
					}
				}
			}
			lengthPP = max_P - min_P + 2;
			icor_ysk[FiveMyid+3] = lengthPP;
			icor_ysk[FiveMyid+4] = min_P;
		}
	}
	icor_ysk[5*nthreads] = lengthAA;
	icor_ysk[5*nthreads+1] = lengthPP;
}

/**
 * \fn static void interp_RS1_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads,, int openmp_holds)
 * \brief Direct interpolation 
 *
 * \param A         pointer to the stiffness matrix
 * \param vertices  pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param Ptr       pointer to the dCSRmat matrix of resulted interpolation
 * \param param     pointer to AMG parameters
 * \param  nthreads  number of threads
 * \param openmp_holds threshold of parallelization
 *
 * Refer to P479, U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. 
 *	          Academic Press Inc., San Diego, CA, 2001. 
 *          With contributions by A. Brandt, P. Oswald and K. St¨uben.
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012 Modified by FENG Chunsheng for omp gcc
 */
static void interp_RS1_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int *icor_ysk, int nthreads, int openmp_holds)
{
	double epsilon_tr = param->truncation_threshold;
	double amN, amP, apN, apP;
	double alpha, beta, aii=0;
	int *vec = vertices->val;
	int countPplus, diagindex;
	
	unsigned int i,j,k,l,index=0;
	int begin_row, end_row;
	int myid;
	int mybegin;
	int myend;
	int stride_i;
	int *indexs = NULL;
	int shift;
	int nbl_ysk, nbr_ysk;
	
	/** Generate interpolation P */
	dCSRmat P=fasp_dcsr_create(Ptr->row,Ptr->col,Ptr->nnz);
	
	/** step 1: Find first the structure IA of P */
	fasp_iarray_cp_omp(P.row+1, Ptr->IA, P.IA, nthreads,openmp_holds);
	
	/** step 2: Find the structure JA of P */
	fasp_iarray_cp_omp(P.nnz, Ptr->JA, P.JA, nthreads,openmp_holds);
	
	/** step 3: Fill the data of P */
	if (A->row > openmp_holds) {
#pragma omp parallel for private(myid,mybegin,myend,i,begin_row,end_row,diagindex,aii,amN,amP,apN,apP,countPplus,j,k,alpha,beta,l)
			for (myid = 0; myid < nthreads; myid++ )
			{
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
	
	int *CoarseIndex;

	CoarseIndex=(int*)fasp_mem_calloc(A->row, sizeof(int));

#if CHMEM_MODE
	total_alloc_mem += (A->row)*sizeof(int);
#endif
	
	// The following is one of OPTIMAL parts ...0802...
	// Generate CoarseIndex in parallel
	if (A->row > openmp_holds) {
#pragma omp master
	{	
		indexs = (int *)fasp_mem_calloc(nthreads, sizeof(int));
	}
#pragma omp parallel for private(myid, mybegin, myend, index, i)
		for (myid = 0; myid < nthreads; myid ++)
		{
			FASP_GET_START_END(myid, nthreads, A->row, mybegin, myend);
			index = 0;
			for (i=mybegin;i<myend;++i) {
				if(vec[i]==1)
				{
					CoarseIndex[i]=index;
					index++;
				}
			}
			indexs[myid] = index;
		}
		for (i = 1; i < nthreads; i ++) {
			indexs[i] += indexs[i-1];
		}
#pragma omp parallel for private(myid, mybegin, myend, shift, i)
		for (myid = 0; myid < nthreads; myid ++)
		{
			FASP_GET_START_END(myid, nthreads, A->row, mybegin, myend);
			shift = 0;
			if (myid > 0) {
				shift = indexs[myid-1];
			}
			for (i=mybegin;i<myend;++i) {
				if(vec[i]==1)
				{
					CoarseIndex[i] += shift;
				}
			}
		}
		fasp_mem_free(indexs);
	}
	else {
		index=0;
		for(i=0;i<A->row;++i) {
			if(vec[i]==1)
			{
				CoarseIndex[i]=index;
				index++;
			}
		}
	}
	
	if (P.IA[P.row] > openmp_holds) {
#pragma omp for parallel private(myid, mybegin,myend,i,j) 
			for (myid = 0; myid < nthreads; myid++ )
			{
				FASP_GET_START_END(myid, nthreads,P.IA[P.row], mybegin, myend);
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
	if (Ptr->col > openmp_holds) {
		// The following is another OPTIMAL part ...0802...
	////////////////////////////////////////////////////////////////////////////////
		fasp_get_nbl_nbr_ysk_omp(A, &nbl_ysk, &nbr_ysk, nthreads, openmp_holds);
		fasp_mod_coarse_index_omp(A->row, CoarseIndex, nthreads, openmp_holds);
		fasp_get_icor_ysk_omp(A->row, Ptr->col, CoarseIndex, nbl_ysk, nbr_ysk, vec, icor_ysk, nthreads, openmp_holds);
	////////////////////////////////////////////////////////////////////////////////
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
 * \fn static void interp_RS2_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads,, int openmp_holds)
 * \brief Direct interpolation 
 *
 * \param A         pointer to the stiffness matrix
 * \param vertices  pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param Ptr       pointer to the dCSRmat matrix of resulted interpolation
 * \param param     pointer to AMG parameters
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
static void interp_RS2_omp(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param, int nthreads, int openmp_holds)
{
	double epsilon_tr = param->truncation_threshold;
	double amN, amP, apN, apP;
	double alpha, beta, aii=0;
	int *vec = vertices->val;
	int countPplus, diagindex;
	
	unsigned int i,j,k,l,index=0;
	int begin_row, end_row;
	int myid;
	int mybegin;
	int myend;
	int stride_i;
	int *indexs = NULL;
	int shift;
	
	/** Generate interpolation P */
	dCSRmat P=fasp_dcsr_create(Ptr->row,Ptr->col,Ptr->nnz);
	
	/** step 1: Find first the structure IA of P */
	fasp_iarray_cp_omp(P.row+1, Ptr->IA, P.IA, nthreads,openmp_holds);
	
	/** step 2: Find the structure JA of P */
	fasp_iarray_cp_omp(P.nnz, Ptr->JA, P.JA, nthreads,openmp_holds);
	
	/** step 3: Fill the data of P */
	if (A->row > openmp_holds) {
		stride_i = A->row/nthreads;
#pragma omp parallel private(myid,mybegin,myend,i,begin_row,end_row,diagindex,aii,amN,amP,apN,apP,countPplus,j,k,alpha,beta,l) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1) myend = mybegin+stride_i;
			else myend = A->row;
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
#if 1
	int *CoarseIndex=(int*)fasp_mem_calloc(A->row, sizeof(int));
#else
	int *CoarseIndex=(int*)fasp_mem_calloc(vertices->row, sizeof(int));
#endif
	
#if CHMEM_MODE
#if 1
	total_alloc_mem += (A->row)*sizeof(int);
#else
	total_alloc_mem += (vertices->row)*sizeof(int);
#endif
#endif
	
	index=0;
#if 0
#if 1
	for(i=0;i<A->row;++i)
#else
	for(i=0;i<vertices->row;++i)
#endif
	{
		if(vec[i]==1)
		{
			CoarseIndex[i]=index;
			index++;
		}
	}
#else
	// The following is one of OPTIMAL parts ...0802...
	// Generate CoarseIndex in parallel
	if (A->row > openmp_holds) {
		indexs = (int *)fasp_mem_calloc(nthreads, sizeof(int));
#pragma omp parallel for private(myid, mybegin, myend, index, i)
		for (myid = 0; myid < nthreads; myid ++)
		{
			FASP_GET_START_END(myid, nthreads, A->row, mybegin, myend);
			index = 0;
			for (i=mybegin;i<myend;++i) {
				if(vec[i]==1)
				{
					CoarseIndex[i]=index;
					index++;
				}
			}
			indexs[myid] = index;
		}
		for (i = 1; i < nthreads; i ++) {
			indexs[i] += indexs[i-1];
		}
#pragma omp parallel for private(myid, mybegin, myend, shift, i)
		for (myid = 0; myid < nthreads; myid ++)
		{
			FASP_GET_START_END(myid, nthreads, A->row, mybegin, myend);
			shift = 0;
			if (myid > 0) {
				shift = indexs[myid-1];
			}
			for (i=mybegin;i<myend;++i) {
				if(vec[i]==1)
				{
					CoarseIndex[i] += shift;
				}
			}
		}
		fasp_mem_free(indexs);
	}
	else {
		index=0;
		for(i=0;i<A->row;++i) {
			if(vec[i]==1)
			{
				CoarseIndex[i]=index;
				index++;
			}
		}
	}
#endif
	if (P.IA[P.row] > openmp_holds) {
		stride_i = P.IA[P.row]/nthreads;
#pragma omp parallel private(myid,mybegin,myend,i,j) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1) myend = mybegin+stride_i;
			else myend = P.IA[P.row];
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
