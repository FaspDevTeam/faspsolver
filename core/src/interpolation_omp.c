/*
 *  interpolation_omp.c
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
