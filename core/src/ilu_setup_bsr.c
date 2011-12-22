/*! \file ilu_setup_bsr.c
 *  \brief Setup incomplete LU decomposition.
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/* ilu.for */
#ifdef __cplusplus 
extern "C" {void symbfactor_(const int *n,int *colind,int *rwptr,const int *levfill,const int *nzmax,int *nzlu,int *ijlu,int *uptr,int *ierr);}
#else
extern void symbfactor_(const int *n,int *colind,int *rwptr,const int *levfill,const int *nzmax,int *nzlu,int *ijlu,int *uptr,int *ierr);
#endif

static int numfac_bsr(dBSRmat *A, double *luval, int *jlu, int *uptr);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_ilu_dbsr_setup(dBSRmat *A, ILU_data *iludata, ILU_param *param)
 * \brief Get ILU decoposition of a BSR matrix A
 *
 * \param *A pointer to bSR matrir of double type
 * \param *iludata pointer to ILU_data
 * \param *param pointer to ILU parameters
 *
 * \note: works for general nb (Xiaozhe)
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 11/08/2010
 */
int fasp_ilu_dbsr_setup (dBSRmat *A, 
                         ILU_data *iludata, 
												 ILU_param *param)
{
	const int print_level=param->print_level;
	const int m=A->ROW, n=A->COL, nnz=A->NNZ, nb=A->nb, nb2=nb*nb;
	
	int lfil=param->ILU_lfil;
	int ierr, iwk, nzlu, nwork, *ijlu, *uptr;
	
	clock_t setup_start, setup_end;
	double setup_duration;
	int status = SUCCESS;
	
#if DEBUG_MODE
	printf("fasp_ilu_dbsr_setup ...... [Start]\n");
#endif
	
	if (param->print_level>8) printf("fasp_ilu_dbsr_setup: m=%d, n=%d, nnz=%d\n",m,n,nnz);
	
	setup_start=clock();
	
	// Expected amount of memory for ILU needed and allocate memory 
	iwk=(lfil+2)*nnz;
	
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_dbsr_setup: fill-in=%d, iwk=%d, nwork=%d\n", lfil, iwk, nwork);
#endif
	
	// setup preconditioner
	iludata->row=iludata->col=n;	
	iludata->nb=nb;
	
	ijlu=(int*)fasp_mem_calloc(iwk,sizeof(int));
	uptr=(int*)fasp_mem_calloc(A->ROW,sizeof(int));
	
#if CHMEM_MODE
	printf("ilu_setup: memory usage after before_ILU_data: \n");
	fasp_mem_usage();
#endif
		
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_dbsr_setup: symbolic factorization ... ");
#endif

	// ILU decomposition	
	// (1) symbolic factoration
	symbfactor_(&A->ROW,A->JA,A->IA,&lfil,&iwk,&nzlu,ijlu,uptr,&ierr);	
	
#if DEBUG_MODE
	printf("Done!!!\n");
#endif
	
	iludata->luval=(double*)fasp_mem_calloc(nzlu*nb2,sizeof(double));
	
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_dbsr_setup: numerical factorization ... ");
#endif
	
	// (2) numerical factoration 
	numfac_bsr(A, iludata->luval, ijlu, uptr);
	
#if DEBUG_MODE
	printf("Done!!!\n");
#endif
	
#if CHMEM_MODE
	printf("fasp_ilu_dbsr_setup: memory usage after ILU setup: \n");
	fasp_mem_usage();
#endif
	
	nwork=6*nzlu*nb;
	iludata->nzlu=nzlu;
	iludata->nwork=nwork;
	iludata->ijlu=(int*)fasp_mem_calloc(nzlu,sizeof(int));
	
	memcpy(iludata->ijlu,ijlu,nzlu*sizeof(int));
	iludata->work=(double*)fasp_mem_calloc(nwork, sizeof(double));  // Xiaozhe: Is the work space too large?
	
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_dbsr_setup: iwk=%d, nzlu=%d\n",iwk,nzlu);
#endif	
	
	if (ierr!=0) {
		printf("Error: ILU setup failed (ierr=%d)!\n", ierr);
		status = ERROR_SOLVER_ILUSETUP;
		goto FINISHED;
	}
	
	if (iwk<nzlu) {
		printf("Error: Need more memory for ILU %d!\n", iwk-nzlu);
		status = ERROR_SOLVER_ILUSETUP;
		goto FINISHED;
	}
	
	if (print_level>0) {
		setup_end=clock();
		setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
		
		printf("BSR ILU(%d) setup costs %f seconds.\n", lfil,setup_duration);	
	}
	
FINISHED: 	
	fasp_mem_free(ijlu);
	fasp_mem_free(uptr);
	
#if DEBUG_MODE
	printf("fasp_ilu_dbsr_setup ...... [Finish]\n");
#endif
	
	return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static int numfac_bsr(dBSRmat *A, double *luval, int *jlu, int *uptr)
 * \brief Get numerical ILU decoposition of a BSR matrix A
 *
 * \param *A      pointer to BSR matrir of double type
 * \param *luval  pointer to numerical value of ILU
 * \param *jlu    pointer to the nonzero pattern of ILU
 * \param *uptr   pointer to the diagnal position of ILU
 *
 * \note: works for general nb (Xiaozhe)
 * \author Shiquan Zhang
 * \date 11/08/2010
 */
static int numfac_bsr(dBSRmat *A, double *luval, int *jlu, int *uptr)
{
	int n=A->ROW,nb=A->nb, nb2=nb*nb, ib, ibstart,ibstart1; 
	int k, indj, inds, indja,jluj, jlus, ijaj;
	double  *mult,*mult1; 
	int *colptrs;
	int status=SUCCESS;
	
	colptrs=(int*)fasp_mem_calloc(n,sizeof(int));
	mult=(double*)fasp_mem_calloc(nb2,sizeof(double));
	mult1=(double*)fasp_mem_calloc(nb2,sizeof(double));
	
	/**  
	 *     colptrs is used to hold the indices of entries in LU of row k.  
	 *     It is initialized to zero here, and then reset after each row's work.
	 *     The first segment of the loop on indj effectively solves
	 *     the transposed upper triangular system
	 *     U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'
	 *     via sparse saxpy operations, throwing away disallowed fill.
	 *     When the loop index indj reaches the k-th column (i.e., the
	 *     diagonal entry), then the innermost sparse saxpy operation 
	 *     effectively is applying the previous updates to the corresponding 
	 *     part of U via sparse vector*matrix, discarding disallowed fill-in
	 *     entries.  That operation is 
	 *     U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
	 */    
	
	for (k=0;k<n;k++) colptrs[k]=0;
	
	switch (nb) {
			
		case 1:
			
			for (k = 0; k < n; ++k)
			{
				
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj){
					colptrs[jlu[indj]] = indj;
					ibstart=indj*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
				}
				
				colptrs[k] =  k;
				
				for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja){
					ijaj = A->JA[indja];
					ibstart=colptrs[ijaj]*nb2;
					ibstart1=indja*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
				}
				
				for (indj = jlu[k]; indj < uptr[k]; ++indj)
				{
					
					jluj = jlu[indj];
					
					luval[indj] = luval[indj]*luval[jluj];
					mult[0] = luval[indj];
					
					for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds)
					{
						jlus = jlu[inds];
						if (colptrs[jlus] != 0)
							luval[colptrs[jlus]] = luval[colptrs[jlus]] - mult[0]*luval[inds];
					}
					
				}
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
				
				colptrs[k] =  0;
				luval[k] = 1.0/luval[k];
			} 
			break;
			
		case 3:
			
			for (k = 0; k < n; ++k)
			{
				
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj){
					colptrs[jlu[indj]] = indj;
					ibstart=indj*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
				}
				
				colptrs[k] =  k;
				
				for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja){
					ijaj = A->JA[indja];
					ibstart=colptrs[ijaj]*nb2;
					ibstart1=indja*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
				}
				
				for (indj = jlu[k]; indj < uptr[k]; ++indj)
				{
					jluj = jlu[indj];
					
					ibstart=indj*nb2;
					fasp_blas_smat_mul_nc3(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
					
					for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds)
					{
						jlus = jlu[inds];
						if (colptrs[jlus] != 0)
						{
							fasp_blas_smat_mul_nc3(mult,&(luval[inds*nb2]),mult1);
							ibstart=colptrs[jlus]*nb2;
							for (ib=0;ib<nb2;++ib) luval[ibstart+ib]-=mult1[ib];
						}
					}
					
				}
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
				
				colptrs[k] =  0;
				
				fasp_blas_smat_inv_nc3(&(luval[k*nb2]));
				
			}
			break;
			
		case 5:
			
			for (k = 0; k < n; ++k)
			{
				
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj){
					colptrs[jlu[indj]] = indj;
					ibstart=indj*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
				}
				
				colptrs[k] =  k;
				
				for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja){
					ijaj = A->JA[indja];
					ibstart=colptrs[ijaj]*nb2;
					ibstart1=indja*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
				}
				
				for (indj = jlu[k]; indj < uptr[k]; ++indj)
				{
					jluj = jlu[indj];
					
					ibstart=indj*nb2;
					fasp_blas_smat_mul_nc5(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
					
					for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds)
					{
						jlus = jlu[inds];
						if (colptrs[jlus] != 0)
						{
							fasp_blas_smat_mul_nc5(mult,&(luval[inds*nb2]),mult1);
							ibstart=colptrs[jlus]*nb2;
							for (ib=0;ib<nb2;++ib) luval[ibstart+ib]-=mult1[ib];
						}
					}
					
				}
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
				
				colptrs[k] =  0;
				
				fasp_blas_smat_inv_nc5(&(luval[k*nb2]));
			}
			break;
			
		case 7:
			
			for (k = 0; k < n; ++k)
			{
				
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj){
					colptrs[jlu[indj]] = indj;
					ibstart=indj*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
				}
				
				colptrs[k] =  k;
				
				for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja){
					ijaj = A->JA[indja];
					ibstart=colptrs[ijaj]*nb2;
					ibstart1=indja*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
				}
				
				for (indj = jlu[k]; indj < uptr[k]; ++indj)
				{
					jluj = jlu[indj];
					
					ibstart=indj*nb2;
					fasp_blas_smat_mul_nc7(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
					
					for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds)
					{
						jlus = jlu[inds];
						if (colptrs[jlus] != 0)
						{
							fasp_blas_smat_mul_nc7(mult,&(luval[inds*nb2]),mult1);
							ibstart=colptrs[jlus]*nb2;
							for (ib=0;ib<nb2;++ib) luval[ibstart+ib]-=mult1[ib];
						}
					}
					
				}
				for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
				
				colptrs[k] =  0;
				
				fasp_blas_smat_inv(&(luval[k*nb2]),nb);
			}
			break;
			
		default:
			
			for (k=0;k<n;k++)
			{
				
				for (indj = jlu[k];indj<jlu[k+1];++indj){
					colptrs[jlu[indj]] = indj;
					ibstart=indj*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
				}
				
				colptrs[k] =  k;
				
				for (indja = A->IA[k]; indja < A->IA[k+1];indja++){
					ijaj = A->JA[indja];
					ibstart=colptrs[ijaj]*nb2;
					ibstart1=indja*nb2;
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
				}
				
				for (indj = jlu[k]; indj < uptr[k]; ++indj)
				{
					jluj = jlu[indj];
					
					ibstart=indj*nb2;
					fasp_blas_smat_mul(&(luval[ibstart]),&(luval[jluj*nb2]),mult,nb);
					for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
					
					for (inds = uptr[jluj]; inds < jlu[jluj+1]; inds++)
					{
						jlus = jlu[inds];
						if (colptrs[jlus] != 0)
						{
							fasp_blas_smat_mul(mult,&(luval[inds*nb2]),mult1,nb);
							ibstart=colptrs[jlus]*nb2;
							for (ib=0;ib<nb2;++ib) luval[ibstart+ib]-=mult1[ib];
						}
					}
					
				}
				
				for (indj = jlu[k];indj<jlu[k+1];++indj)
					colptrs[jlu[indj]] = 0;
				
				colptrs[k] =  0;
				
				fasp_blas_smat_inv(&(luval[k*nb2]),nb);
				
			}
	}
	
	fasp_mem_free(colptrs);
	fasp_mem_free(mult);
	fasp_mem_free(mult1);
	
	return status;
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
