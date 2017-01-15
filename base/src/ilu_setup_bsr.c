/*! \file ilu_setup_bsr.c
 *
 *  \brief Setup incomplete LU decomposition for dBSRmat matrices
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/* The following functions are defined in ilu.f */
#ifdef __cplusplus 
extern "C" {void symbfactor_(const INT *n,INT *colind,INT *rwptr,
                             const INT *levfill,const INT *nzmax,
                             INT *nzlu,INT *ijlu,INT *uptr,INT *ierr);}
#else
extern void symbfactor_(const INT *n,INT *colind,INT *rwptr,const INT *levfill,
                        const INT *nzmax,INT *nzlu,INT *ijlu,INT *uptr,INT *ierr);
#endif

static INT numfac_bsr(dBSRmat *A, REAL *luval, INT *jlu, INT *uptr);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_ilu_dbsr_setup (dBSRmat *A, ILU_data *iludata, ILU_param *iluparam)
 *
 * \brief Get ILU decoposition of a BSR matrix A
 *
 * \param A         Pointer to dBSRmat matrix
 * \param iludata   Pointer to ILU_data
 * \param iluparam  Pointer to ILU_param
 * 
 * \return          FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   11/08/2010
 *
 * \note Works for general nb (Xiaozhe)
 * \note Change the size of work space by Zheng Li 04/26/2015.
 */
SHORT fasp_ilu_dbsr_setup (dBSRmat    *A,
                           ILU_data   *iludata,
                           ILU_param  *iluparam)
{
        
    const SHORT  prtlvl = iluparam->print_level;
    const INT    n = A->COL, nnz = A->NNZ, nb = A->nb, nb2 = nb*nb;
    
    // local variables
    INT lfil=iluparam->ILU_lfil;
    INT ierr, iwk, nzlu, nwork, *ijlu, *uptr;
    
    REAL    setup_start, setup_end, setup_duration;
    SHORT   status = FASP_SUCCESS;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n", A->ROW, n, nnz);
#endif
    
    fasp_gettime(&setup_start);

    // Expected amount of memory for ILU needed and allocate memory 
    iwk = (lfil+2)*nnz;
    
    // setup preconditioner
    iludata->row = iludata->col=n;
    iludata->nb  = nb;
    iludata->ilevL = iludata->jlevL=NULL;
    iludata->ilevU = iludata->jlevU=NULL;
    
    ijlu = (INT*)fasp_mem_calloc(iwk,sizeof(INT));
    uptr = (INT*)fasp_mem_calloc(A->ROW,sizeof(INT));
    
#if DEBUG_MODE > 1
    printf("### DEBUG: symbolic factorization ... \n ");
#endif
    
    // ILU decomposition    
    // (1) symbolic factoration
    // previous Fortran version
    // symbfactor_(&A->ROW,A->JA,A->IA,&lfil,&iwk,&nzlu,ijlu,uptr,&ierr);
    fasp_symbfactor(A->ROW,A->JA,A->IA,lfil,iwk,&nzlu,ijlu,uptr,&ierr);
    
    iludata->luval = (REAL*)fasp_mem_calloc(nzlu*nb2,sizeof(REAL));
    
#if DEBUG_MODE > 1
    printf("### DEBUG: numerical factorization ... \n ");
#endif
    
    // (2) numerical factoration 
    numfac_bsr(A, iludata->luval, ijlu, uptr);
    
    //nwork = 6*nzlu*nb;
    nwork = 20*A->ROW*A->nb;
    iludata->nzlu  = nzlu;
    iludata->nwork = nwork;
    iludata->ijlu  = (INT*)fasp_mem_calloc(nzlu,sizeof(INT));
    
    memcpy(iludata->ijlu,ijlu,nzlu*sizeof(INT));
    iludata->work = (REAL*)fasp_mem_calloc(nwork, sizeof(REAL));
    // Check: Is the work space too large? --Xiaozhe
    
#if DEBUG_MODE > 1
    printf("### DEBUG: fill-in = %d, nwork = %d\n", lfil, nwork);
    printf("### DEBUG: iwk = %d, nzlu = %d\n",iwk,nzlu);
#endif
    
    if ( ierr != 0 ) {
        printf("### ERROR: ILU setup failed (ierr=%d)!\n", ierr);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if ( iwk < nzlu ) {
        printf("### ERROR: Need more memory for ILU %d!\n", iwk-nzlu);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;    
        printf("BSR ILU(%d) setup costs %f seconds.\n", lfil,setup_duration);    
    }
    
 FINISHED:     
    fasp_mem_free(ijlu);
    fasp_mem_free(uptr);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static INT numfac_bsr (dBSRmat *A, REAL *luval, INT *jlu, INT *uptr)
 * \brief Get numerical ILU decoposition of a BSR matrix A
 *
 * \param A        Pointer to dBSRmat matrix
 * \param luval    Pointer to numerical value of ILU
 * \param jlu      Pointer to the nonzero pattern of ILU
 * \param uptr     Pointer to the diagnal position of ILU
 *
 * \author Shiquan Zhang
 * \date 11/08/2010
 *
 * \note Works for general nb (Xiaozhe)
 */
static INT numfac_bsr (dBSRmat   *A,
                       REAL      *luval,
                       INT       *jlu,
                       INT       *uptr)
{
    INT n=A->ROW,nb=A->nb, nb2=nb*nb, ib, ibstart,ibstart1;
    INT k, indj, inds, indja,jluj, jlus, ijaj;
    REAL  *mult,*mult1;
    INT *colptrs;
    INT status=FASP_SUCCESS;
    
    colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
    mult=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
    mult1=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
    
    /**
     *     colptrs is used to hold the indices of entries in LU of row k.
     *     It is initialized to zero here, and then reset after each row's
     *     work. The first segment of the loop on indj effectively solves
     *     the transposed upper triangular system
     *            U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'
     *     via sparse saxpy operations, throwing away disallowed fill.
     *     When the loop index indj reaches the k-th column (i.e., the diag
     *     entry), then the innermost sparse saxpy operation effectively is
     *     applying the previous updates to the corresponding part of U via
     *     sparse vec*mat, discarding disallowed fill-in entries, i.e.
     *            U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
     */    
    
    //for (k=0;k<n;k++) colptrs[k]=0;
    memset(colptrs, 0, sizeof(INT)*n);

    switch (nb) {
    
    case 1:
    
        for (k = 0; k < n; ++k) {
    
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
    
            colptrs[k] =  k;
    
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
    
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
    
                jluj = jlu[indj];
    
                luval[indj] = luval[indj]*luval[jluj];
                mult[0] = luval[indj];
    
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
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
    
        for (k = 0; k < n; ++k) {
    
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
    
            colptrs[k] =  k;
    
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
    
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
    
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc3(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
    
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
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
    
        for (k = 0; k < n; ++k) {
    
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
    
            colptrs[k] =  k;
    
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
    
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
    
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc5(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
    
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
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
    
        for (k = 0; k < n; ++k) {
    
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
    
            colptrs[k] =  k;
    
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
    
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
    
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc7(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
    
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
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
    
        for (k=0;k<n;k++) {
    
            for (indj = jlu[k];indj<jlu[k+1];++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
    
            colptrs[k] =  k;
    
            for (indja = A->IA[k]; indja < A->IA[k+1];indja++) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
    
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
    
                ibstart=indj*nb2;
                fasp_blas_smat_mul(&(luval[ibstart]),&(luval[jluj*nb2]),mult,nb);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
    
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; inds++) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
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


/**
 * \fn static INT numfac_bsr_mc_omp (dBSRmat *A, REAL *luval, INT *jlu, 
 *                                   INT *uptr, INT ncolors, INT *ic, INT *icmap)
 * \brief Multi-thread ILU decoposition of a BSR matrix A based on graph coloring
 *
 * \param A        Pointer to dBSRmat matrix
 * \param luval    Pointer to numerical value of ILU
 * \param jlu      Pointer to the nonzero pattern of ILU
 * \param uptr     Pointer to the diagnal position of ILU
 * \param ncolors  Number of colors of adjacency graph of A 
 * \param ic       Pointer to number of vertices in each color  
 * \param icmap    Mapping 
 * 
 * \author Zheng Li
 * \date 12/04/2016
 *
 * \note Only works for 1, 2, 3 nb (Zheng)
 */
static INT numfac_bsr_mc_omp (dBSRmat   *A,
                              REAL      *luval,
                              INT       *jlu,
                              INT       *uptr,
                              INT        ncolors,
                              INT       *ic,
                              INT       *icmap)
{
    INT status = FASP_SUCCESS;

#ifdef _OPENMP
    INT n = A->ROW, nb = A->nb, nb2 = nb*nb;
    INT ib, ibstart,ibstart1;
    INT k, i, indj, inds, indja,jluj, jlus, ijaj, tmp;
    REAL  *mult,*mult1;
    INT *colptrs;
    
    /**
     *     colptrs is used to hold the indices of entries in LU of row k.
     *     It is initialized to zero here, and then reset after each row's
     *     work. The first segment of the loop on indj effectively solves
     *     the transposed upper triangular system
     *            U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'
     *     via sparse saxpy operations, throwing away disallowed fill.
     *     When the loop index indj reaches the k-th column (i.e., the diag
     *     entry), then the innermost sparse saxpy operation effectively is
     *     applying the previous updates to the corresponding part of U via
     *     sparse vec*mat, discarding disallowed fill-in entries, i.e.
     *            U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
     */    

    switch (nb) {
    
    case 1:
        for (i = 0; i < ncolors; ++i) {
#pragma omp parallel private(k,indj,ibstart,ib,indja,ijaj,ibstart1,jluj,inds,jlus,colptrs,tmp)
{
        colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
        memset(colptrs, 0, sizeof(INT)*n);
#pragma omp for 
        for (k = ic[i]; k < ic[i+1]; ++k) {
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
            colptrs[k] =  k;
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
                luval[indj] = luval[indj]*luval[jluj];
                tmp = luval[indj];
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0)
                        luval[colptrs[jlus]] = luval[colptrs[jlus]] - tmp*luval[inds];
                }
    
            }
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
            colptrs[k] =  0;
            luval[k] = 1.0/luval[k];
        } 
        fasp_mem_free(colptrs);
        }
 }

        break;
            
    case 2:

        for (i = 0; i < ncolors; ++i) {
#pragma omp parallel private(k,indj,ibstart,ib,indja,ijaj,ibstart1,jluj,inds,jlus,mult,mult1,colptrs)
{
        colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
        memset(colptrs, 0, sizeof(INT)*n);
        mult=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
        mult1=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
#pragma omp for 
        for (k = ic[i]; k < ic[i+1]; ++k) {
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
            colptrs[k] =  k;
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc2(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
                        fasp_blas_smat_mul_nc2(mult,&(luval[inds*nb2]),mult1);
                        ibstart=colptrs[jlus]*nb2;
                        for (ib=0;ib<nb2;++ib) luval[ibstart+ib]-=mult1[ib];
                    }
                }
            }
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
            colptrs[k] =  0;
            fasp_blas_smat_inv_nc2(&(luval[k*nb2]));
        }
       fasp_mem_free(colptrs);
       fasp_mem_free(mult);
       fasp_mem_free(mult1);
}
        }
        break;    
    
    case 3:

        for (i = 0; i < ncolors; ++i) {
#pragma omp parallel private(k,indj,ibstart,ib,indja,ijaj,ibstart1,jluj,inds,jlus,mult,mult1,colptrs)
{
        colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
        memset(colptrs, 0, sizeof(INT)*n);
        mult=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
        mult1=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
#pragma omp for 
        for (k = ic[i]; k < ic[i+1]; ++k) {
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
            colptrs[k] =  k;
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc3(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
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
        fasp_mem_free(colptrs);
        fasp_mem_free(mult);
        fasp_mem_free(mult1);
}
        }
        break;    

        default:
        {
            if (nb > 3) printf("Multi-thread ILU numerical decomposition for %d\
                                components has not been implemented!!!", nb);
            exit(0);
        }
    }
    
#endif

    return status;
}

/**
 * \fn static INT numfac_bsr_levsch_omp (dBSRmat *A, REAL *luval, INT *jlu, 
 *                                       INT *uptr, INT ncolors, INT *ic, INT *icmap)
 * \brief Multi-thread ILU decoposition of a BSR matrix A based on level schedule strategy
 *
 * \param A        Pointer to dBSRmat matrix
 * \param luval    Pointer to numerical value of ILU
 * \param jlu      Pointer to the nonzero pattern of ILU
 * \param uptr     Pointer to the diagnal position of ILU
 * \param ncolors  Number of colors of adjacency graph of A 
 * \param ic       Pointer to number of vertices in each color  
 * \param icmap    Mapping 
 * 
 * \author Zheng Li
 * \date 12/04/2016
 *
 * \note Only works for 1, 2, 3 nb (Zheng)
 */
static INT numfac_bsr_levsch_omp (dBSRmat *A,
                                  REAL *luval,
                                  INT *jlu,
                                  INT *uptr,
                                  INT ncolors,
                                  INT *ic,
                                  INT *icmap)
{
    INT status = FASP_SUCCESS;

#ifdef _OPENMP
    INT n = A->ROW, nb = A->nb, nb2 = nb*nb;
    INT ib, ibstart,ibstart1;
    INT k, i, indj, inds, indja, jluj, jlus, ijaj, tmp, ii;
    REAL *mult, *mult1;
    INT  *colptrs;
    
    /**
     *     colptrs is used to hold the indices of entries in LU of row k.
     *     It is initialized to zero here, and then reset after each row's
     *     work. The first segment of the loop on indj effectively solves
     *     the transposed upper triangular system
     *            U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'
     *     via sparse saxpy operations, throwing away disallowed fill.
     *     When the loop index indj reaches the k-th column (i.e., the diag
     *     entry), then the innermost sparse saxpy operation effectively is
     *     applying the previous updates to the corresponding part of U via
     *     sparse vec*mat, discarding disallowed fill-in entries, i.e.
     *            U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
     */    

    switch (nb) {
    
    case 1:
        for (i = 0; i < ncolors; ++i) {
#pragma omp parallel private(k,indj,ibstart,ib,indja,ijaj,ibstart1,jluj,inds,jlus,colptrs,tmp)
{
        colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
        memset(colptrs, 0, sizeof(INT)*n);
#pragma omp for 
        for (k = ic[i]; k < ic[i+1]; ++k) {
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
            colptrs[k] =  k;
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
                luval[indj] = luval[indj]*luval[jluj];
                tmp = luval[indj];
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0)
                        luval[colptrs[jlus]] = luval[colptrs[jlus]] - tmp*luval[inds];
                }
    
            }
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
            colptrs[k] =  0;
            luval[k] = 1.0/luval[k];
        } 
        fasp_mem_free(colptrs);
        }
 }

        break;
    case 2:

        for (i = 0; i < ncolors; ++i) {
#pragma omp parallel private(k,indj,ibstart,ib,indja,ijaj,ibstart1,jluj,inds,jlus,mult,mult1,colptrs,ii)
{
        colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
        memset(colptrs, 0, sizeof(INT)*n);
        mult=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
        mult1=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
#pragma omp for 
        for (ii = ic[i]; ii < ic[i+1]; ++ii) {
            k = icmap[ii];
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
            colptrs[k] =  k;
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc2(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
                        fasp_blas_smat_mul_nc2(mult,&(luval[inds*nb2]),mult1);
                        ibstart=colptrs[jlus]*nb2;
                        for (ib=0;ib<nb2;++ib) luval[ibstart+ib]-=mult1[ib];
                    }
                }
            }
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) colptrs[jlu[indj]] = 0;
            colptrs[k] =  0;
            fasp_blas_smat_inv_nc2(&(luval[k*nb2]));
        }
        fasp_mem_free(colptrs);
        fasp_mem_free(mult);
        fasp_mem_free(mult1);
}
        }
        break;    
    
    case 3:

        for (i = 0; i < ncolors; ++i) {
#pragma omp parallel private(k,indj,ibstart,ib,indja,ijaj,ibstart1,jluj,inds,jlus,mult,mult1,colptrs,ii)
{
        colptrs=(INT*)fasp_mem_calloc(n,sizeof(INT));
        memset(colptrs, 0, sizeof(INT)*n);
        mult=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
        mult1=(REAL*)fasp_mem_calloc(nb2,sizeof(REAL));
#pragma omp for 
        for (ii = ic[i]; ii < ic[i+1]; ++ii) {
            k = icmap[ii];
            for (indj = jlu[k]; indj < jlu[k+1]; ++indj) {
                colptrs[jlu[indj]] = indj;
                ibstart=indj*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = 0;
            }
            colptrs[k] =  k;
            for (indja = A->IA[k]; indja < A->IA[k+1]; ++indja) {
                ijaj = A->JA[indja];
                ibstart=colptrs[ijaj]*nb2;
                ibstart1=indja*nb2;
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib] = A->val[ibstart1+ib];
            }
            for (indj = jlu[k]; indj < uptr[k]; ++indj) {
                jluj = jlu[indj];
                ibstart=indj*nb2;
                fasp_blas_smat_mul_nc3(&(luval[ibstart]),&(luval[jluj*nb2]),mult);
                for (ib=0;ib<nb2;++ib) luval[ibstart+ib]=mult[ib];
                for (inds = uptr[jluj]; inds < jlu[jluj+1]; ++inds) {
                    jlus = jlu[inds];
                    if (colptrs[jlus] != 0) {
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
        fasp_mem_free(colptrs);
        fasp_mem_free(mult);
        fasp_mem_free(mult1);
}
        }
        break;    

        default:
        {
            if (nb > 3) printf("Multi-thread ILU numerical decomposition for %d \
                                components has not been implemented!!!", nb);
                exit(0);
            break;
        }
    }
    
#endif

    return status;
}

/**
 * \fn SHORT fasp_ilu_dbsr_setup_levsch_omp (dBSRmat *A, ILU_data *iludata, ILU_param *iluparam)
 *
 * \brief Get ILU decoposition of a BSR matrix A based on level schedule strategy
 *
 * \param A         Pointer to dBSRmat matrix
 * \param iludata   Pointer to ILU_data
 * \param iluparam  Pointer to ILU_param
 * 
 * \return          FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Zheng Li
 * \date   12/04/2016
 *
 * \note Only works for 1, 2, 3 nb (Zheng)
 */
SHORT fasp_ilu_dbsr_setup_levsch_omp (dBSRmat    *A,
                                      ILU_data   *iludata,
                                      ILU_param  *iluparam)
{
    const SHORT  prtlvl = iluparam->print_level;
    const INT    n = A->COL, nnz = A->NNZ, nb = A->nb, nb2 = nb*nb;
    
    // local variables
    INT lfil=iluparam->ILU_lfil;
    INT ierr, iwk, nzlu, nwork, *ijlu, *uptr;
    
    REAL    setup_start, setup_end, setup_duration;
    REAL    symbolic_start, symbolic_end, numfac_start, numfac_end;
    
    SHORT   status = FASP_SUCCESS;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",A->ROW,n,nnz);
#endif
    
    fasp_gettime(&setup_start);

    // Expected amount of memory for ILU needed and allocate memory 
    iwk = (lfil+2)*nnz;
    
    // setup preconditioner
    iludata->row = iludata->col=n;
    iludata->nb  = nb;
    
    ijlu = (INT*)fasp_mem_calloc(iwk,sizeof(INT));
    uptr = (INT*)fasp_mem_calloc(A->ROW,sizeof(INT));

        
#if DEBUG_MODE > 1
    printf("### DEBUG: symbolic factorization ... \n ");
#endif

    fasp_gettime(&symbolic_start);
    
    // ILU decomposition    
    // (1) symbolic factoration
    fasp_symbfactor(A->ROW,A->JA,A->IA,lfil,iwk,&nzlu,ijlu,uptr,&ierr);
    
    fasp_gettime(&symbolic_end);

    printf("symbolic time=%f\n", symbolic_end-symbolic_start);
    
    nwork = 5*A->ROW*A->nb; 
    iludata->nzlu  = nzlu;
    iludata->nwork = nwork;
    iludata->ijlu  = (INT*)fasp_mem_calloc(nzlu,sizeof(INT));
    iludata->luval = (REAL*)fasp_mem_calloc(nzlu*nb2,sizeof(REAL));
    iludata->work = (REAL*)fasp_mem_calloc(nwork, sizeof(REAL));
    memcpy(iludata->ijlu,ijlu,nzlu*sizeof(INT));
    fasp_array_set(nzlu*nb2, iludata->luval, 0.0);
    iludata->uptr = NULL,iludata->ic = NULL, iludata->icmap = NULL;

    fasp_topological_sorting_ilu(iludata);

#if DEBUG_MODE > 1
    printf("### DEBUG: numerical factorization ... \n ");
#endif
    
    fasp_gettime(&numfac_start);
    
    // (2) numerical factoration 
    numfac_bsr_levsch_omp(A, iludata->luval, ijlu, uptr, iludata->nlevL, iludata->ilevL, iludata->jlevL);
        
    fasp_gettime(&numfac_end);

    printf("numfac time =%f\n", numfac_end-numfac_start);
    
#if DEBUG_MODE > 1
    printf("### DEBUG: fill-in = %d, nwork = %d\n", lfil, nwork);
    printf("### DEBUG: iwk = %d, nzlu = %d\n",iwk,nzlu);
#endif
    
    if ( ierr != 0 ) {
        printf("### ERROR: ILU setup failed (ierr=%d)!\n", ierr);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if ( iwk < nzlu ) {
        printf("### ERROR: Need more memory for ILU %d!\n", iwk-nzlu);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;    
        printf("BSR ILU(%d) setup costs %f seconds.\n", lfil,setup_duration);    
    }
    
 FINISHED:     
    fasp_mem_free(ijlu);
    fasp_mem_free(uptr);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/**
 * \fn SHORT fasp_ilu_dbsr_setup_omp (dBSRmat *A, ILU_data *iludata, ILU_param *iluparam)
 *
 * \brief Multi-threads parallel ILU decoposition of a BSR matrix A based on graph coloring
 *
 * \param A         Pointer to dBSRmat matrix
 * \param iludata   Pointer to ILU_data
 * \param iluparam  Pointer to ILU_param
 * 
 * \return          FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Zheng Li
 * \date   12/04/2016
 *
 * \note Only works for 1, 2, 3 nb (Zheng)
 */
SHORT fasp_ilu_dbsr_setup_omp (dBSRmat    *A,
                               ILU_data   *iludata,
                               ILU_param  *iluparam)
{
    
    const SHORT  prtlvl = iluparam->print_level;
    const INT    n = A->COL, nnz = A->NNZ, nb = A->nb, nb2 = nb*nb;
    
    // local variables
    INT lfil=iluparam->ILU_lfil;
    INT ierr, iwk, nzlu, nwork, *ijlu, *uptr;
    
    REAL    setup_start, setup_end, setup_duration;
    
    SHORT   status = FASP_SUCCESS;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",A->ROW,n,nnz);
#endif
    
    fasp_gettime(&setup_start);

    // Expected amount of memory for ILU needed and allocate memory 
    iwk = (lfil+2)*nnz;
    
    // setup preconditioner
    iludata->row = iludata->col=n;
    iludata->nb  = nb;
    
    ijlu = (INT*)fasp_mem_calloc(iwk,sizeof(INT));
    uptr = (INT*)fasp_mem_calloc(A->ROW,sizeof(INT));
        
#if DEBUG_MODE > 1
    printf("### DEBUG: symbolic factorization ... \n ");
#endif
    
    // ILU decomposition    
    // (1) symbolic factoration
    fasp_symbfactor(A->ROW,A->JA,A->IA,lfil,iwk,&nzlu,ijlu,uptr,&ierr);
    
    nwork = 5*A->ROW*A->nb; 
    iludata->nzlu  = nzlu;
    iludata->nwork = nwork;
    iludata->ijlu  = (INT*)fasp_mem_calloc(nzlu,sizeof(INT));
    iludata->luval = (REAL*)fasp_mem_calloc(nzlu*nb2,sizeof(REAL));
    iludata->work = (REAL*)fasp_mem_calloc(nwork, sizeof(REAL));
    memcpy(iludata->ijlu,ijlu,nzlu*sizeof(INT));
    fasp_array_set(nzlu*nb2, iludata->luval, 0.0);


#if DEBUG_MODE > 1
    printf("### DEBUG: numerical factorization ... \n ");
#endif
    
    // (2) numerical factoration 
    numfac_bsr_mc_omp(A, iludata->luval, ijlu, uptr, iludata->nlevL, iludata->ilevL, iludata->jlevL);

    //printf("numfac time =%f\n", numfac_end-numfac_start);

#if DEBUG_MODE > 1
    printf("### DEBUG: fill-in = %d, nwork = %d\n", lfil, nwork);
    printf("### DEBUG: iwk = %d, nzlu = %d\n",iwk,nzlu);
#endif
    
    if ( ierr != 0 ) {
        printf("### ERROR: ILU setup failed (ierr=%d)!\n", ierr);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if ( iwk < nzlu ) {
        printf("### ERROR: Need more memory for ILU %d!\n", iwk-nzlu);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;    
        printf("BSR ILU(%d) setup costs %f seconds.\n", lfil,setup_duration);    
    }
    
 FINISHED:     
    fasp_mem_free(ijlu);
    fasp_mem_free(uptr);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/**
 * \fn SHORT fasp_ilu_dbsr_setup_mc_omp (dBSRmat *A, dCSRmat *Ap, ILU_data *iludata, ILU_param *iluparam)
 *
 * \brief Multi-threads parallel ILU decoposition of a BSR matrix A based on graph coloring
 *
 * \param A         Pointer to dBSRmat matrix
 * \param Ap        Pointer to dCSRmat matrix and provide sparsity pattern 
 * \param iludata   Pointer to ILU_data
 * \param iluparam  Pointer to ILU_param
 * 
 * \return          FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Zheng Li
 * \date   12/04/2016
 *
 * \note Only works for 1, 2, 3 nb (Zheng)
 */
SHORT fasp_ilu_dbsr_setup_mc_omp (dBSRmat    *A,
                                  dCSRmat    *Ap,
                                  ILU_data   *iludata,
                                  ILU_param  *iluparam)
{
     INT status;
     AMG_data *mgl=fasp_amg_data_create(1);
     dCSRmat pp, Ap1;
     dBSRmat A_LU;

    if (iluparam->ILU_lfil==0) {
        mgl[0].A = fasp_dcsr_sympart(Ap);  //for ILU0
    }
    else if (iluparam->ILU_lfil==1) {  // for ILU1
        Ap1 = fasp_dcsr_create(Ap->row,Ap->col, Ap->nnz);
        fasp_dcsr_cp(Ap, &Ap1);
        fasp_blas_dcsr_mxm (Ap,&Ap1,&pp);
        mgl[0].A = fasp_dcsr_sympart(&pp);
        fasp_dcsr_free(&Ap1);
        fasp_dcsr_free(&pp);
    }

    mgl->num_levels = 20;

    fasp_multicolors_independent_set(mgl, 1);

    A_LU = fasp_dbsr_perm(A, mgl[0].icmap);
    
    // hold color info with nlevl, ilevL and jlevL.
    iludata->nlevL = mgl[0].colors;
    iludata->ilevL = mgl[0].ic;
    iludata->jlevL = mgl[0].icmap;
    iludata->nlevU = 0;
    iludata->ilevU = NULL;
    iludata->jlevU = NULL;

    status = fasp_ilu_dbsr_setup_omp(&A_LU,iludata,iluparam);

    fasp_dcsr_free(&mgl[0].A);
    fasp_dbsr_free(&A_LU);

    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
