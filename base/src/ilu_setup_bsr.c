/*! \file ilu_setup_bsr.c
 *
 *  \brief Setup Incomplete LU decomposition for dBSRmat matrices
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/* The following functions are defined in ilu.for */
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
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   11/08/2010
 *
 * \note Works for general nb (Xiaozhe)
 */
SHORT fasp_ilu_dbsr_setup (dBSRmat *A, 
                           ILU_data *iludata, 
                           ILU_param *iluparam)
{
    
    const SHORT  print_level=iluparam->print_level;
    const INT    n=A->COL, nnz=A->NNZ, nb=A->nb, nb2=nb*nb;
    
    // local variables
    INT lfil=iluparam->ILU_lfil;
    INT ierr, iwk, nzlu, nwork, *ijlu, *uptr;
    
    REAL    setup_start, setup_end, setup_duration;
    SHORT   status = FASP_SUCCESS;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_ilu_dbsr_setup ...... [Start]\n");
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",A->ROW,n,nnz);
#endif
    
    fasp_gettime(&setup_start);

    // Expected amount of memory for ILU needed and allocate memory 
    iwk=(lfil+2)*nnz;
    
    // setup preconditioner
    iludata->row=iludata->col=n;    
    iludata->nb=nb;
    
    ijlu=(INT*)fasp_mem_calloc(iwk,sizeof(INT));
    uptr=(INT*)fasp_mem_calloc(A->ROW,sizeof(INT));
        
#if DEBUG_MODE
    printf("### DEBUG: symbolic factorization ... \n ");
#endif
    
    // ILU decomposition    
    // (1) symbolic factoration
    symbfactor_(&A->ROW,A->JA,A->IA,&lfil,&iwk,&nzlu,ijlu,uptr,&ierr);
    
    iludata->luval=(REAL*)fasp_mem_calloc(nzlu*nb2,sizeof(REAL));
    
#if DEBUG_MODE
    printf("### DEBUG: numerical factorization ... \n ");
#endif
    
    // (2) numerical factoration 
    numfac_bsr(A, iludata->luval, ijlu, uptr);
        
    nwork=6*nzlu*nb;
    iludata->nzlu=nzlu;
    iludata->nwork=nwork;
    iludata->ijlu=(INT*)fasp_mem_calloc(nzlu,sizeof(INT));
    
    memcpy(iludata->ijlu,ijlu,nzlu*sizeof(INT));
    iludata->work=(REAL*)fasp_mem_calloc(nwork, sizeof(REAL));  // Xiaozhe: Is the work space too large?

#if DEBUG_MODE
    printf("### DEBUG: fill-in=%d, nwork=%d\n", lfil, nwork);
    printf("### DEBUG: iwk=%d, nzlu=%d\n",iwk,nzlu);
#endif
    
    if (ierr!=0) {
        printf("### ERROR: ILU setup failed (ierr=%d)!\n", ierr);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if (iwk<nzlu) {
        printf("### ERROR: Need more memory for ILU %d!\n", iwk-nzlu);
        status = ERROR_SOLVER_ILUSETUP;
        goto FINISHED;
    }
    
    if (print_level>PRINT_NONE) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;    
        printf("BSR ILU(%d) setup costs %f seconds.\n", lfil,setup_duration);    
    }
    
 FINISHED:     
    fasp_mem_free(ijlu);
    fasp_mem_free(uptr);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_ilu_dbsr_setup ...... [Finish]\n");
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
static INT numfac_bsr (dBSRmat *A,
                       REAL *luval,
                       INT *jlu,
                       INT *uptr)
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
