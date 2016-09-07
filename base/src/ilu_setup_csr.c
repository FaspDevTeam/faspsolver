/*! \file ilu_setup_csr.c
 *
 *  \brief Setup incomplete LU decomposition for dCSRmat matrices
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#define ILU_C_VERSION TRUE /**< Use the C version of ILU functions */

/* declarations for ilu.for */
#ifdef __cplusplus 
extern "C" {void iluk_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,REAL *alu,
                       INT *jlu,INT *iwk,INT *ierr,INT *nzlu);}
extern "C" {void ilut_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,
                       const REAL *droptol,REAL *alu,INT *jlu,INT *iwk,
                       INT *ierr,INT *nz);}
extern "C" {void ilutp_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,
                        const REAL *droptol,const REAL *permtol,const INT *mbloc,
                        REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);}
#else
extern void iluk_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,REAL *alu,
                  INT *jlu,INT *iwk,INT *ierr,INT *nzlu);
extern void ilut_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,
                  REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);
extern void ilutp_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,
                   const REAL *permtol,const INT *mbloc,REAL *alu,INT *jlu,INT *iwk,
                   INT *ierr,INT *nz);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_ilu_dcsr_setup (dCSRmat *A, ILU_data *iludata, ILU_param *iluparam)
 *
 * \brief Get ILU decomposition of a CSR matrix A
 *
 * \param A         Pointer to dCSRmat matrix
 * \param iludata   Pointer to ILU_data
 * \param iluparam  Pointer to ILU_param
 *
 * \return          FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Shiquan Zhang Xiaozhe Hu
 * \date   12/27/2009
 */
SHORT fasp_ilu_dcsr_setup (dCSRmat *A, 
                           ILU_data *iludata, 
                           ILU_param *iluparam)
{
#if FASP_USE_ILU || ILU_C_VERSION
    const INT   type=iluparam->ILU_type, print_level=iluparam->print_level;
    const INT   n=A->col, nnz=A->nnz, mbloc=n;
    const REAL  ILU_droptol=iluparam->ILU_droptol;
    const REAL  permtol=iluparam->ILU_permtol;
    
    // local variable
    INT    lfil=iluparam->ILU_lfil, lfilt=iluparam->ILU_lfil;
    INT    ierr, iwk, nzlu, nwork, *ijlu;
    REAL  *luval;
    
    REAL   setup_start, setup_end, setup_duration;
    SHORT  status = FASP_SUCCESS;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n", A->row, n, nnz);
#endif
    
    fasp_gettime(&setup_start);
    
    // Expected amount of memory for ILU needed and allocate memory 
    switch (type) {
        case ILUt:
            iwk=10*nnz;     // iwk is the maxim possible nnz for ILU
            lfilt=floor(n*0.5)+1;
            break;
        case ILUtp:
            iwk=10*nnz;     // iwk is the maxim possible nnz for ILU
            lfilt=floor(n*0.5)+1;
            break;
        default: // ILUk
            if (lfil == 0) iwk=nnz+500;
            else iwk=(lfil+5)*nnz;
            break;
    } 
    
    nwork  = 4*n;
    
#if DEBUG_MODE > 1
    printf("### DEBUG: fill-in = %d, iwk = %d, nwork = %d\n", lfil, iwk, nwork);
#endif
    
    // setup ILU preconditioner
    iludata->row=iludata->col=n;    
    fasp_ilu_data_alloc(iwk, nwork, iludata);
    
#if CHMEM_MODE
    printf("### DEBUG: memory usage after %s: \n", __FUNCTION__);
    fasp_mem_usage();
#endif
    
    // ILU decomposition
    ijlu=iludata->ijlu;
    luval=iludata->luval;
    
    switch (type) {
        case ILUt:
#if DEBUG_MODE > 0
    printf("### DEBUG: %s (ILUt) ...... [Start]\n", __FUNCTION__);
#endif

#if	ILU_C_VERSION
            fasp_ilut(n,A->val,A->JA,A->IA,lfilt,ILU_droptol,luval,ijlu,iwk,&ierr,&nzlu);
#else // Fortran version
            ilut_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,luval,ijlu,&iwk,&ierr,&nzlu);
#endif
            break;
            
        case ILUtp:
#if DEBUG_MODE > 0
    printf("### DEBUG: %s (ILUp) ...... [Start]\n", __FUNCTION__);
#endif
            
#if	ILU_C_VERSION
            fasp_ilutp(n,A->val,A->JA,A->IA, lfilt, ILU_droptol, permtol,
                       mbloc,luval,ijlu,iwk,&ierr,&nzlu);
#else // Fortran version
            ilutp_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,&permtol,
                   &mbloc,luval,ijlu,&iwk,&ierr,&nzlu);
#endif
	 
            break;
            
        default: // ILUk
#if DEBUG_MODE > 0
    printf("### DEBUG: %s (ILUk) ...... [Start]\n", __FUNCTION__);
#endif
            
#if	ILU_C_VERSION
            fasp_iluk(n,A->val,A->JA,A->IA,lfil,luval,ijlu,iwk,&ierr,&nzlu);
#else // Fortran version
            iluk_(&n,A->val,A->JA,A->IA,&lfil,luval,ijlu,&iwk,&ierr,&nzlu);
#endif
            break;
    } 
    
    fasp_dcsr_shift(A, -1);
    
#if CHMEM_MODE
    printf("### DEBUG: memory usage after ILU setup: \n");
    fasp_mem_usage();
#endif
    
    iludata->nzlu=nzlu;
    iludata->nwork=nwork;
    
#if DEBUG_MODE > 1
    printf("### DEBUG: iwk = %d, nzlu = %d\n",iwk,nzlu);
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
        
        switch (type) {
            case ILUt:
                printf("ILUt setup costs %f seconds.\n", setup_duration);    
                break;
            case ILUtp:
                printf("ILUtp setup costs %f seconds.\n", setup_duration);    
                break;
            default: // ILUk
                printf("ILUk setup costs %f seconds.\n", setup_duration);    
                break;
        }     
    }
    
FINISHED:     
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
    
#else // WITH_ILU
    
    printf("### ERROR: ILU is not enabled!\n");
    fasp_chkerr(ERROR_MISC, __FUNCTION__);
    
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
