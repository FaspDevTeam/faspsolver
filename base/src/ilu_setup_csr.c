/*! \file ilu_setup_csr.c
 *  \brief Interface between FASP and ILU packages: setup ILU.
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/* ilu.for */
#ifdef __cplusplus 
extern "C" {void iluk_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nzlu);}
extern "C" {void ilut_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);}
extern "C" {void ilutp_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,const REAL *permtol,const INT *mbloc,REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);}
#else
extern void iluk_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nzlu);
extern void ilut_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);
extern void ilutp_(const INT *n,REAL *a,INT *ja,INT *ia,INT *lfil,const REAL *droptol,const REAL *permtol,const INT *mbloc,REAL *alu,INT *jlu,INT *iwk,INT *ierr,INT *nz);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_ilu_dcsr_setup (dCSRmat *A, ILU_data *iludata, ILU_param *iluparam)
 *
 * \brief Get ILU decoposition of a CSR matrix A
 *
 * \param A         Pointer to dCSRmat matrix
 * \param iludata   Pointer to ILU_data
 * \param iluparam  Pointer to ILU_param
 *
 * \author Shiquan Zhang
 * \date   12/27/2009
 */
SHORT fasp_ilu_dcsr_setup (dCSRmat *A, 
                           ILU_data *iludata, 
                           ILU_param *iluparam)
{
#if FASP_USE_ILU
    const INT   type=iluparam->ILU_type, print_level=iluparam->print_level;
    const INT   n=A->col, nnz=A->nnz, mbloc=n;
    const REAL  ILU_droptol=iluparam->ILU_droptol;
    const REAL  permtol=iluparam->ILU_permtol;
    
    // local variable
    INT    lfil=iluparam->ILU_lfil, lfilt=iluparam->ILU_lfil;
    INT    ierr, iwk, nzlu, nwork, *ijlu;
    REAL  *luval;
    
    clock_t  setup_start, setup_end;
    REAL     setup_duration;
    SHORT    status = SUCCESS;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_ilu_dcsr_setup ...... [Start]\n");
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",m,n,nnz);
#endif
    
    setup_start=clock();
    
    // Expected amount of memory for ILU needed and allocate memory 
    switch (type) {
        case ILUt:
            iwk=3*nnz;     // iwk is the maxim possible nnz for ILU    
            lfilt=floor(n*0.5)+1;
            break;
        case ILUtp:
            iwk=2*nnz;     // iwk is the maxim possible nnz for ILU    
            lfilt=floor(n*0.5)+1;
            break;
        default: // ILUk
            if (lfil == 0) iwk=nnz+500;
            else iwk=(lfil+2)*nnz;
            break;
    } 
    
    nwork  = 4*n;
    
#if DEBUG_MODE
    printf("### DEBUG: fill-in=%d, iwk=%d, nwork=%d\n", lfil, iwk, nwork);
#endif
    
    // setup ILU preconditioner
    iludata->row=iludata->col=n;    
    fasp_ilu_data_alloc(iwk, nwork, iludata);
    
#if CHMEM_MODE
    printf("### DEBUG: memory usage after fasp_ilu_data_alloc: \n");
    fasp_mem_usage();
#endif
    
    // ILU decomposition
    ijlu=iludata->ijlu;
    luval=iludata->luval;
    
    switch (type) {
        case ILUt:
            ilut_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,luval,ijlu,&iwk,&ierr,&nzlu);
            break;
        case ILUtp:
            ilutp_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,&permtol,
                   &mbloc,luval,ijlu,&iwk,&ierr,&nzlu);
            break;
        default: // ILUk
            iluk_(&n,A->val,A->JA,A->IA,&lfil,luval,ijlu,&iwk,&ierr,&nzlu);
            break;
    } 
    
    fasp_dcsr_shift(A, -1);
    
#if CHMEM_MODE
    printf("### DEBUG: memory usage after ILU setup: \n");
    fasp_mem_usage();
#endif
    
    iludata->nzlu=nzlu;
    iludata->nwork=nwork;
    
#if DEBUG_MODE
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
        setup_end=clock();
        setup_duration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
        
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
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_ilu_dcsr_setup ...... [Finish]\n");
#endif
    
    return status;
    
#else // WITH_ILU
    
    printf("### ERROR: ILU is not enabled!!!\n");
    exit(ERROR_MISC);
    
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
