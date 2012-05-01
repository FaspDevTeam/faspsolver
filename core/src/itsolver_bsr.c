/*! \file itsolver_bsr.c
 *  \brief Iterative solvers for BSR matrices(main file)
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dbsr_itsolver(dBSRmat *A, dvector *b, dvector *x, 
 *                                   precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param b        Pointer to the dvector of right hand side
 * \param x        Pointer to the dvector of dofs
 * \param pc       Pointer to the preconditioner data
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Zhiyang Zhou, Xiaozhe Hu
 * \date   2010/10/26
 */
INT fasp_solver_dbsr_itsolver (dBSRmat *A, 
                               dvector *b, 
                               dvector *x, 
                               precond *pc, 
                               itsolver_param *itparam)
{
    const SHORT print_level = itparam->print_level;    
    const SHORT itsolver_type = itparam->itsolver_type;
    const SHORT stop_type = itparam->stop_type;
    const SHORT restart = itparam->restart;
    const INT   MaxIt = itparam->maxit;
    const REAL  tol = itparam->tol; 
    
    INT iter;
    clock_t solver_start=clock();
    
    switch (itsolver_type) {
    
    case SOLVER_BiCGstab:
    if (print_level>0) printf("Calling BiCGstab solver (BSR format) ...\n");
    iter=fasp_solver_dbsr_pbcgs(A, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
    
    case SOLVER_GMRES:
    if (print_level>0) printf("Calling GMRES solver (BSR format) ...\n");
    iter=fasp_solver_dbsr_pgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);    
            break;    
    
    case SOLVER_VGMRES:
    if (print_level>0) printf("Calling vGMRES solver (BSR format) ...\n");
    iter=fasp_solver_dbsr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;    
            
        case SOLVER_VFGMRES: 
    if (print_level>0) printf("Calling vFGMRes solver (BSR format) ...\n");    
    iter = fasp_solver_dbsr_pvfgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
    
    default:
    printf("### ERROR: Wrong itertive solver type %d!\n", itsolver_type);
    iter = ERROR_SOLVER_TYPE;
    
    }
    
    if ((print_level>PRINT_MIN) && (iter >= 0)) {
    clock_t solver_end=clock();    
    REAL solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
    printf("Iterative solver costs %f seconds.\n", solver_duration);
    }
    
    return iter;
}

/**
 * \fn INT fasp_solver_dbsr_krylov (dBSRmat *A, dvector *b, dvector *x, 
 *                                  itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A            Pointer to the dCSRmat matrix
 * \param b            Pointer to the dvector of right hand side
 * \param x            Pointer to the dvector of dofs
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return          Number of iterations if succeed
 *
 * \author Zhiyang Zhou, Xiaozhe Hu
 * \date   2010/10/26
 */
INT fasp_solver_dbsr_krylov (dBSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             itsolver_param *itparam)
{
    const SHORT print_level = itparam->print_level;
    INT status = SUCCESS;
    clock_t solver_start, solver_end;
    REAL solver_duration;
    
    // solver part
    solver_start=clock();
    status=fasp_solver_dbsr_itsolver(A,b,x,NULL,itparam);
    solver_end=clock();
    
    if (print_level>PRINT_NONE) {
    solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
    printf("Solver costs %f seconds.\n", solver_duration);
    }
    
    return status;
}

/**
 * \fn INT fasp_solver_dbsr_krylov_diag (dBSRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 *
 * \param A            Pointer to the dBSRmat matrix
 * \param b            Pointer to the dvector of right hand side
 * \param x            Pointer to the dvector of dofs
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return the number of iterations
 *
 * \author Zhiyang Zhou, Xiaozhe Hu
 * \date   2010/10/26
 */
INT fasp_solver_dbsr_krylov_diag (dBSRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam)
{
    const SHORT print_level = itparam->print_level;
    INT status = SUCCESS;
    clock_t solver_start, solver_end;
    REAL solver_duration;
    
    INT nb=A->nb,i,k;
    INT nb2=nb*nb;
    INT ROW=A->ROW;
    
    // setup preconditioner
    precond_diagbsr diag;
    fasp_dvec_alloc(ROW*nb2, &(diag.diag));
    
    // get all the diagonal sub-blocks   
    for (i = 0; i < ROW; ++i)
    {
    for (k = A->IA[i]; k < A->IA[i+1]; ++k)
    {
    if (A->JA[k] == i)
    memcpy(diag.diag.val+i*nb2, A->val+k*nb2, nb2*sizeof(REAL));
    }
    }
    
    diag.nb=nb;
    
    for (i=0; i<ROW; ++i) fasp_blas_smat_inv(&(diag.diag.val[i*nb2]), nb);
    
    precond *pc = (precond *)fasp_mem_calloc(1,sizeof(precond));    
    pc->data = &diag;
    pc->fct  = fasp_precond_dbsr_diag;
    
    // solver part
    solver_start=clock();
    status=fasp_solver_dbsr_itsolver(A,b,x,pc,itparam);
    solver_end=clock();
    
    if (print_level>0) {
    solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
    printf("Solver costs %f seconds.\n", solver_duration);
    }    
    
    // clean up
    fasp_dvec_free(&(diag.diag));

    return status;
}

/**
 * \fn INT fasp_solver_dbsr_krylov_ilu (dBSRmat *A, dvector *b, dvector *x, 
 *                                      itsolver_param *itparam, ILU_param *iluparam)
 *
 * \brief Solve Ax=b by ILUs preconditioned Krylov methods
 *
 * \param A            Pointer to dBSRmat matrix
 * \param b            Pointer to dvector of right hand side
 * \param x            Pointer to dvector of dofs
 * \param itparam   Pointer to parameters for iterative solvers
 * \param iluparam  Pointer to parameters of ILU
 *
 * \return          Number of iterations if succeed
 *
 * \author Shiquang Zhang, Xiaozhe Hu
 * \date   10/26/2010
 */
INT fasp_solver_dbsr_krylov_ilu (dBSRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 ILU_param *iluparam)
{
    const SHORT print_level = itparam->print_level;    
    clock_t solver_start, solver_end;
    REAL solver_duration;
    INT status = SUCCESS;
    
    ILU_data LU; 
    precond pc; 
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dbsr_krylov_ilu ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->ROW, A->COL, A->NNZ);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    // ILU setup for whole matrix
    if ( (status = fasp_ilu_dbsr_setup(A,&LU,iluparam))<0 ) goto FINISHED;
    
    // check iludata
    if ( (status = fasp_mem_iludata_check(&LU))<0 ) goto FINISHED;
    
    // set preconditioner 
    pc.data = &LU; pc.fct = fasp_precond_dbsr_ilu;
    
    // solve
    solver_start=clock();
    status=fasp_solver_dbsr_itsolver(A,b,x,&pc,itparam);
    solver_end=clock();    
    
    if (print_level>PRINT_NONE) {
    solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
    printf("ILU(%d)_Krylov solver costs %f seconds.\n", iluparam->ILU_lfil,solver_duration);
    }
    
FINISHED: 
    fasp_ilu_data_free(&LU);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dbsr_krylov_ilu ...... [Finish]\n");
#endif
    
    if (status < 0) goto MEMORY_ERROR;
    return status;
    
MEMORY_ERROR:
    fasp_chkerr (status, "fasp_solver_dbsr_krylov_ilu");
}

/**
 * \fn INT fasp_solver_dbsr_krylov_amg (dBSRmat *A, dvector *b, dvector *x, 
 *                                      itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve Ax=b by AMG preconditioned Krylov methods 
 *
 * \param A            Pointer to dBSRmat matrix
 * \param b            Pointer to dvector of right hand side
 * \param x            Pointer to dvector of dofs
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters of AMG
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date   03/16/2012
 */
INT fasp_solver_dbsr_krylov_amg(dBSRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam, 
                                  AMG_param *amgparam 
                                )
{
    //--------------------------------------------------------------
    // Part 1: prepare
    // --------------------------------------------------------------
    //! parameters of iterative method
    const SHORT print_level = itparam->print_level;
    const SHORT max_levels = amgparam->max_levels;
    
    //! return variable
    INT status = SUCCESS;
    
    // data of AMG 
    AMG_data_bsr *mgl=fasp_amg_data_bsr_create(max_levels);
    
    // timing
    clock_t setup_start, setup_end, solver_start, solver_end;
    REAL solver_duration, setup_duration;
    
    //--------------------------------------------------------------
    //Part 2: set up the preconditioner
    //--------------------------------------------------------------
    setup_start=clock();
    
    
    //initialize A, b, x for mgl[0]
    mgl[0].A = fasp_dbsr_create(A->ROW, A->COL, A->NNZ, A->nb, A->storage_manner);
    fasp_dbsr_cp(A,  &(mgl[0].A));
    mgl[0].b = fasp_dvec_create(mgl[0].A.ROW*mgl[0].A.nb); 
    mgl[0].x = fasp_dvec_create(mgl[0].A.COL*mgl[0].A.nb); 
    
    status = fasp_amg_setup_ua_bsr(mgl, amgparam);
    if (status < 0) goto FINISHED;
    
    setup_end=clock();
    
    setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
    
    precond_data_bsr precdata;
    precdata.maxit = amgparam->maxit;
    precdata.tol = amgparam->tol;
    precdata.cycle_type = amgparam->cycle_type;
    precdata.smoother = amgparam->smoother;
    precdata.presmooth_iter = amgparam->presmooth_iter;
    precdata.postsmooth_iter = amgparam->postsmooth_iter;
    precdata.coarsening_type = amgparam->coarsening_type;
    precdata.relaxation = amgparam->relaxation;
    precdata.coarse_scaling = amgparam->coarse_scaling;
    precdata.amli_degree = amgparam->amli_degree;
    precdata.amli_coef = amgparam->amli_coef;
    precdata.tentative_smooth = amgparam->tentative_smooth;
    precdata.max_levels = mgl[0].num_levels;
    precdata.mgl_data = mgl;
    precdata.A = A;
    
    if (status < 0) goto FINISHED;
    
    precond prec; 
    prec.data = &precdata; 
    switch (amgparam->cycle_type) {
        case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
            prec.fct = fasp_precond_dbsr_nl_amli; break;
        default: // V,W-Cycle AMG
            prec.fct = fasp_precond_dbsr_amg; break;
    }
    
    //--------------------------------------------------------------
    // Part 3: solver
    //--------------------------------------------------------------
    solver_start=clock();
    status=fasp_solver_dbsr_itsolver(A,b,x,&prec,itparam);
    solver_end=clock();
    
    solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
    
    if (print_level>0) {
    printf("Setup costs %f seconds.\n", setup_duration);
    printf("Iterative solver costs %f seconds.\n", solver_duration);
    printf ("BSR_CPR_krylov method totally costs %f seconds.\n", setup_duration + solver_duration);
    }
    
FINISHED:
    //fasp_amg_data_free(mgl);    // Xiaozhe: need to be added
    if (status == ERROR_ALLOC_MEM) goto MEMORY_ERROR;
    return status;
    
MEMORY_ERROR:
    printf("krylov_CPR_bsr: Cannot allocate memory!\n");
    exit(status);    
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
