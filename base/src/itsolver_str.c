/*! \file itsolver_str.c
 *
 *  \brief Iterative solvers for dSTRmat matrices
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dstr_itsolver (dSTRmat *A, dvector *b, dvector *x, 
 *                                    precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A        Pointer to the coeff matrix in dSTRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date   09/25/2009 
 */
INT fasp_solver_dstr_itsolver (dSTRmat *A, 
                               dvector *b, 
                               dvector *x, 
                               precond *pc, 
                               itsolver_param *itparam)
{
    const SHORT  print_level = itparam->print_level;
    const SHORT  itsolver_type = itparam->itsolver_type;
    const SHORT  stop_type = itparam->stop_type;
    const SHORT  restart = itparam->restart;
    const INT    MaxIt = itparam->maxit;
    const REAL   tol = itparam->tol; 

    // local variables
    INT iter;    
    REAL solver_start, solver_end, solver_duration;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);
#endif

    fasp_gettime(&solver_start);

    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );
    
    switch (itsolver_type) {
    
    case SOLVER_CG: 
        if ( print_level>PRINT_NONE ) printf("\nCalling CG solver (STR) ...\n");
        iter=fasp_solver_dstr_pcg(A, b, x, pc, tol, MaxIt, stop_type, print_level); 
        break;
    
    case SOLVER_BiCGstab:
        if ( print_level>PRINT_NONE ) printf("\nCalling BiCGstab solver (STR) ...\n");
        iter=fasp_solver_dstr_pbcgs(A, b, x, pc, tol, MaxIt, stop_type, print_level); 
        break;
    
    case SOLVER_GMRES:
        if ( print_level>PRINT_NONE ) printf("\nCalling GMRES solver (STR) ...\n");
        iter=fasp_solver_dstr_pgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);    
        break;    
    
    case SOLVER_VGMRES:
        if ( print_level>PRINT_NONE ) printf("\nCalling vGMRES solver (STR) ...\n");
        iter=fasp_solver_dstr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);    
        break;    
    
    default:
        printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
        iter=ERROR_SOLVER_TYPE;
    
    }
    
    if ( (print_level>PRINT_MIN) && (iter >= 0) ) {
        fasp_gettime(&solver_end);    
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return iter;
}    

/**
 * \fn INT fasp_solver_dstr_krylov (dSTRmat *A, dvector *b, dvector *x, 
 *                                  itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A         Pointer to the coeff matrix in dSTRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return          Number of iterations if succeed
 *
 * \author Zhiyang Zhou
 * \date   04/25/2010
 */
INT fasp_solver_dstr_krylov (dSTRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             itsolver_param *itparam)
{
    const INT print_level = itparam->print_level;
    INT status = FASP_SUCCESS;
    REAL solver_start, solver_end, solver_duration;
    
    // solver part
    fasp_gettime(&solver_start);

    status=fasp_solver_dstr_itsolver(A,b,x,NULL,itparam);

    fasp_gettime(&solver_end);
    
    if (print_level>=PRINT_MIN) {
        solver_duration = solver_end - solver_start;
        print_cputime("Krylov method totally", solver_duration);
    }
    
    return status;
}

/**
 * \fn INT fasp_solver_dstr_krylov_diag (dSTRmat *A, dvector *b, dvector *x, 
 *                                       itsolver_param *itparam)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 *
 * \param A         Pointer to the coeff matrix in dSTRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return          Number of iterations if succeed
 *
 * \author Zhiyang Zhou
 * \date   4/23/2010
 */
INT fasp_solver_dstr_krylov_diag (dSTRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam)
{
    const INT print_level = itparam->print_level;
    INT status = FASP_SUCCESS;
    REAL solver_start, solver_end, solver_duration;
    INT nc=A->nc,i;
    INT nc2=nc*nc;
    INT ngrid=A->ngrid;
    
    // setup preconditioner
    precond_diagstr diag;
    fasp_dvec_alloc(ngrid*nc2, &(diag.diag));
    fasp_array_cp(ngrid*nc2, A->diag, diag.diag.val);
    
    diag.nc=nc;
    
    for (i=0;i<ngrid;++i) fasp_blas_smat_inv(&(diag.diag.val[i*nc2]),nc);
    
    precond *pc = (precond *)fasp_mem_calloc(1,sizeof(precond));    
    
    pc->data = &diag;
    pc->fct  = fasp_precond_dstr_diag;
    
    // solver part
    fasp_gettime(&solver_start);

    status=fasp_solver_dstr_itsolver(A,b,x,pc,itparam);

    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if (print_level>=PRINT_MIN) 
        print_cputime("Diag_Krylov method totally", solver_duration);
    
    return status;
}

/**
 * \fn INT fasp_solver_dstr_krylov_ilu(dSTRmat *A, dvector *b, dvector *x, 
 *                                 itsolver_param *itparam, ILU_param *iluparam)
 *
 * \brief Solve Ax=b by structured ILU preconditioned Krylov methods 
 *
 * \param A         Pointer to the coeff matrix in dSTRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param iluparam  Pointer to parameters for ILU
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date   05/01/2010
 */
INT fasp_solver_dstr_krylov_ilu (dSTRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 ILU_param *iluparam)
{
    const INT print_level = itparam->print_level;
    const INT ILU_lfil = iluparam->ILU_lfil;
    INT status = FASP_SUCCESS;
    REAL setup_start, setup_end, solver_start, solver_end, solver_duration, setup_duration;
    
    //set up
    dSTRmat LU;
    
    fasp_gettime(&setup_start);

    if (ILU_lfil == 0) {
        fasp_ilu_dstr_setup0(A,&LU);
    }
    else if (ILU_lfil == 1) {
        fasp_ilu_dstr_setup1(A,&LU);
    }
    else  {    
        printf("### ERROR: Illegal level of fill-in for structured ILU (lfil>=2)!!\n");
        return 0;
    }

    fasp_gettime(&setup_end);
    
    setup_duration = setup_end - setup_start;
    
    if (print_level>PRINT_NONE) 
        printf("structrued ILU(%d) setup costs %f seconds.\n", ILU_lfil, setup_duration);
    
    precond pc; pc.data=&LU;
    if (ILU_lfil == 0) {
        pc.fct = fasp_precond_dstr_ilu0;
    }
    else if (ILU_lfil == 1) {
        pc.fct = fasp_precond_dstr_ilu1;
    }
    else {
        printf("### ERROR: Illegal level of fill-in for structured ILU (lfil>=2)!!\n");
        return 0;
    }
    
    // solver part
    fasp_gettime(&solver_start);

    status=fasp_solver_dstr_itsolver(A,b,x,&pc,itparam);

    fasp_gettime(&solver_end);
        
    if (print_level>=PRINT_MIN) {
        solver_duration = solver_end - solver_start;
        printf("Iterative solver costs %f seconds.\n", solver_duration);
        print_cputime("ILU_Krylov method totally", setup_duration+solver_duration);
    }
    
    fasp_dstr_free(&LU);
    
    return status;
}

/**
 * \fn INT fasp_solver_dstr_krylov_blockgs (dSTRmat *A, dvector *b, dvector *x, 
 *                                          itsolver_param *itparam, ivector *neigh,
 *                                          ivector *order)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 *
 * \param A         Pointer to the coeff matrix in dSTRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param neigh     Pointer to neighbor vector
 * \param order     Pointer to solver ordering
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date   10/10/2010
 */
INT fasp_solver_dstr_krylov_blockgs (dSTRmat *A, 
                                     dvector *b, 
                                     dvector *x, 
                                     itsolver_param *itparam, 
                                     ivector *neigh, 
                                     ivector *order)
{
    // Parameter for iterative method
    const INT print_level = itparam->print_level;
    
    // Information about matrices
    INT ngrid=A->ngrid;
    
    // return parameter
    INT status = FASP_SUCCESS;
    
    // local parameter
    REAL solver_start, solver_end, setup_start, setup_end, solver_duration, setup_duration;
    
    dvector *diaginv;
    ivector *pivot;
    
    // setup preconditioner
    fasp_gettime(&setup_start);
    
    diaginv = (dvector *)fasp_mem_calloc(ngrid, sizeof(dvector));
    pivot = (ivector *)fasp_mem_calloc(ngrid, sizeof(ivector));
    fasp_generate_diaginv_block(A, neigh, diaginv, pivot);
    
    precond_data_str pcdata;
    pcdata.A_str = A;
    pcdata.diaginv = diaginv;
    pcdata.pivot = pivot;
    pcdata.order = order;
    pcdata.neigh = neigh;
    
    precond pc; pc.data = &pcdata; pc.fct = fasp_precond_dstr_blockgs;
    
    fasp_gettime(&setup_end);
    
    if (print_level>PRINT_NONE) {
        setup_duration = setup_end - setup_start;      
        printf("Preconditioner setup costs %f seconds.\n", setup_duration);
    }
    
    // solver part
    fasp_gettime(&solver_start);

    status=fasp_solver_dstr_itsolver(A,b,x,&pc,itparam);

    fasp_gettime(&solver_end);
    
    if (print_level>=PRINT_MIN) {
        solver_duration = solver_end - solver_start;
        printf("Iterative solver costs %f seconds.\n", solver_duration);
        print_cputime("BlockGS_Krylov method totally", setup_duration+solver_duration);
    }
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
