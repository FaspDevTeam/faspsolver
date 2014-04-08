/*! \file itsolver_bcsr.c
 *
 *  \brief Iterative solvers for block_dCSRmat matrices
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_bdcsr_itsolver (block_dCSRmat *A, dvector *b, dvector *x, 
 *                                     precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods 
 *
 * \param A        Pointer to the coeff matrix in block_dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date   11/25/2010
 */
INT fasp_solver_bdcsr_itsolver (block_dCSRmat *A, 
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

    REAL  solver_start, solver_end, solver_duration;        
    
    fasp_gettime(&solver_start);

    INT iter;
    
    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );

    switch (itsolver_type) {
            
    case SOLVER_BiCGstab:
        if ( print_level>PRINT_NONE ) printf("\nCalling BiCGstab solver (Block CSR) ...\n");
        iter=fasp_solver_bdcsr_pbcgs(A, b, x, pc, tol, MaxIt, stop_type, print_level); 
        break;    
            
    case SOLVER_MinRes:
        if ( print_level>PRINT_NONE ) printf("\nCalling MinRes solver (Block CSR) ...\n");
        iter=fasp_solver_bdcsr_pminres(A, b, x, pc, tol, MaxIt, stop_type, print_level); 
        break;    
            
    case SOLVER_GMRES:
        if ( print_level>PRINT_NONE ) printf("\nCalling GMRES solver (Block CSR) ...\n");
        iter=fasp_solver_bdcsr_pgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level); 
        break;
            
    case SOLVER_VGMRES:
        if (print_level>0) printf("Calling vGMRES solver (Block CSR format) ...\n");
        iter=fasp_solver_bdcsr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
        break;
        
    case SOLVER_VFGMRES:
        if (print_level>0) printf("Calling FGMRES solver (Block CSR format) ...\n");
        iter=fasp_solver_bdcsr_pvfgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
        break;
            
    default:
        printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
        iter = ERROR_SOLVER_TYPE;
            
    }
    
    if ( (print_level>=PRINT_MIN) && (iter >= 0) ) {
        fasp_gettime(&solver_end);    
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
    return iter;
}    

/**
 * \fn INT fasp_solver_bdcsr_krylov (block_dCSRmat *A, dvector *b, dvector *x, 
 *                                   itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date   07/18/2010
 */
INT fasp_solver_bdcsr_krylov (block_dCSRmat *A, 
                              dvector *b, 
                              dvector *x, 
                              itsolver_param *itparam)
{
    const INT print_level = itparam->print_level;
    INT status=SUCCESS;
    REAL solver_start, solver_end, solver_duration;
    
    // solver part
    fasp_gettime(&solver_start);

    status=fasp_solver_bdcsr_itsolver(A,b,x,NULL,itparam);
    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( print_level>=PRINT_MIN ) 
        print_cputime("Krylov method totally", solver_duration);
    
    return status;
}

/**
 * \fn INT fasp_solver_bdcsr_krylov_block(block_dCSRmat *A, dvector *b, dvector *x,
 *                                   itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date   04/07/2014
 */
INT fasp_solver_bdcsr_krylov_block(block_dCSRmat *A,
                              dvector *b,
                              dvector *x,
                              itsolver_param *itparam,
                                AMG_param *amgparam)
{
    const INT print_level = itparam->print_level;
    const INT precond_type = itparam->precond_type;
    INT status=SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;
    
    dCSRmat *A11, *A22, *A33;
    
    const INT max_levels = amgparam->max_levels;
    INT m, n, nnz;
    
    /* setup preconditioner */
    fasp_gettime(&setup_start);
    
    // AMG for A11 block
    AMG_data *mgl1=fasp_amg_data_create(max_levels);
    A11 = A->blocks[0];
    m = A11->row; n = A11->col; nnz = A11->nnz;
    mgl1[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A11,&mgl1[0].A);
    mgl1[0].b=fasp_dvec_create(n); mgl1[0].x=fasp_dvec_create(n);
    
    // setup preconditioner
    switch (amgparam->AMG_type) {
        case SA_AMG: // Smoothed Aggregation AMG
            status = fasp_amg_setup_sa(mgl1, amgparam); break;
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = fasp_amg_setup_ua(mgl1, amgparam); break;
        default: // Classical AMG
            status = fasp_amg_setup_rs(mgl1, amgparam); break;
    }
    
    // AMG for A22 block
    AMG_data *mgl2=fasp_amg_data_create(max_levels);
    A22 = A->blocks[4];
    m = A22->row; n = A22->col; nnz = A22->nnz;
    mgl2[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A22,&mgl2[0].A);
    mgl2[0].b=fasp_dvec_create(n); mgl2[0].x=fasp_dvec_create(n);
    
    // setup preconditioner
    switch (amgparam->AMG_type) {
        case SA_AMG: // Smoothed Aggregation AMG
            status = fasp_amg_setup_sa(mgl2, amgparam); break;
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = fasp_amg_setup_ua(mgl2, amgparam); break;
        default: // Classical AMG
            status = fasp_amg_setup_rs(mgl2, amgparam); break;
    }
    
    // AMG for A33 block
    AMG_data *mgl3=fasp_amg_data_create(max_levels);
    A33 = A->blocks[8];
    m = A33->row; n = A33->col; nnz = A33->nnz;
    mgl3[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A33,&mgl3[0].A);
    mgl3[0].b=fasp_dvec_create(n); mgl3[0].x=fasp_dvec_create(n);
    
    // setup preconditioner
    switch (amgparam->AMG_type) {
        case SA_AMG: // Smoothed Aggregation AMG
            status = fasp_amg_setup_sa(mgl3, amgparam); break;
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = fasp_amg_setup_ua(mgl3, amgparam); break;
        default: // Classical AMG
            status = fasp_amg_setup_rs(mgl3, amgparam); break;
    }
    
    precond_block_data_3 precdata;
    precdata.Abcsr = A;
    precdata.amgparam = amgparam;
    precdata.mgl1 = mgl1;
    precdata.mgl2 = mgl3;
    precdata.mgl3 = mgl3;
    precdata.r = fasp_dvec_create(b->row);
    
    precond prec; prec.data = &precdata;
    
    switch (precond_type)
    {
        case 21:
            prec.fct = fasp_precond_block_diag;
            break;
            
        case 22:
            prec.fct = fasp_precond_block_lower;
            break;
    }
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }

    
    // solver part
    fasp_gettime(&solver_start);
    
    status=fasp_solver_bdcsr_itsolver(A,b,x, &prec,itparam);
    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( print_level>=PRINT_MIN )
        print_cputime("Krylov method totally", solver_duration);
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
