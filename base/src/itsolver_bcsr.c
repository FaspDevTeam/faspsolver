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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
