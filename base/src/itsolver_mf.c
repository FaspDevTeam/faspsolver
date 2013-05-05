/*! \file itsolver_mf.c
 *  \brief Matrix-free iterative solvers
 */

#include <time.h>
#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_itsolver (mxv_matfree *mf, dvector *b, dvector *x, 
 *                               precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods for CSR matrices
 *
 * \param mf       Pointer to mxv_matfree
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date   09/25/2009 
 *
 * \note This is an abstract interface for iterative methods.
 *
 * Modified by Feiteng Huang on 09/19/2012: matrix free
 */
INT fasp_solver_itsolver (mxv_matfree *mf, 
                          dvector *b, 
                          dvector *x, 
                          precond *pc, 
                          itsolver_param *itparam)
{
    const SHORT  print_level   = itparam->print_level;    
    const SHORT  itsolver_type = itparam->itsolver_type;
    const SHORT  stop_type     = itparam->stop_type;
    const SHORT  restart       = itparam->restart;
    const INT    MaxIt         = itparam->maxit;
    const REAL   tol           = itparam->tol; 
    
    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT iter;
    
    fasp_gettime(&solver_start);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_itsolver ...... [Start]\n");
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );
    
    /* Choose a desirable Krylov iterative solver */
    switch ( itsolver_type ) {
        case SOLVER_CG:
            if (print_level>PRINT_NONE) printf("\nCalling PCG solver (matrix-free) ...\n");
            iter = fasp_solver_pcg(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        case SOLVER_BiCGstab:
            if (print_level>PRINT_NONE) printf("\nCalling BiCGstab solver (matrix-free) ...\n");
            iter = fasp_solver_pbcgs(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        case SOLVER_MinRes:
            if (print_level>PRINT_NONE) printf("\nCalling MinRes solver (matrix-free) ...\n");
            iter = fasp_solver_pminres(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        case SOLVER_GMRES:
            if (print_level>PRINT_NONE) printf("\nCalling GMRes solver (matrix-free) ...\n");
            iter = fasp_solver_pgmres(mf, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_VGMRES: 
            if (print_level>PRINT_NONE) printf("\nCalling vGMRes solver (matrix-free) ...\n");
            iter = fasp_solver_pvgmres(mf, b, x, pc, tol, MaxIt, restart, stop_type, print_level);    
            break;
            
        case SOLVER_VFGMRES: 
            if (print_level>PRINT_NONE) printf("\nCalling vFGMRes solver (matrix-free) ...\n");
            iter = fasp_solver_pvfgmres(mf, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_GCG:
            if (print_level>PRINT_NONE) printf("\nCalling GCG solver (matrix-free) ...\n");
            iter = fasp_solver_pgcg(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
            return ERROR_SOLVER_TYPE;
            
    } 
    
    if ( (print_level>=PRINT_SOME) && (iter >= 0) ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_itsolver ...... [Finish]\n");
#endif
    
    return iter;
}    

/**
 * \fn INT fasp_solver_krylov (mxv_matfree *mf, dvector *b, dvector *x, 
 *                             itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods -- without preconditioner 
 *
 * \param mf       Pointer to mxv_matfree
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date   09/25/2009 
 *
 * Modified by Feiteng Huang on 09/20/2012: matrix free
 */
INT fasp_solver_krylov (mxv_matfree *mf, 
                        dvector *b, 
                        dvector *x, 
                        itsolver_param *itparam)
{
    const SHORT print_level = itparam->print_level;
    
    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_krylov ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);
    
    status = fasp_solver_itsolver(mf,b,x,NULL,itparam);
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Krylov method totally", solver_duration);
    }    
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_krylov ...... [Finish]\n");
#endif
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
