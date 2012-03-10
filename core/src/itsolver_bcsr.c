/*! \file itsolver_bcsr.c
 *  \brief Iterative solvers for block_CSR matrices (main file)
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_bdcsr_itsolver (block_dCSRmat *A, dvector *b, dvector *x, 
 *                                    precond *prec, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A        pointer to the block dCSRmat matrix
 * \param b        pointer to the dvector of right hand side
 * \param x        pointer to the dvector of dofs
 * \param prec     pointer to the preconditioner data
 * \param itparam  pointer to parameters for iterative solvers
 *
 * \return          the number of iterations
 *
 * \author Chensong Zhang
 * \date 11/25/2010
 */
INT fasp_solver_bdcsr_itsolver (block_dCSRmat *A, 
                                dvector *b, 
                                dvector *x, 
                                precond *prec, 
                                itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;
	const INT itsolver_type = itparam->itsolver_type;
	const INT stop_type = itparam->stop_type;
	const REAL tol = itparam->tol; 
	const INT MaxIt = itparam->maxit;
	const INT restart = itparam->restart;
	
	clock_t solver_start=clock();
	INT iter;
	
	switch (itsolver_type) {
            
        case SOLVER_BiCGstab:
            if (print_level>PRINT_NONE) printf("BiCGstab method (Block CSR format) ...\n");
            iter=fasp_solver_bdcsr_pbcgs(A, b, x, MaxIt, tol, prec, print_level, stop_type);			
            
		case SOLVER_MinRes:
			if (print_level>PRINT_NONE) printf("Calling MinRes solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pminres(A, b, x, MaxIt, tol, prec, print_level, stop_type); break;		
            
		case SOLVER_GMRES:
			if (print_level>PRINT_NONE) printf("Calling GMRES solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart); break;			
            
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
 * \fn INT fasp_solver_bdcsr_krylov (block_dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 * \param A:	pointer to the block_dCSRmat matrix
 * \param b:	pointer to the dvector of right hand side
 * \param x:	pointer to the dvector of dofs
 * \param itparam: pointer to parameters for iterative solvers
 *
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \date 07/18/2010
 */
INT fasp_solver_bdcsr_krylov (block_dCSRmat *A, 
                              dvector *b, 
                              dvector *x, 
                              itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;
	INT status=SUCCESS;
	clock_t solver_start, solver_end;
	REAL solver_duration;
	
	// solver part
	solver_start=clock();
	status=fasp_solver_bdcsr_itsolver(A,b,x,NULL,itparam);
	solver_end=clock();
	
	solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
	
	if (print_level>PRINT_NONE) printf("Solver costs %f seconds.\n", solver_duration);
	
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
