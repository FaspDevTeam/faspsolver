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
 * \fn int fasp_solver_dcsr_itsolver(dCSRmat *A, dvector *b, dvector *x, 
 *                                   precond *prec, itsolver_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param *A        pointer to the block dCSRmat matrix
 * \param *b        pointer to the dvector of right hand side
 * \param *x        pointer to the dvector of dofs
 * \param *prec     pointer to the preconditioner data
 * \param *itparam  pointer to parameters for iterative solvers
 * \return          the number of iterations
 *
 * \author Chensong Zhang
 * \date 11/25/2010
 */
int fasp_solver_bdcsr_itsolver(block_dCSRmat *A, 
															 dvector *b, 
															 dvector *x, 
															 precond *prec, 
															 itsolver_param *itparam)
{
	const int print_level = itparam->print_level;
	const int itsolver_type = itparam->itsolver_type;
	const int stop_type = itparam->stop_type;
	const double tol = itparam->tol; 
	const int MaxIt = itparam->maxit;
	const int restart = itparam->restart;
	
	clock_t solver_start=clock();
	int iter;
	
	switch (itsolver_type) {
	
	 if (itsolver_type == SOLVER_BiCGstab) {
		if (print_level>0) printf("BiCGstab method (Block CSR format) ...\n");
		iter=fasp_solver_bdcsr_pbcgs(A, b, x, MaxIt, tol, prec, print_level, stop_type);			
	}

		case SOLVER_MinRes:
			if (print_level>0) printf("Calling MinRes solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pminres(A, b, x, MaxIt, tol, prec, print_level, stop_type); break;		

		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRES solver (Block CSR format) ...\n");
			iter=fasp_solver_bdcsr_pgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart); break;			

		default:
			printf("Error: wrong itertive solver type %d!\n", itsolver_type);
			iter = ERROR_SOLVER_TYPE;

	}
	
	if ((print_level>1) && (iter >= 0)) {
		clock_t solver_end=clock();	
		double solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
	
	return iter;
}	

/**
 * \fn int fasp_solver_bdcsr_krylov (block_dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods 
 * \param *A:	pointer to the block_dCSRmat matrix
 * \param *b:	pointer to the dvector of right hand side
 * \param *x:	pointer to the dvector of dofs
 * \param *itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \date 07/18/2010
 */
int fasp_solver_bdcsr_krylov (block_dCSRmat *A, 
															dvector *b, 
															dvector *x, 
															itsolver_param *itparam)
{
	const int print_level = itparam->print_level;
	int status=SUCCESS;
	clock_t solver_start, solver_end;
	double solver_duration;
	
	// solver part
	solver_start=clock();
	status=fasp_solver_bdcsr_itsolver(A,b,x,NULL,itparam);
	solver_end=clock();
	
	solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
	
	if (print_level>0) printf("Solver costs %f seconds.\n", solver_duration);
	
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
