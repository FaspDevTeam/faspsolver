/*! \file itsolver_bsr.c
 *  \brief Iterative solvers for BSR matrices(main file)
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
/*-----------------------------------omp--------------------------------------*/

#if FASP_USE_OPENMP

/**
 * \fn void fasp_set_GS_threads_omp(int threads,int its)
 * \brief set threads for CPR. Please add it at the begin of Krylov openmp method function and after iter++.
 *
 * \param threads       total threads of sovler
 * \param its           curent its of the Krylov methods
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/20/2011
 */
void fasp_set_GS_threads_omp(int mythreads, int its)
{
#if 1

	if (its <=8) {
		THDs_AMG_GS =  mythreads;
		THDs_CPR_lGS = mythreads ;
		THDs_CPR_gGS = mythreads ;
	}
	else if (its <=12) {
		THDs_AMG_GS =  mythreads;
		THDs_CPR_lGS = (6 < mythreads) ? 6 : mythreads;
		THDs_CPR_gGS = (4 < mythreads) ? 4 : mythreads;
	}
	else if (its <=15) {
		THDs_AMG_GS =  (4 < mythreads) ? 4 : mythreads;
		THDs_CPR_lGS = (4 < mythreads) ? 4 : mythreads;
		THDs_CPR_gGS = (2 < mythreads) ? 2 : mythreads;
     }
	else if (its <=18) {
		THDs_AMG_GS =  (2 < mythreads) ? 2 : mythreads;
		THDs_CPR_lGS = (2 < mythreads) ? 2 : mythreads;
		THDs_CPR_gGS = (1 < mythreads) ? 1 : mythreads;
	}
	else {
		THDs_AMG_GS =  1;
		THDs_CPR_lGS = 1;
		THDs_CPR_gGS = 1;
	}

#else

	THDs_AMG_GS =  mythreads;
	THDs_CPR_lGS = mythreads ;
	THDs_CPR_gGS = mythreads ;

#endif
}

/**
 * \fn int fasp_solver_dbsr_itsolver_omp(dBSRmat *A, dvector *b, dvector *x, 
 *                                   precond *prec, itsolver_param *itparam, int nthreads, int openmp_holds)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A        pointer to the dBSRmat matrix
 * \param b        pointer to the dvector of right hand side
 * \param x        pointer to the dvector of dofs
 * \param prec     pointer to the preconditioner data
 * \param itparam  pointer to parameters for iterative solvers
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return          the number of iterations
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
int fasp_solver_dbsr_itsolver_omp(dBSRmat *A,
                                  dvector *b,
                                  dvector *x,
                                  precond *prec,
                                  itsolver_param *itparam,
                                  int nthreads,
                                  int openmp_holds)
{
	int iter = 0;
	
#if FASP_USE_OPENMP
	const int print_level = itparam->print_level;	
	const int itsolver_type = itparam->itsolver_type;
	const int stop_type = itparam->stop_type;
	const int MaxIt = itparam->maxit;
	const int restart = itparam->restart;
	const double tol = itparam->tol; 
	
	double solver_start=omp_get_wtime();
	
	switch (itsolver_type) {
			
		case SOLVER_BiCGstab:
			if (print_level>0) printf("Calling BiCGstab solver (BSR format) ...\n");
			iter=fasp_solver_dbsr_pbcgs(A, b, x, MaxIt, tol, prec, print_level, stop_type); break;
			
		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRES solver (BSR format) ...\n");
			iter=fasp_solver_dbsr_pgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	break;		
			
		case SOLVER_VGMRES:
			if (print_level>0) printf("Calling vGMRES solver (BSR format) ...\n");
			iter=fasp_solver_dbsr_pvgmres_omp(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart, nthreads, openmp_holds); break;
			
		default:
			printf("Error: wrong itertive solver type %d!\n", itsolver_type);
			iter = ERROR_SOLVER_TYPE;
			
	}
	
	if ((print_level>1) && (iter >= 0)) {
		double solver_end=omp_get_wtime();	
		double solver_duration = solver_end - solver_start;
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
	
#endif
	return iter;
}

#endif // FASP_USE_OPENMP

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
