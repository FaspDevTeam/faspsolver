/*! \file itsolver_omp.c
 *  \brief Iterative solvers (main file)
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*---------------------------------omp----------------------------------------*/

/**
 * \fn int fasp_solver_dcsr_itsolver_omp(dCSRmat *A, dvector *b, dvector *x, 
 *                                   precond *prec, itsolver_param *itparam, int nthreads, int openmp_holds)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param *A        pointer to the dCSRmat matrix
 * \param *b        pointer to the dvector of right hand side
 * \param *x        pointer to the dvector of dofs
 * \param *prec     pointer to the preconditioner data
 * \param *itparam  pointer to parameters for iterative solvers
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return          number of iterations
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
int fasp_solver_dcsr_itsolver_omp( dCSRmat *A,
                                   dvector *b,
                                   dvector *x,
                                   precond *prec,
                                   itsolver_param *itparam,
                                   int nthreads,
                                   int openmp_holds )
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
			
		case SOLVER_CG:
			if (print_level>0) printf("Calling PCG solver ...\n");
			iter=fasp_solver_dcsr_pcg_omp(A, b, x, MaxIt, tol, prec, print_level, stop_type, nthreads, openmp_holds); break;
			
		case SOLVER_BiCGstab:
			if (print_level>0) printf("Calling BiCGstab solver ...\n");
			iter=fasp_solver_dcsr_pbcgs(A, b, x, MaxIt, tol, prec, print_level, stop_type); break;
			
		case SOLVER_MinRes:
			if (print_level>0) printf("Calling MinRes solver ...\n");		
			iter=fasp_solver_dcsr_pminres(A, b, x, MaxIt, tol, prec, print_level, stop_type); break;
			
		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRes solver ...\n");		
			iter=fasp_solver_dcsr_pgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	break;
			
		case SOLVER_VGMRES: 
			if (print_level>0) printf("Calling vGMRes solver ...\n");		
			iter=fasp_solver_dcsr_pvgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	break;
			
		default:
			printf("Error: wrong itertive solver type %d!\n", itsolver_type);
			return ERROR_SOLVER_TYPE;
			
	} 
	
	if ((print_level>1) && (iter >= 0)) {
		double solver_end=omp_get_wtime();	
		double solver_duration = solver_end - solver_start;
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
#endif
	return iter;
}

/**
 * \fn int fasp_solver_dcsr_krylov_amg_omp(dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam, AMG_param *amgparam, int nthreads, int openmp_holds)
 * \brief Solve Ax=b by preconditioned Krylov methods with AMG as precondition
 *
 * \param *A        pointer to the dCSRmat matrix
 * \param *b        pointer to the dvector of right hand side
 * \param *x        pointer to the dvector of dofs
 * \param *itparam  pointer to parameters for iterative solvers
 * \param *amgparam pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return          number of iterations
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
int fasp_solver_dcsr_krylov_amg_omp (dCSRmat *A,
                                     dvector *b,
                                     dvector *x,
                                     itsolver_param *itparam,
                                     AMG_param *amgparam,
                                     int nthreads,
                                     int openmp_holds)
{
	int status = SUCCESS;
#if FASP_USE_OPENMP
	const int print_level = itparam->print_level;
	const int max_levels = amgparam->max_levels;
	const int nnz=A->nnz, m=A->row, n=A->col;	
	double solver_start, solver_end;
	double solver_duration;
	
#if DEBUG_MODE
	printf("krylov_amg ...... [Start]\n");
	printf("krylov_amg: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("krylov_amg: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
	
	solver_start=omp_get_wtime();
	
	// initialize A, b, x for mgl[0]		
	AMG_data *mgl=fasp_amg_data_create(max_levels);
	mgl[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp_omp(A, &mgl[0].A, nthreads, openmp_holds);
	mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n); 
	
	// setup preconditioner  
	switch (amgparam->AMG_type) {
		default: // Classical AMG
			status = fasp_amg_setup_rs_omp(mgl, amgparam, nthreads, openmp_holds);
			break;
	}
	
#if CHMEM_MODE	
	fasp_mem_usage();
#endif
	
	if (status < 0) goto FINISHED;
	
	// solver part
	precond_data precdata;
	precdata.max_iter = amgparam->max_iter;
	precdata.tol = amgparam->tol;
	precdata.cycle_type = amgparam->cycle_type;
	precdata.smoother = amgparam->smoother;
	precdata.smooth_order = amgparam->smooth_order;
	precdata.presmooth_iter  = amgparam->presmooth_iter;
	precdata.postsmooth_iter = amgparam->postsmooth_iter;
	precdata.coarsening_type = amgparam->coarsening_type;
	precdata.relaxation = amgparam->relaxation;
	precdata.coarse_scaling = amgparam->coarse_scaling;
	precdata.amli_degree = amgparam->amli_degree;
	precdata.amli_coef = amgparam->amli_coef;
	precdata.tentative_smooth = amgparam->tentative_smooth;
	precdata.max_levels = mgl[0].num_levels;
	precdata.mgl_data = mgl;
	
	precond prec;
	prec.data = &precdata; 
	if (itparam->precond_type == PREC_FMG)
	{
		//prec.fct = fasp_precond_famg;
	}
	else{
		if (amgparam->cycle_type == AMLI_CYCLE)
		{
//			prec.fct = fasp_precond_amli;
		}
		else
		{
			prec.fct_omp = fasp_precond_amg_omp;
		}
	}
	
	status=fasp_solver_dcsr_itsolver_omp(A,b,x,&prec,itparam,nthreads,openmp_holds);
	
	if (print_level>0) {
		solver_end=omp_get_wtime();	
		solver_duration = solver_end - solver_start;
		printf("AMG_Krylov method totally costs %f seconds.\n", solver_duration);
	}
	
FINISHED:
	fasp_amg_data_free(mgl);
	
#if DEBUG_MODE
	printf("krylov_amg ...... [Finish]\n");
#endif
	
#endif
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
