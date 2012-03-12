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
 * \fn int fasp_solver_dbsr_itsolver(dBSRmat *A, dvector *b, dvector *x, 
 *                                   precond *pc, itsolver_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A        pointer to the dBSRmat matrix
 * \param b        pointer to the dvector of right hand side
 * \param x        pointer to the dvector of dofs
 * \param pc     pointer to the preconditioner data
 * \param itparam  pointer to parameters for iterative solvers
 * \return          the number of iterations
 *
 * \author Zhiyang Zhou, Xiaozhe Hu
 * \date 2010/10/26
 */
int fasp_solver_dbsr_itsolver (dBSRmat *A, 
                               dvector *b, 
                               dvector *x, 
                               precond *pc, 
                               itsolver_param *itparam)
{
	const int print_level = itparam->print_level;	
	const int itsolver_type = itparam->itsolver_type;
	const int stop_type = itparam->stop_type;
	const int MaxIt = itparam->maxit;
	const int restart = itparam->restart;
	const double tol = itparam->tol; 
	
	int iter;
	clock_t solver_start=clock();
	
	switch (itsolver_type) {
			
		case SOLVER_BiCGstab:
			if (print_level>0) printf("Calling BiCGstab solver (BSR format) ...\n");
			iter=fasp_solver_dbsr_pbcgs(A, b, x, MaxIt, tol, pc, print_level, stop_type); break;
			
		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRES solver (BSR format) ...\n");
			iter=fasp_solver_dbsr_pgmres(A, b, x, MaxIt, tol, pc, print_level, stop_type, restart);	break;		
			
		case SOLVER_VGMRES:
			if (print_level>0) printf("Calling vGMRES solver (BSR format) ...\n");
			iter=fasp_solver_dbsr_pvgmres(A, b, x, MaxIt, tol, pc, print_level, stop_type, restart); break;	
            
        case SOLVER_VFGMRES: 
			if (print_level>0) printf("Calling vFGMRes solver (BSR format) ...\n");		
			iter = fasp_solver_dbsr_pvfgmres(A, b, x, MaxIt, tol, pc, print_level, stop_type, restart);	break;
			
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
 * \fn int fasp_solver_dbsr_krylov (dBSRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods 
 * \param A:	pointer to the dCSRmat matrix
 * \param b:	pointer to the dvector of right hand side
 * \param x:	pointer to the dvector of dofs
 * \param itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 * \author Zhiyang Zhou, Xiaozhe Hu
 * \date 2010/10/26
 */
int fasp_solver_dbsr_krylov (dBSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             itsolver_param *itparam)
{
	const int print_level = itparam->print_level;
	int status = SUCCESS;
	clock_t solver_start, solver_end;
	double solver_duration;
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dbsr_itsolver(A,b,x,NULL,itparam);
	solver_end=clock();
	
	if (status == ERROR_ALLOC_MEM) goto MEMORY_ERROR;
	
	if (print_level>0) {
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Solver costs %f seconds.\n", solver_duration);
	}
	
	return status;
	
MEMORY_ERROR:
	printf("krylov_bsr: Cannot allocate memory!\n");
	exit(status);	
}

/**
 * \fn int fasp_solver_dbsr_krylov_diag (dBSRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 * \param A:	pointer to the dBSRmat matrix
 * \param b:	pointer to the dvector of right hand side
 * \param x:	pointer to the dvector of dofs
 * \param itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 *\author Zhiyang Zhou, Xiaozhe Hu
 * \date 2010/10/26
 */
int fasp_solver_dbsr_krylov_diag (dBSRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam)
{
	const int print_level = itparam->print_level;
	int status = SUCCESS;
	clock_t solver_start, solver_end;
	double solver_duration;
	
	int nb=A->nb,i,k;
	int nb2=nb*nb;
	int ROW=A->ROW;
	
	// setup preconditioner
	precond_diagbsr diag;
	fasp_dvec_alloc(ROW*nb2, &(diag.diag));
	
	//! get all the diagonal sub-blocks   
	for (i = 0; i < ROW; ++i)
	{
		for (k = A->IA[i]; k < A->IA[i+1]; ++k)
		{
			if (A->JA[k] == i)
				memcpy(diag.diag.val+i*nb2, A->val+k*nb2, nb2*sizeof(double));
		}
	}
	
	diag.nb=nb;
	
	for (i=0; i<ROW; ++i) fasp_blas_smat_inv(&(diag.diag.val[i*nb2]), nb);
	
	precond *pc = (precond *)fasp_mem_calloc(1,sizeof(precond));	
	if (status<0) goto FINISHED;
	
	pc->data = &diag;
	pc->fct  = fasp_precond_dbsr_diag;
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dbsr_itsolver(A,b,x,pc,itparam);
	solver_end=clock();
	
	if (print_level>0) {
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Solver costs %f seconds.\n", solver_duration);
	}	
	
	// clean up
	fasp_dvec_free(&(diag.diag));
	
FINISHED:
	if (status == ERROR_ALLOC_MEM) goto MEMORY_ERROR;
	return status;
	
MEMORY_ERROR:
	printf("krylov_diag_bsr: Cannot allocate memory!\n");
	exit(status);	
	return status;
	
}

/**
 * \fn int fasp_solver_dbsr_krylov_ilu(dBSRmat *A, dvector *b, dvector *x, itsolver_param *itparam, ILU_param *iluparam)
 * \brief Solve Ax=b by ILUs preconditioned Krylov methods
 * \param A	pointer to the dBSRmat matrix
 * \param b	pointer to the dvector of right hand side
 * \param x	pointer to the dvector of dofs
 * \param itparam pointer to parameters for iterative solvers
 * \param iluparam pointer to parameters of ILU
 * \return the number of iterations
 *
 * \author Shiquang Zhang, Xiaozhe Hu
 * \date 10/26/2010
 */
int fasp_solver_dbsr_krylov_ilu (dBSRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 ILU_param *iluparam)
{
	const int print_level = itparam->print_level;	
	clock_t solver_start, solver_end;
	double solver_duration;
	int status = SUCCESS;
	
	ILU_data LU; 
	precond pc; 
	
#if DEBUG_MODE
	printf("krylov_ilu ...... [Start]\n");
	printf("krylov_ilu: matrix size: %d %d %d\n", A->ROW, A->COL, A->NNZ);
	printf("krylov_ilu: rhs/sol size: %d %d\n", b->row, x->row);		
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
	
	if (print_level>0) {
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("ILU(%d)_Krylov solver costs %f seconds.\n", iluparam->ILU_lfil,solver_duration);
	}
	
FINISHED: 
	fasp_ilu_data_free(&LU);
	
#if DEBUG_MODE
	printf("krylov_ilu ...... [Finish]\n");
#endif
	
	if (status == ERROR_ALLOC_MEM) goto MEMORY_ERROR;
	return status;
	
MEMORY_ERROR:
	printf("krylov_ilu_bsr: Cannot allocate memory!\n");
	exit(status);	
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
