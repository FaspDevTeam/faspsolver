/*! \file itsolver_str.c
 *  \brief Iterative solvers (main file)
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dstr_itsolver(dSTRmat *A, dvector *b, dvector *x, 
 *                                   precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A        pointer to the dSTRmat matrix
 * \param b        pointer to the dvector of right hand side
 * \param x        pointer to the dvector of dofs
 * \param pc     pointer to the preconditioner data
 * \param itparam  pointer to parameters for iterative solvers
 *
 * \return          the number of iterations
 *
 * \author Chensong Zhang
 * \date 09/25/2009 
 */
INT fasp_solver_dstr_itsolver(dSTRmat *A, 
                              dvector *b, 
                              dvector *x, 
                              precond *pc, 
                              itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;
	const INT itsolver_type = itparam->itsolver_type;
	const INT stop_type = itparam->stop_type;
	const REAL tol = itparam->tol; 
	const INT MaxIt = itparam->maxit;
	const INT restart = itparam->restart;
	INT iter;
	
	clock_t solver_start=clock();
	
	switch (itsolver_type) {
			
		case SOLVER_CG: 
			if (print_level>PRINT_NONE) printf("Calling CG solver (STR format) ...\n");
			iter=fasp_solver_dstr_pcg(A, b, x, MaxIt, tol, pc, print_level, stop_type); 
			break;
			
		case SOLVER_BiCGstab:
			if (print_level>PRINT_NONE) printf("Calling BiCGstab solver (STR format) ...\n");
			iter=fasp_solver_dstr_pbcgs(A, b, x, MaxIt, tol, pc, print_level, stop_type); 
			break;
			
		case SOLVER_GMRES:
			if (print_level>PRINT_NONE) printf("Calling GMRES solver (STR format) ...\n");
			iter=fasp_solver_dstr_pgmres(A, b, x, MaxIt, tol, pc, print_level, stop_type, restart);	
			break;		
			
		case SOLVER_VGMRES:
			if (print_level>PRINT_NONE) printf("Calling vGMRES solver (STR format) ...\n");
			iter=fasp_solver_dstr_pvgmres(A, b, x, MaxIt, tol, pc, print_level, stop_type, restart);	
			break;	
			
		default:
			printf("Error: wrong itertive solver type %d!\n", itsolver_type);
			iter=ERROR_SOLVER_TYPE;
			
	}
	
	if ((print_level>PRINT_MIN) && (iter >= 0)) {
		clock_t solver_end=clock();	
		REAL solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
	
	return iter;
}	

/**
 * \fn INT fasp_solver_dstr_krylov (dSTRmat *A, dvector *b, dvector *x, 
 *                                  itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A	pointer to the dCSRmat matrix
 * \param b	pointer to the dvector of right hand side
 * \param x	pointer to the dvector of dofs
 * \param itparam pointer to parameters for iterative solvers
 *
 * \return the number of iterations
 *
 * \author Zhiyang Zhou
 * \date 4/25/2010
 */
INT fasp_solver_dstr_krylov (dSTRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;
	INT status = SUCCESS;
	clock_t solver_start, solver_end;
	REAL solver_duration;
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,NULL,itparam);
	solver_end=clock();
	
	solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
	
	if (print_level>PRINT_NONE) printf("Solver costs %f seconds.\n", solver_duration);
	
	return status;
}

/**
 * \fn INT fasp_solver_dstr_krylov_diag (dSTRmat *A, dvector *b, dvector *x, 
 *                                       itsolver_param *itparam)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 *
 * \param A:	pointer to the dSTRmat matrix
 * \param b:	pointer to the dvector of right hand side
 * \param x:	pointer to the dvector of dofs
 * \param itparam: pointer to parameters for iterative solvers
 *
 * \return the number of iterations
 *
 * \author Zhiyang Zhou
 * \date 4/23/2010
 */
INT fasp_solver_dstr_krylov_diag (dSTRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;
	INT status = SUCCESS;
	clock_t solver_start, solver_end;
	REAL solver_duration;
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
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,pc,itparam);
	solver_end=clock();
	
	solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
	
	if (print_level>PRINT_NONE) printf("Solver costs %f seconds.\n", solver_duration);
	
	return status;
}

/**
 * \fn INT fasp_solver_dstr_krylov_ilu(dSTRmat *A, dvector *b, dvector *x, 
 *                                 itsolver_param *itparam, ILU_param *iluparam)
 *
 * \brief Solve Ax=b by structured ILU preconditioned Krylov methods 
 *
 * \param A	pointer to the dSTRmat matrix
 * \param b	pointer to the dvector of right hand side
 * \param x	pointer to the dvector of dofs
 * \param itparam pointer to parameters for iterative solvers
 * \param iluparam pointer to parameters for ILU
 *
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \date 05/01/2010
 */
INT fasp_solver_dstr_krylov_ilu (dSTRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 ILU_param *iluparam)
{
	const INT print_level = itparam->print_level;
	const INT ILU_lfil = iluparam->ILU_lfil;
	INT status = SUCCESS;
	clock_t setup_start, setup_end, solver_start, solver_end;
	REAL solver_duration, setup_duration;
	
	//set up
	dSTRmat LU;
	
	setup_start=clock();
	if (ILU_lfil == 0) 
	{
		fasp_ilu_dstr_setup0(A,&LU);
	}
	else if (ILU_lfil == 1)		
	{
		fasp_ilu_dstr_setup1(A,&LU);
	}
	else 
	{		
		printf("krylov_ilu_str: illegal level of fill-in for structured ILU (lfil>=2)!!\n");
		return 0;
	}
	setup_end=clock();
	
	setup_duration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
	
	if (print_level>PRINT_NONE) printf("structrued ILU(%d) setup costs %f seconds.\n", ILU_lfil, setup_duration);
	
	precond pc; pc.data=&LU;
	if (ILU_lfil == 0)
	{
		pc.fct = fasp_precond_dstr_ilu0;
	}
	else if (ILU_lfil == 1)
	{
		pc.fct = fasp_precond_dstr_ilu1;
	}
	else
	{
		printf("krylov_ilu_str: illegal level of fill-in for structured ILU (lfil>=2)!!\n");
		return 0;
	}
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,&pc,itparam);
	solver_end=clock();
	
	solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
	
	if (print_level>PRINT_NONE) printf("Iterative solver costs %f seconds.\n", solver_duration);
	
	if (print_level>PRINT_NONE) printf("krylov_ilu_str method totally costs %f seconds.\n", setup_duration + solver_duration);
	
	fasp_dstr_free(&LU);
	
	return status;
}

/**
 * \fn INT fasp_solver_dstr_krylov_blockgs(dSTRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 *
 * \param A:	pointer to the dSTRmat matrix
 * \param b:	pointer to the dvector of right hand side
 * \param x:	pointer to the dvector of dofs
 * \param itparam: pointer to parameters for iterative solvers
 *
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \date 10/10/2010
 */
INT fasp_solver_dstr_krylov_blockgs (dSTRmat *A, 
                                     dvector *b, 
                                     dvector *x, 
                                     itsolver_param *itparam, 
                                     ivector *neigh, 
                                     ivector *order)
{
	//! Parameter for iterative method
	const INT print_level = itparam->print_level;
	
	//! Information about matrices
	INT ngrid=A->ngrid;
	
	//! return parameter
	INT status = SUCCESS;
	
	//! local parameter
	clock_t solver_start, solver_end, setup_start, setup_end;
	REAL solver_duration, setup_duration;
	
	dvector *diaginv;
	ivector *pivot;
	
	// setup preconditioner
	setup_start=clock();
	
	diaginv = (dvector *)fasp_mem_calloc(ngrid, sizeof(dvector));
	pivot = (ivector *)fasp_mem_calloc(ngrid, sizeof(ivector));
	fasp_generate_diaginv_block(A, neigh, diaginv, pivot);
	
	precond_data_str pcdata;
	pcdata.A_str = A;
	pcdata.diaginv = diaginv;
	pcdata.pivot = pivot;
	pcdata.order = order;
	pcdata.neigh = neigh;
	
	precond pc; pc.data = &pcdata; pc.fct  = fasp_precond_dstr_blockgs;
	
	setup_end=clock();
	
	if (print_level>PRINT_NONE) {
		setup_duration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);	  
		printf("setup costs %f seconds.\n", setup_duration);
	}
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,&pc,itparam);
	solver_end=clock();
	
	
	if (print_level>PRINT_NONE) {
		solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
		printf("solve costs %f seconds.\n", solver_duration);
		printf("krylov_blockGS_str costs %f seconds.\n", solver_duration+setup_duration);
	}
	
	return status;
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
