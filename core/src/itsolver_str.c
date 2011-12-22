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
 * \fn int fasp_solver_dstr_itsolver(dSTRmat *A, dvector *b, dvector *x, 
 *                                   precond *prec, itsolver_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param *A        pointer to the dSTRmat matrix
 * \param *b        pointer to the dvector of right hand side
 * \param *x        pointer to the dvector of dofs
 * \param *prec     pointer to the preconditioner data
 * \param *itparam  pointer to parameters for iterative solvers
 * \return          the number of iterations
 *
 * \author Chensong Zhang
 * \date 09/25/2009 
 */
int fasp_solver_dstr_itsolver(dSTRmat *A, 
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
	int iter;
	
	clock_t solver_start=clock();
	
	switch (itsolver_type) {
			
		case SOLVER_CG: 
			if (print_level>0) printf("Calling CG solver (STR format) ...\n");
			iter=fasp_solver_dstr_pcg(A, b, x, MaxIt, tol, prec, print_level, stop_type); 
			break;
			
		case SOLVER_BiCGstab:
			if (print_level>0) printf("Calling BiCGstab solver (STR format) ...\n");
			iter=fasp_solver_dstr_pbcgs(A, b, x, MaxIt, tol, prec, print_level, stop_type); 
			break;
			
		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRES solver (STR format) ...\n");
			iter=fasp_solver_dstr_pgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	
			break;		
			
		case SOLVER_VGMRES:
			if (print_level>0) printf("Calling vGMRES solver (STR format) ...\n");
			iter=fasp_solver_dstr_pvgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	
			break;	
			
		default:
			printf("Error: wrong itertive solver type %d!\n", itsolver_type);
			iter=ERROR_SOLVER_TYPE;
			
	}
	
	if ((print_level>1) && (iter >= 0)) {
		clock_t solver_end=clock();	
		double solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
	
	return iter;
}	

/**
 * \fn int fasp_solver_dstr_krylov (dSTRmat *A, dvector *b, dvector *x, 
 *                                  itsolver_param *itparam)
 * \brief Solve Ax=b by standard Krylov methods 
 * \param *A:	pointer to the dCSRmat matrix
 * \param *b:	pointer to the dvector of right hand side
 * \param *x:	pointer to the dvector of dofs
 * \param *itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 * \author Zhiyang Zhou
 * \data 4/25/2010
 */
int fasp_solver_dstr_krylov (dSTRmat *A, 
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
	status=fasp_solver_dstr_itsolver(A,b,x,NULL,itparam);
	solver_end=clock();
	
	solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
	
	if (print_level>0) printf("Solver costs %f seconds.\n", solver_duration);
	
	return status;
}

/**
 * \fn int fasp_solver_dstr_krylov_diag (dSTRmat *A, dvector *b, dvector *x, 
 *                                       itsolver_param *itparam)
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 * \param *A:	pointer to the dSTRmat matrix
 * \param *b:	pointer to the dvector of right hand side
 * \param *x:	pointer to the dvector of dofs
 * \param *itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 * \author Zhiyang Zhou
 * \data 4/23/2010
 */
int fasp_solver_dstr_krylov_diag (dSTRmat *A, 
																	dvector *b, 
																	dvector *x, 
																	itsolver_param *itparam)
{
	const int print_level = itparam->print_level;
	int status = SUCCESS;
	clock_t solver_start, solver_end;
	double solver_duration;
	int nc=A->nc,i;
	int nc2=nc*nc;
	int ngrid=A->ngrid;
	
	// setup preconditioner
	precond_diagstr diag;
	fasp_dvec_alloc(ngrid*nc2, &(diag.diag));
	fasp_array_cp(ngrid*nc2, A->diag, diag.diag.val);
	
	diag.nc=nc;
	
	for (i=0;i<ngrid;++i) fasp_blas_smat_inv(&(diag.diag.val[i*nc2]),nc);
	
	precond *prec = (precond *)fasp_mem_calloc(1,sizeof(precond));	
	
	prec->data = &diag;
	prec->fct  = fasp_precond_dstr_diag;
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,prec,itparam);
	solver_end=clock();
	
	solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
	
	if (print_level>0) printf("Solver costs %f seconds.\n", solver_duration);
	
	return status;
}

/**
 * \fn int fasp_solver_dstr_krylov_ilu(dSTRmat *A, dvector *b, dvector *x, 
 *                                 itsolver_param *itparam, ILU_param *iluparam)
 * \brief Solve Ax=b by structured ILU preconditioned Krylov methods 
 * \param *A:	pointer to the dSTRmat matrix
 * \param *b:	pointer to the dvector of right hand side
 * \param *x:	pointer to the dvector of dofs
 * \param *itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \data 05/01/2010
 */
int fasp_solver_dstr_krylov_ilu (dSTRmat *A, 
																 dvector *b, 
																 dvector *x, 
																 itsolver_param *itparam, 
																 ILU_param *iluparam)
{
	const int print_level = itparam->print_level;
	const int ILU_lfil = iluparam->ILU_lfil;
	int status = SUCCESS;
	clock_t setup_start, setup_end, solver_start, solver_end;
	double solver_duration, setup_duration;
	
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
	
	setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
	
	if (print_level>0) printf("structrued ILU(%d) setup costs %f seconds.\n", ILU_lfil, setup_duration);
	
	precond prec; prec.data=&LU;
	if (ILU_lfil == 0)
	{
		prec.fct = fasp_precond_dstr_ilu0;
	}
	else if (ILU_lfil == 1)
	{
		prec.fct = fasp_precond_dstr_ilu1;
	}
	else
	{
		printf("krylov_ilu_str: illegal level of fill-in for structured ILU (lfil>=2)!!\n");
		return 0;
	}
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,&prec,itparam);
	solver_end=clock();
	
	solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
	
	if (print_level>0) printf("Iterative solver costs %f seconds.\n", solver_duration);
	
	if (print_level>0) printf ("Struct_ILU_krylov method totally costs %f seconds.\n", setup_duration + solver_duration);
	
	fasp_dstr_free(&LU);
	
	return status;
}

/**
 * \fn int fasp_solver_dstr_krylov_blockgs(dSTRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 * \param *A:	pointer to the dSTRmat matrix
 * \param *b:	pointer to the dvector of right hand side
 * \param *x:	pointer to the dvector of dofs
 * \param *itparam: pointer to parameters for iterative solvers
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \data 10/10/2010
 */
int fasp_solver_dstr_krylov_blockgs (dSTRmat *A, 
																		 dvector *b, 
																		 dvector *x, 
																		 itsolver_param *itparam, 
																		 ivector *neigh, 
																		 ivector *order)
{
	//! Parameter for iterative method
	const int print_level = itparam->print_level;
	
	//! Information about matrices
	int ngrid=A->ngrid;
	
	//! return parameter
	int status = SUCCESS;
	
	//! local parameter
	clock_t solver_start, solver_end, setup_start, setup_end;
	double solver_duration, setup_duration;
	
	dvector *diaginv;
	ivector *pivot;
	
	// setup preconditioner
	setup_start=clock();
	
	diaginv = (dvector *)fasp_mem_calloc(ngrid, sizeof(dvector));
	pivot = (ivector *)fasp_mem_calloc(ngrid, sizeof(ivector));
	fasp_generate_diaginv_block(A, neigh, diaginv, pivot);
	
	precond_data_str precdata;
	precdata.A_str = A;
	precdata.diaginv = diaginv;
	precdata.pivot = pivot;
	precdata.order = order;
	precdata.neigh = neigh;
	
	precond prec; prec.data = &precdata; prec.fct  = fasp_precond_dstr_blockgs;
	
	setup_end=clock();
	
	if (print_level>0) {
		setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);	  
		printf("setup costs %f seconds.\n", setup_duration);
	}
	
	// solver part
	solver_start=clock();
	status=fasp_solver_dstr_itsolver(A,b,x,&prec,itparam);
	solver_end=clock();
	
	
	if (print_level>0) {
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("solve costs %f seconds.\n", solver_duration);
		printf("krylov_blockGS_str costs %f seconds.\n", solver_duration+setup_duration);
	}
	
	return status;
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
