/*! \file itsolver_Stokes.c
 *  \brief Iterative solvers for Stokes-type matrices (main file)
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
 * \fn INT fasp_solver_bdcsr_krylov_stokes (block_dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam, precond_Stokes_param *param, precond_Stokes_data *pcdata)
 * \brief Solve Ax=b by standard Krylov methods 
 *
 * \param A	       pointer to the block_dCSRmat matrix
 * \param b	       pointer to the dvector of right hand side
 * \param x	       pointer to the dvector of dofs
 * \param itparam  pointer to parameters for iterative solvers
 * \param param    parameters to Stokes problems
 * \param pcdata   pionter to preconditioner data for Stokes
 *
 * \return         number of iterations
 *
 * \author Chensong Zhang
 * \date   11/25/2010
 */
INT fasp_solver_bdcsr_krylov_stokes (block_dCSRmat *Mat, 
                                     dvector *b, 
                                     dvector *x, 
                                     itsolver_param *itparam,
                                     precond_Stokes_param *param, 
                                     precond_Stokes_data *pcdata)
{
	// parameters
	const SHORT print_level  = itparam->print_level;
	const SHORT precond_type = itparam->precond_type;
	
	// Stokes matrix 
	dCSRmat *A = Mat->blocks[0];
	dCSRmat *B = Mat->blocks[1];
	const INT n = A->row, nnzA = A->nnz, m = B->row;	
	
	// preconditioner data
	dCSRmat *M = pcdata->M;
	precond pc;
	AMG_param amgparam;
	dvector diag_M;	
	
	// local variable
	clock_t solver_start, solver_end, setup_start, setup_end;
	REAL solver_duration, setup_duration;
	INT status=SUCCESS;
	
	// initialize preconditioner 
	pc.data = &pcdata; 
	switch (precond_type) {
		case 1:
			pc.fct = fasp_precond_stokes_bdiag;
			break;
		default:
			printf("### ERROR: Unknown preconditioner type!\n");
			exit(ERROR_SOLVER_PRECTYPE);
	}
	
	// AMG parameters
	amgparam.print_level = param->print_level;
	amgparam.max_levels = param->max_levels;
	amgparam.AMG_type = param->AMG_type;
	
	//------ setup phase ------//
	setup_start = clock();
	
	pcdata->colA = n;
	pcdata->colB = m;
	pcdata->col  = n+m;	
	
	// setup work array space
	pcdata->w = (REAL *)fasp_mem_calloc(pcdata->col,sizeof(REAL));
	
	// initialize AMG for A
	AMG_data *mgl=fasp_amg_data_create(amgparam.max_levels);
    mgl[0].A=fasp_dcsr_create(n,n,nnzA); fasp_dcsr_cp(A,&mgl[0].A);
	mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n);
	
	// setup AMG
	switch (amgparam.AMG_type) {
		case CLASSIC_AMG:
			fasp_amg_setup_rs(mgl, &amgparam);
			break;
		case SA_AMG:
			fasp_amg_setup_sa(mgl, &amgparam);
			break;
		default:
			printf("### ERROR Wrong AMG type %d!\n",amgparam.AMG_type);
			exit(ERROR_INPUT_PAR);
	}	
	pcdata->max_levels = mgl[0].num_levels;
	pcdata->mgl_data = mgl;
	
	// setup diagonal for M
	fasp_dcsr_getdiag(0,M,&diag_M);	
	pcdata->diag_M = &diag_M;
	
	setup_end = clock();
	
	if (print_level>PRINT_NONE) {
		setup_duration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Setup costs %f.\n", setup_duration);
	}
	
	//------ solver phase ------// 
	solver_start=clock();
	status=fasp_solver_bdcsr_itsolver(Mat,b,x,&pc,itparam);
	solver_end=clock();
	
	if (print_level>PRINT_NONE) {
		solver_duration = (REAL)(solver_end - solver_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Solver costs %f seconds.\n", solver_duration);	
		printf("Total costs %f seconds.\n", setup_duration + solver_duration);
	}
	
	// clean up memory
	if (mgl) fasp_amg_data_free(mgl);
	fasp_mem_free(pcdata->w);
	
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
