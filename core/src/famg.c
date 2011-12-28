/*! \file famg.c
 *  \brief full AMG method as an iterative solver (main file)
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_solver_famg(dCSRmat *A, dvector *b, dvector *x, AMG_param *param)
 * \brief Solve Ax=b by full AMG.
 *
 * \param *A      pointer to the coefficient matrix
 * \param *b      pointer to the dvector of right hand side
 * \param *x      pointer to the dvector of dofs
 * \param *param  pointer to AMG parameters
 *
 * \author Xiaozhe Hu
 * \date 02/27/2011
 */
int fasp_solver_famg (dCSRmat *A, 
                      dvector *b, 
											dvector *x, 
											AMG_param *param)
{
	const int nnz=A->nnz, m=A->row, n=A->col;
	const int max_levels=param->max_levels;
	const int print_level=param->print_level;
	const int amg_type=param->AMG_type;
	
	clock_t AMG_start, AMG_end;
	double AMG_duration;
	int status = SUCCESS;
	
#if DEBUG_MODE
	printf("famg ...... [Start]\n");
	printf("famg: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	AMG_start=clock();
	
	// initialize A, b, x for mgl[0]	
	AMG_data *mgl=fasp_amg_data_create(max_levels);	
	mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
	mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp(b,&mgl[0].b); 	
	mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp(x,&mgl[0].x); 
	
	// AMG setup	
	switch (amg_type) {
		case CLASSIC_AMG:
			if ( (status=fasp_amg_setup_rs(mgl, param))<0 ) goto FINISHED;
			break;
		case SA_AMG:
			if ( (status=fasp_amg_setup_sa(mgl, param))<0 ) goto FINISHED;
			break;
        case UA_AMG:
			if ( (status=fasp_amg_setup_ua(mgl, param))<0 ) goto FINISHED;
			break;
		default:
			printf("Error: Wrong AMG type %d!\n", amg_type); goto FINISHED;
	}
#if 0	
	// full AMG solve
	if ( (status=famg_solve(mgl, param)) < 0 ) goto FINISHED;
#endif	
	fasp_dvec_cp(&mgl[0].x, x); // save solution vector
	
	if (print_level>0) {
		AMG_end=clock();		
		AMG_duration = (double)(AMG_end - AMG_start)/(double)(CLOCKS_PER_SEC);		
		printf("full AMG totally costs %f seconds.\n", AMG_duration);
	}	
	
FINISHED:	
	fasp_amg_data_free(mgl);	// clean-up memory	
	if (status<0) printf("Error: full AMG solver quit unexpectedly!\n");
	
#if DEBUG_MODE
	printf("famg ...... [Finish]\n");
#endif
	
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
