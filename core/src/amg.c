/*! \file amg.c 
 *  \brief AMG method as an iterative solver (main file)
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_amg (dCSRmat *A, dvector *b, dvector *x, 
 *                           AMG_param *param)
 *
 * \brief Solve Ax=b by Algebaric MultiGrid Method.
 *
 * \param A      pointer to the coefficient matrix
 * \param b      pointer to the dvector of right hand side
 * \param x      pointer to the dvector of dofs
 * \param param  pointer to AMG parameters
 *
 * \note
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Chensong Zhang
 * \date 04/06/2010
 *
 * \note Modified by Chensong Zhang on 01/10/2012
 */
void fasp_solver_amg (dCSRmat *A, 
                      dvector *b, 
                      dvector *x, 
                      AMG_param *param)
{
    const SHORT   max_levels=param->max_levels;
    const SHORT   print_level=param->print_level;
    const SHORT   amg_type=param->AMG_type;
    const SHORT   cycle_type=param->cycle_type;    
    const INT     nnz=A->nnz, m=A->row, n=A->col;
	
    // local variables
    clock_t       AMG_start, AMG_end;
    REAL          AMG_duration;
    SHORT         status = SUCCESS;
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amg ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
    AMG_start = clock();
	
    // initialize mgl[0] with A, b, x	
    AMG_data *mgl=fasp_amg_data_create(max_levels);
    mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
    mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp(b,&mgl[0].b); 	
    mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp(x,&mgl[0].x); 
	
    // AMG setup phase
    switch (amg_type) {
        case CLASSIC_AMG: // Classical AMG setup phase
            if ( (status=fasp_amg_setup_rs(mgl, param))<0 ) goto FINISHED;
            break;
        case SA_AMG: // Smoothed Aggregation AMG setup phase
            if ( (status=fasp_amg_setup_sa(mgl, param))<0 ) goto FINISHED;
            break;
        case UA_AMG: // Unsmoothed Aggregation AMG setup phase
            if ( (status=fasp_amg_setup_ua(mgl, param))<0 ) goto FINISHED;
            break;
        default: // Unknown setup type
            status=ERROR_SOLVER_TYPE; goto FINISHED;
    }
	
    // AMG solve phase
    switch (cycle_type)
    {
        case AMLI_CYCLE: // call AMLI-cycle
            if ( (status=fasp_amg_solve_amli(mgl, param)) < 0 ) goto FINISHED;
            break;
        case NL_AMLI_CYCLE: // call Nonlinear AMLI-cycle
            if ( (status=fasp_amg_solve_nl_amli(mgl, param)) < 0 ) goto FINISHED;
            break;
        default: // call classical V,W-cycles
            if ( (status=fasp_amg_solve(mgl, param)) < 0 ) goto FINISHED;
            break;
    }
	
    // save solution vector
    fasp_dvec_cp(&mgl[0].x, x); 
	
    // print out CPU time when needed
    if (print_level>PRINT_NONE) {
        AMG_end = clock();		
        AMG_duration = (double)(AMG_end - AMG_start)/(double)(CLOCKS_PER_SEC);		
        printf("AMG totally costs %f seconds.\n", AMG_duration);
    }	
	
FINISHED:	
    fasp_amg_data_free(mgl); // clean-up memory	
    
    fasp_chkerr(status, "fasp_solver_amg");
    
#if DEBUG_MODE
    printf("fasp_solver_amg ...... [Finish]\n");
#endif
	
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
