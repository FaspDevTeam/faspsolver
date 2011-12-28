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
 * \fn INT fasp_solver_amg (dCSRmat *A, dvector *b, dvector *x, 
 *                          AMG_param *param)
 *
 * \brief Solve Ax=b by Algebaric MultiGrid Method.
 *
 * \param * A      pointer to the coefficient matrix
 * \param * b      pointer to the dvector of right hand side
 * \param * x      pointer to the dvector of dofs
 * \param * param  pointer to AMG parameters
 *
 * \note
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Chensong Zhang
 * \date 04/06/2010
 */
INT fasp_solver_amg (dCSRmat *A, 
                     dvector *b, 
                     dvector *x, 
                     AMG_param *param)
{
    const INT nnz=A->nnz, m=A->row, n=A->col;
    const INT max_levels=param->max_levels;
    const INT print_level=param->print_level;
    const INT amg_type=param->AMG_type;
    const INT cycle_type=param->cycle_type;
	
    // local variables
    clock_t AMG_start, AMG_end;
    double AMG_duration;
    INT status = SUCCESS;
	
#if DEBUG_MODE
    printf("### DEBUG: amg ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
    AMG_start = clock();
	
    // initialize mgl[0] with A, b, x	
    AMG_data *mgl=fasp_amg_data_create(max_levels);
    mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
    mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp(b,&mgl[0].b); 	
    mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp(x,&mgl[0].x); 
	
    // AMG setup	
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
        default:
            printf("### ERROR: Unknown AMG type %d!\n", amg_type); goto FINISHED;
    }
	
    // AMG solve
    switch (cycle_type)
    {
        case AMLI_CYCLE: // call AMLI-cycle
            if ( (status=fasp_amg_solve_amli(mgl, param)) < 0 ) goto FINISHED;
            break;
        default:
            if ( (status=fasp_amg_solve(mgl, param)) < 0 ) goto FINISHED;
            break;
    }
	
    fasp_dvec_cp(&mgl[0].x, x); // save solution vector
	
    if (print_level>PRINT_NONE) {
        AMG_end = clock();		
        AMG_duration = (double)(AMG_end - AMG_start)/(double)(CLOCKS_PER_SEC);		
        printf("AMG totally costs %f seconds.\n", AMG_duration);
    }	
	
FINISHED:	
    fasp_amg_data_free(mgl); // clean-up memory	
    
    if (status<0) {
        printf("### ERROR: AMG solver quit unexpectedly: %d!\n", status);
    }
    
#if DEBUG_MODE
    printf("amg ...... [Finish]\n");
#endif
	
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
