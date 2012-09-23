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
 * \param A      Pointer to the coefficient matrix
 * \param b      Pointer to the dvector of right hand side
 * \param x      Pointer to the dvector of dofs
 * \param param  Pointer to AMG parameters
 *
 * \author Chensong Zhang
 * \date   04/06/2010
 *
 * \note Refter to "Multigrid"
 *       by U. Trottenberg, C. W. Oosterlee and A. Schuller
 *       Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *       Academic Press Inc., San Diego, CA, 2001.
 *
 * Modified by Chensong Zhang on 01/10/2012
 */

void fasp_solver_amg (dCSRmat *A,
                      dvector *b,
                      dvector *x,
                      AMG_param *param)
{
    const SHORT   max_levels  = param->max_levels;
    const SHORT   print_level = param->print_level;
    const SHORT   amg_type    = param->AMG_type;
    const SHORT   cycle_type  = param->cycle_type;
    const INT     nnz = A->nnz, m = A->row, n = A->col;
    
    // local variables
    REAL          AMG_start, AMG_end;
    REAL          AMG_duration = 0.0;
    SHORT         status = SUCCESS;
    AMG_data *    mgl = fasp_amg_data_create(max_levels);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amg ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if ( print_level > PRINT_NONE ) {
#ifdef _OPENMP
        AMG_start = omp_get_wtime();
#else
        AMG_start = clock();
#endif
    }
    
    // Step 0: initialize mgl[0] with A, b and x
    mgl[0].A = fasp_dcsr_create(m, n, nnz);
    fasp_dcsr_cp(A, &mgl[0].A);
    
    mgl[0].b = fasp_dvec_create(n);
    fasp_dvec_cp(b, &mgl[0].b);
    
    mgl[0].x = fasp_dvec_create(n);
    fasp_dvec_cp(x, &mgl[0].x);

    
    // Step 1: AMG setup phase
    switch (amg_type) {
            
        // Smoothed Aggregation AMG setup phase
        case SA_AMG:
            if ( (status=fasp_amg_setup_sa(mgl, param)) < 0 ) goto FINISHED;
            break;

        // Unsmoothed Aggregation AMG setup phase
        case UA_AMG:
            if ( (status=fasp_amg_setup_ua(mgl, param)) < 0 ) goto FINISHED;
            break;

        // Classical AMG setup phase
        default:
#ifdef _OPENMP // omp version RS coarsening
            if ( (status=fasp_amg_setup_rs_omp(mgl, param)) < 0 ) goto FINISHED;
#else
            if ( (status=fasp_amg_setup_rs(mgl, param)) < 0 ) goto FINISHED;
#endif
            break;
            
    }
    
    // Step 2: AMG solve phase
    switch (cycle_type) {

        // call AMLI-cycle
        case AMLI_CYCLE:
            if ( (status=fasp_amg_solve_amli(mgl, param)) < 0 ) goto FINISHED;
            break;

        // call Nonlinear AMLI-cycle
        case NL_AMLI_CYCLE:
            if ( (status=fasp_amg_solve_nl_amli(mgl, param)) < 0 ) goto FINISHED;
            break;

        // call classical V,W-cycles (determined by param)
        default:
            if ( (status=fasp_amg_solve(mgl, param)) < 0 ) goto FINISHED;
            break;
            
    }
    
    // Step 3: Save solution vector and return
    fasp_dvec_cp(&mgl[0].x, x);
    
    // print out CPU time if needed
    if ( print_level > PRINT_NONE ) {
#ifdef _OPENMP
        AMG_end = omp_get_wtime();
        AMG_duration = AMG_end - AMG_start;
#else
        AMG_end = clock();
        AMG_duration = (REAL)(AMG_end - AMG_start)/(REAL)(CLOCKS_PER_SEC);
#endif
        print_cputime("AMG totally", AMG_duration);
    }
    
FINISHED:
    fasp_amg_data_free(mgl); // clean-up memory
    fasp_chkerr(status, "fasp_solver_amg");
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amg ...... [Finish]\n");
#endif
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
