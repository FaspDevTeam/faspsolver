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
 * \fn void fasp_solver_famg(dCSRmat *A, dvector *b, dvector *x, AMG_param *param)
 *
 * \brief Solve Ax=b by full AMG.
 *
 * \param A      Pointer to the coefficient matrix
 * \param b      Pointer to the dvector of right hand side
 * \param x      Pointer to the dvector of dofs
 * \param param  Pointer to AMG parameters
 *
 * \author Xiaozhe Hu
 * \date   02/27/2011
 *
 * Modified by Chensong Zhang on 01/10/2012
 * Modified by Chensong Zhang on 05/05/2013: Remove error handling for AMG setup
 */
void fasp_solver_famg (dCSRmat *A,
                       dvector *b,
                       dvector *x,
                       AMG_param *param)
{
    const SHORT   max_levels  = param->max_levels;
    const SHORT   print_level = param->print_level;
    const SHORT   amg_type    = param->AMG_type;
    const INT     nnz = A->nnz, m = A->row, n = A->col;
    
    // local variables
    REAL          FMG_start, FMG_end;
    SHORT         status = SUCCESS;
    AMG_data *    mgl = fasp_amg_data_create(max_levels);
    
#if DEBUG_MODE
    printf("###DEBUG: fasp_solver_famg ...... [Start]\n");
    printf("###DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if ( print_level > PRINT_NONE ) {
        fasp_gettime(&FMG_start);
    }
    
    // Step 0: initialize mgl[0] with A, b and x
    mgl[0].A = fasp_dcsr_create(m,n,nnz);
    fasp_dcsr_cp(A,&mgl[0].A);
    
    mgl[0].b = fasp_dvec_create(n);
    fasp_dvec_cp(b,&mgl[0].b);
    
    mgl[0].x = fasp_dvec_create(n);
    fasp_dvec_cp(x,&mgl[0].x);
    
    // Step 1: AMG setup phase
    switch (amg_type) {
            
            // Smoothed Aggregation AMG setup phase
        case SA_AMG:
            if ( print_level > PRINT_NONE ) printf("\nCalling SA FAMG solver ...\n");
            fasp_amg_setup_sa(mgl, param);
            break;
            
            // Unsmoothed Aggregation AMG setup phase
        case UA_AMG:
            if ( print_level > PRINT_NONE ) printf("\nCalling UA FAMG solver ...\n");
            fasp_amg_setup_ua(mgl, param);
            break;
            
            // Classical AMG setup phase
        default:
            if ( print_level > PRINT_NONE ) printf("\nCalling classical FAMG solver ...\n");
#ifdef _OPENMP // omp version RS coarsening
            fasp_amg_setup_rs_omp(mgl, param);
#else
            fasp_amg_setup_rs(mgl, param);
#endif
            break;
            
    }
    
    // Step 2: FAMG solve phase
    fasp_famg_solve(mgl, param);
    
    // Step 3: Save solution vector and return
    fasp_dvec_cp(&mgl[0].x, x);
    
    // print out CPU time if needed
    if (print_level>PRINT_NONE) {
        fasp_gettime(&FMG_end);
        printf("FMG totally costs %f seconds.\n", FMG_end - FMG_start);
    }
    
FINISHED:
    fasp_amg_data_free(mgl);    // clean-up memory
    
    fasp_chkerr(status, "fasp_solver_famg");
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_famg ...... [Finish]\n");
#endif
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
