/*! \file SolFAMG.c
 *
 *  \brief Full AMG method as an iterative solver
 *
 *  \note This file contains Level-5 (Sol) functions. It requires
 *        AuxMessage.c, AuxTiming.c, AuxVector.c, BlaSparseCSR.c, PreAMGSetupRS.c,
 *        PreAMGSetupSA.c, PreAMGSetupUA.c, PreDataInit.c, and PreMGSolve.c
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_famg (dCSRmat *A, dvector *b, dvector *x, AMG_param *param)
 *
 * \brief Solve Ax=b by full AMG.
 *
 * \param A      Pointer to dCSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param x      Pointer to dvector: the unknowns
 * \param param  Pointer to AMG_param: AMG parameters
 *
 * \author Xiaozhe Hu
 * \date   02/27/2011
 *
 * Modified by Chensong Zhang on 01/10/2012
 * Modified by Chensong Zhang on 05/05/2013: Remove error handling for AMG setup
 */
void fasp_solver_famg (dCSRmat    *A,
                       dvector    *b,
                       dvector    *x,
                       AMG_param  *param)
{
    const SHORT   max_levels  = param->max_levels;
    const SHORT   prtlvl      = param->print_level;
    const SHORT   amg_type    = param->AMG_type;
    const INT     nnz = A->nnz, m = A->row, n = A->col;
    
    // local variables
    AMG_data *    mgl = fasp_amg_data_create(max_levels);
    REAL          FMG_start, FMG_end;
    
#if DEBUG_MODE > 0
    printf("###DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("###DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if ( prtlvl > PRINT_NONE ) fasp_gettime(&FMG_start);
    
    // Step 0: initialize mgl[0] with A, b and x
    mgl[0].A = fasp_dcsr_create(m,n,nnz);
    fasp_dcsr_cp(A,&mgl[0].A);
    
    mgl[0].b = fasp_dvec_create(n);
    fasp_dvec_cp(b,&mgl[0].b);
    
    mgl[0].x = fasp_dvec_create(n);
    fasp_dvec_cp(x,&mgl[0].x);
    
    // Step 1: AMG setup phase
    switch (amg_type) {
            
        case SA_AMG:
            // Smoothed Aggregation AMG setup phase
            if ( prtlvl > PRINT_NONE ) printf("\nCalling SA FAMG ...\n");
            fasp_amg_setup_sa(mgl, param); break;
            
        case UA_AMG:
            // Unsmoothed Aggregation AMG setup phase
            if ( prtlvl > PRINT_NONE ) printf("\nCalling UA FAMG ...\n");
            fasp_amg_setup_ua(mgl, param); break;
            
        default:
            // Classical AMG setup phase
            if ( prtlvl > PRINT_NONE ) printf("\nCalling FAMG ...\n");
            fasp_amg_setup_rs(mgl, param);
            
    }
    
    // Step 2: FAMG solve phase
    fasp_famg_solve(mgl, param);
    
    // Step 3: Save solution vector and return
    fasp_dvec_cp(&mgl[0].x, x);
    
    // clean-up memory
    fasp_amg_data_free(mgl, param);
    
    // print out CPU time if needed
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&FMG_end);
        print_cputime("FAMG totally", FMG_end - FMG_start);
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
