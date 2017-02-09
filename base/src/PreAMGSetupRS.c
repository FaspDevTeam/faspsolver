/*! \file  PreAMGSetupRS.c
 *
 *  \brief Ruge-Stuben AMG: SETUP phase
 *
 *  \note  This file contains Level-4 (Pre) functions. It requires:
 *         AuxMemory.c, AuxMessage.c, AuxTiming.c, AuxVector.c, BlaILUSetupCSR.c,
 *         BlaSchwarzSetup.c, BlaSparseCSR.c, BlaSpmvCSR.c, PreAMGCoarsenRS.c,
 *         PreAMGInterp.c, and PreMGRecurAMLI.c
 *
 *  Reference: 
 *         Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller
 *         Appendix P475 A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *         Academic Press Inc., San Diego, CA, 2001.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_amg_setup_rs (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of Ruge and Stuben's classic AMG
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Chensong Zhang
 * \date   05/09/2010
 *
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle.
 * Modified by Xiaozhe Hu on 04/24/2013: aggressive coarsening.
 * Modified by Chensong Zhang on 09/23/2014: check coarse spaces.
 */
SHORT fasp_amg_setup_rs (AMG_data   *mgl,
                         AMG_param  *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,MIN_CDOF);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT      status = FASP_SUCCESS;
    INT        lvl = 0, max_lvls = param->max_levels;
    REAL       setup_start, setup_end;
    ILU_param  iluparam;
    SWZ_param  swzparam;
    iCSRmat    Scouple; // strong n-couplings

    // level info (fine: 0; coarse: 1)
    ivector    vertices = fasp_ivec_create(m);

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: n = %d, nnz = %d\n", mgl[0].A.row, mgl[0].A.nnz);
#endif

    fasp_gettime(&setup_start);

    // Make sure classical AMG will not call fasp_blas_dcsr_mxv_agg!
    param->tentative_smooth = 1.0;

    // If user want to use aggressive coarsening but did not specify number of
    // levels use aggressive coarsening, make sure apply aggressive coarsening
    // on the finest level only !!!
    if ( param->coarsening_type == COARSE_AC ) {
        param->aggressive_level = MAX(param->aggressive_level, 1);
    }

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)fasp_mem_calloc(amlideg+1,sizeof(REAL));
        fasp_amg_amli_coef(2.0, 0.5, amlideg, param->amli_coef);
    }

    // Initialize ILU parameters
    mgl->ILU_levels = param->ILU_levels;
    if ( param->ILU_levels > 0 ) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }

    // Initialize Schwarz parameters
    mgl->SWZ_levels = param->SWZ_levels;
    if ( param->SWZ_levels > 0 ) {
        swzparam.SWZ_mmsize = param->SWZ_mmsize;
        swzparam.SWZ_maxlvl = param->SWZ_maxlvl;
        swzparam.SWZ_type   = param->SWZ_type;
        swzparam.SWZ_blksolver = param->SWZ_blksolver;
    }

#if DIAGONAL_PREF
    // Reorder each row to keep the diagonal entries appear first !!!
    fasp_dcsr_diagpref(&mgl[0].A);
#endif

    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_lvls-1) ) {

#if DEBUG_MODE > 1
        printf("### DEBUG: level = %d, row = %d, nnz = %d\n",
               lvl, mgl[lvl].A.row, mgl[lvl].A.nnz);
#endif

        /*-- Setup ILU decomposition if needed --*/
        if ( lvl < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[lvl].A, &mgl[lvl].LU, &iluparam);
            if ( status < 0 ) {
                if ( prtlvl > PRINT_MIN ) {
                    printf("### WARNING: ILU setup on level-%d failed!\n", lvl);
                    printf("### WARNING: Disable ILU for level >= %d.\n", lvl);
                }
                param->ILU_levels = lvl;
            }
        }

        /*-- Setup Schwarz smoother if needed --*/
        if ( lvl < param->SWZ_levels ) {
            mgl[lvl].Schwarz.A = fasp_dcsr_sympart(&mgl[lvl].A);
            fasp_dcsr_shift(&(mgl[lvl].Schwarz.A), 1);
            status = fasp_swz_dcsr_setup(&mgl[lvl].Schwarz, &swzparam);
            if ( status < 0 ) {
                if ( prtlvl > PRINT_MIN ) {
                    printf("### WARNING: Schwarz on level-%d failed!\n", lvl);
                    printf("### WARNING: Disable Schwarz for level >= %d.\n", lvl);
                }
                param->SWZ_levels = lvl;
            }
        }

        /*-- Coarsening and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[lvl].A, &vertices, &mgl[lvl].P,
		                                &Scouple, param);

        // Check 1: Did coarsening step succeeded?
        if ( status < 0 ) {
            /*-- Clean up Scouple generated in coarsening --*/
            fasp_mem_free(Scouple.IA);
            fasp_mem_free(Scouple.JA);

            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Could not find any C-variables!\n");
                printf("### WARNING: RS coarsening on level-%d failed!\n", lvl);
            }
            status = FASP_SUCCESS; break;
        }

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) {
            /*-- Clean up Scouple generated in coarsening --*/
            fasp_mem_free(Scouple.IA);
            fasp_mem_free(Scouple.JA);
            break;
        }

        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * 10.0 ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening might be too aggressive!\n");
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
            }

            /*-- Clean up Scouple generated in coarsening --*/
            fasp_mem_free(Scouple.IA);
            fasp_mem_free(Scouple.JA);
            break;
        }

        /*-- Perform aggressive coarsening only up to the specified level --*/
        if ( mgl[lvl].P.col*1.5 > mgl[lvl].A.row ) param->coarsening_type = COARSE_RS;
        if ( lvl == param->aggressive_level ) param->coarsening_type = COARSE_RS;

        /*-- Store the C/F marker --*/
        {
            INT size = mgl[lvl].A.row;
            mgl[lvl].cfmark = fasp_ivec_create(size);
            memcpy(mgl[lvl].cfmark.val, vertices.val, size*sizeof(INT));
        }

        /*-- Form interpolation --*/
        fasp_amg_interp(&mgl[lvl].A, &vertices, &mgl[lvl].P, &Scouple, param);

        /*-- Form coarse level matrix: two RAP routines available! --*/
        fasp_dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        fasp_blas_dcsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);

        /*-- Clean up Scouple generated in coarsening --*/
        fasp_mem_free(Scouple.IA);
        fasp_mem_free(Scouple.JA);

        ++lvl;

#if DIAGONAL_PREF
        // reorder each row to make diagonal appear first
        fasp_dcsr_diagpref(&mgl[lvl].A);
#endif

        // Check 4: Is the coarse matrix too dense?
        if ( mgl[lvl].A.nnz / mgl[lvl].A.row > mgl[lvl].A.col * 0.2 ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarse matrix is too dense!\n");
                printf("### WARNING: m = n = %d, nnz = %d!\n",
                       mgl[lvl].A.col, mgl[lvl].A.nnz);
            }

            break;
        }
        
    } // end of the main while loop

    // Setup coarse level systems for direct solvers
    switch (csolver) {

#if WITH_MUMPS
        case SOLVER_MUMPS: {
            // Setup MUMPS direct solver on the coarsest level
            mgl[lvl].mumps.job = 1;
            fasp_solver_mumps_steps(&mgl[lvl].A, &mgl[lvl].b, &mgl[lvl].x, &mgl[lvl].mumps);
            break;
        }
#endif

#if WITH_UMFPACK
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            dCSRmat Ac_tran;
            Ac_tran = fasp_dcsr_create(mgl[lvl].A.row, mgl[lvl].A.col, mgl[lvl].A.nnz);
            fasp_dcsr_transz(&mgl[lvl].A, NULL, &Ac_tran);
            // It is equivalent to do transpose and then sort
            //     fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
            //     fasp_dcsr_sort(&Ac_tran);
            fasp_dcsr_cp(&Ac_tran, &mgl[lvl].A);
            fasp_dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = fasp_umfpack_factorize(&mgl[lvl].A, 0);
            break;
        }
#endif

#if WITH_PARDISO
        case SOLVER_PARDISO: {
             fasp_dcsr_sort(&mgl[lvl].A);
             fasp_pardiso_factorize(&mgl[lvl].A, &mgl[lvl].pdata, 0);
             break;
         }
#endif

        default:
            // Do nothing!
            break;
    }

    // setup total level number and current level
    mgl[0].num_levels = max_lvls = lvl+1;
    mgl[0].w          = fasp_dvec_create(m);

    for ( lvl = 1; lvl < max_lvls; ++lvl ) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_lvls;
        mgl[lvl].b          = fasp_dvec_create(mm);
        mgl[lvl].x          = fasp_dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!
        mgl[lvl].ILU_levels     = param->ILU_levels - lvl; // initialize ILU levels!
        mgl[lvl].SWZ_levels = param->SWZ_levels -lvl; // initialize Schwarz!

        // allocate work arrays for the solve phase
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = fasp_dvec_create(3*mm);
        else
            mgl[lvl].w = fasp_dvec_create(2*mm);
    }

    fasp_ivec_free(&vertices);

    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        print_amgcomplexity(mgl, prtlvl);
        print_cputime("Classical AMG setup", setup_end - setup_start);
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

