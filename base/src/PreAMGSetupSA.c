/*! \file  PreAMGSetupSA.c
 *
 *  \brief Smoothed aggregation AMG: SETUP phase
 *
 *  \note  This file contains Level-4 (Pre) functions. It requires:
 *         AuxArray.c, AuxMemory.c, AuxMessage.c, AuxThreads.c, AuxTiming.c,
 *         AuxVector.c, BlaILUSetupCSR.c, BlaSchwarzSetup.c, BlaSparseCSR.c, 
 *         BlaSpmvCSR.c, and PreMGRecurAMLI.c
 *
 *  \note  Setup A, P, PT and levels using the unsmoothed aggregation algorithm
 *
 *  Reference: 
 *         P. Vanek, J. Madel and M. Brezina
 *         Algebraic Multigrid on Unstructured Meshes, 1994
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/

#include "PreAMGAggregationCSR.inl"
#include "PreAMGAggregation.inl"

static SHORT amg_setup_smoothP_smoothR (AMG_data *, AMG_param *);
static SHORT amg_setup_smoothP_unsmoothR (AMG_data *, AMG_param *);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_amg_setup_sa (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of smoothed aggregation AMG
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       FASP_SUCCESS if successed; otherwise, error information.
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle.
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
SHORT fasp_amg_setup_sa (AMG_data   *mgl,
                         AMG_param  *param)
{
    SHORT status  = FASP_SUCCESS;

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif

#if TRUE
    status = amg_setup_smoothP_smoothR(mgl, param);
#else // smoothed P, unsmoothed R
    status = amg_setup_smoothP_unsmoothR(mgl, param);
#endif

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static SHORT amg_setup_smoothP_smoothR (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of smoothed aggregation AMG, using smoothed P and smoothed R
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 * Modified by Chensong Zhang on 07/26/2014: handle coarsening errors.
 * Modified by Chensong Zhang on 09/23/2014: check coarse spaces.
 */
static SHORT amg_setup_smoothP_smoothR (AMG_data   *mgl,
                                        AMG_param  *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT       max_levels = param->max_levels, lvl = 0, status = FASP_SUCCESS;
    INT         i, j;
    REAL        setup_start, setup_end;
    ILU_param   iluparam;
    SWZ_param   swzparam;

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    fasp_gettime(&setup_start);

    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));

    // each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));

    // each level stores the information of the strongly coupled neighbourhood
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));

    // each level stores the information of the tentative prolongations
    dCSRmat *tentp = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));

    // Initialize level information
    for ( i = 0; i < max_levels; ++i ) num_aggregations[i] = 0;

    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(REAL*));

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        mgl[0].near_kernel_basis[i] = (REAL *)fasp_mem_calloc(m,sizeof(REAL));
        for ( j = 0; j < m; ++j ) mgl[0].near_kernel_basis[i][j] = 1.0;
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

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)fasp_mem_calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }

#if DIAGONAL_PREF
    fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif

    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    if ( param->aggregation_type == PAIRWISE )
        param->pair_number = MIN(param->pair_number, max_levels);

    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {

#if DEBUG_MODE > 2
        printf("### DEBUG: level = %d, row = %d, nnz = %d\n",
               lvl, mgl[lvl].A.row, mgl[lvl].A.nnz);
#endif

        /*-- setup ILU decomposition if necessary */
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

        /* -- setup Schwarz smoother if necessary */
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

        /*-- Aggregation --*/
        status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param, lvl+1,
                                 &Neighbor[lvl], &num_aggregations[lvl]);

        // Check 1: Did coarsening step succeed?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Forming aggregates on level-%d failed!\n", lvl);
            }
            status = FASP_SUCCESS; 
            fasp_ivec_free(&vertices[lvl]);
            fasp_dcsr_free(&Neighbor[lvl]);
            break;
        }

        /* -- Form Tentative prolongation --*/
        form_tentative_p(&vertices[lvl], &tentp[lvl], mgl[0].near_kernel_basis,
                         lvl+1, num_aggregations[lvl]);

        /* -- Form smoothed prolongation -- */
        smooth_agg(&mgl[lvl].A, &tentp[lvl], &mgl[lvl].P, param, lvl+1,
                   &Neighbor[lvl]);

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) {
            fasp_ivec_free(&vertices[lvl]);
            fasp_dcsr_free(&Neighbor[lvl]);
            fasp_dcsr_free(&tentp[lvl]);
            break;
        }

        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * MAX_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening might be too aggressive!\n");
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
            }
            fasp_ivec_free(&vertices[lvl]);
            fasp_dcsr_free(&Neighbor[lvl]);
            fasp_dcsr_free(&tentp[lvl]);
            break;
        }

        /*-- Form restriction --*/
        fasp_dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);

        /*-- Form coarse level stiffness matrix --*/
        fasp_blas_dcsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);

        fasp_dcsr_free(&Neighbor[lvl]);
        fasp_dcsr_free(&tentp[lvl]);
        fasp_ivec_free(&vertices[lvl]);

        ++lvl;

#if DIAGONAL_PREF
        // reorder each row to make diagonal appear first
        fasp_dcsr_diagpref(&mgl[lvl].A);
#endif

        // Check 4: Is this coarsening ratio too small?
        if ( (REAL)mgl[lvl].P.col > mgl[lvl].P.row * MIN_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening rate is too small!\n");
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
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
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w          = fasp_dvec_create(m);

    for ( lvl = 1; lvl < max_levels; ++lvl) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = fasp_dvec_create(mm);
        mgl[lvl].x          = fasp_dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!
        mgl[lvl].ILU_levels     = param->ILU_levels - lvl; // initialize ILU levels!
        mgl[lvl].SWZ_levels = param->SWZ_levels -lvl; // initialize Schwarz!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = fasp_dvec_create(3*mm);
        else
            mgl[lvl].w = fasp_dvec_create(2*mm);
    }

    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        fasp_amgcomplexity(mgl,prtlvl);
        fasp_cputime("Smoothed aggregation setup", setup_end - setup_start);
    }

    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    fasp_mem_free(tentp);

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn static SHORT amg_setup_smoothP_unsmoothR (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of plain aggregation AMG, using smoothed P and unsmoothed R
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 * Modified by Chensong Zhang on 07/26/2014: handle coarsening errors.
 * Modified by Chensong Zhang on 09/23/2014: check coarse spaces.
 */
static SHORT amg_setup_smoothP_unsmoothR (AMG_data   *mgl,
                                          AMG_param  *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;

    // local variables
    SHORT       max_levels = param->max_levels, lvl = 0, status = FASP_SUCCESS;
    INT         i, j;
    REAL        setup_start, setup_end;
    ILU_param   iluparam;
    SWZ_param swzparam;

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    fasp_gettime(&setup_start);

    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));

    // each level stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));

    // each level stores the information of the strongly coupled neighbourhood
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));

    // each level stores the information of the tentative prolongations
    dCSRmat *tentp = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    dCSRmat *tentr = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));

    for ( i = 0; i < max_levels; ++i ) num_aggregations[i] = 0;

    mgl[0].near_kernel_dim   = 1;

    mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(REAL*));

    for ( i = 0; i < mgl->near_kernel_dim; ++i ) {
        mgl[0].near_kernel_basis[i] = (REAL *)fasp_mem_calloc(m,sizeof(REAL));
        for ( j = 0; j < m; ++j ) mgl[0].near_kernel_basis[i][j] = 1.0;
    }

    // Initialize ILU parameters
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

    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)fasp_mem_calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }

    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_levels-1) ) {

        /*-- setup ILU decomposition if necessary */
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

        /* -- setup Schwarz smoother if necessary */
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

        /*-- Aggregation --*/
        status = aggregation_vmb(&mgl[lvl].A, &vertices[lvl], param, lvl+1,
                                 &Neighbor[lvl], &num_aggregations[lvl]);

        // Check 1: Did coarsening step succeeded?
        if ( status < 0 ) {
            // When error happens, stop at the current multigrid level!
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Forming aggregates on level-%d failed!\n", lvl);
            }
            status = FASP_SUCCESS; break;
        }

        /* -- Form Tentative prolongation --*/
        form_tentative_p(&vertices[lvl], &tentp[lvl], mgl[0].near_kernel_basis,
                         lvl+1, num_aggregations[lvl]);

        /* -- Form smoothed prolongation -- */
        smooth_agg(&mgl[lvl].A, &tentp[lvl], &mgl[lvl].P, param, lvl+1,
                   &Neighbor[lvl]);

        // Check 2: Is coarse sparse too small?
        if ( mgl[lvl].P.col < MIN_CDOF ) break;

        // Check 3: Does this coarsening step too aggressive?
        if ( mgl[lvl].P.row > mgl[lvl].P.col * MAX_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening might be too aggressive!\n");
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
            }
            break;
        }

        // Check 4: Is this coarsening ratio too small?
        if ( (REAL)mgl[lvl].P.col > mgl[lvl].P.row * MIN_CRATE ) {
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Coarsening rate is too small!\n");
                printf("### WARNING: Fine level = %d, coarse level = %d. Discard!\n",
                       mgl[lvl].P.row, mgl[lvl].P.col);
            }
            break;
        }

        /*-- Form restriction --*/
        fasp_dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);
        fasp_dcsr_trans(&tentp[lvl], &tentr[lvl]);

        /*-- Form coarse level stiffness matrix --*/
        fasp_blas_dcsr_rap_agg(&tentr[lvl], &mgl[lvl].A, &tentp[lvl], &mgl[lvl+1].A);

        fasp_dcsr_free(&Neighbor[lvl]);
        fasp_dcsr_free(&tentp[lvl]);
        fasp_ivec_free(&vertices[lvl]);

        ++lvl;
    }

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
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w          = fasp_dvec_create(m);

    for ( lvl = 1; lvl < max_levels; ++lvl) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = fasp_dvec_create(mm);
        mgl[lvl].x          = fasp_dvec_create(mm);

        mgl[lvl].cycle_type     = cycle_type; // initialize cycle type!
        mgl[lvl].ILU_levels     = param->ILU_levels - lvl; // initialize ILU levels!
        mgl[lvl].SWZ_levels = param->SWZ_levels -lvl; // initialize Schwarz!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = fasp_dvec_create(3*mm);
        else
            mgl[lvl].w = fasp_dvec_create(2*mm);
    }

    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        fasp_amgcomplexity(mgl,prtlvl);
        fasp_cputime("Smoothed aggregation 1/2 setup", setup_end - setup_start);
    }

    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    fasp_mem_free(tentp);
    fasp_mem_free(tentr);

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
