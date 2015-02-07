/*! \file amg_setup_sa.c
 *
 *  \brief Smoothed aggregation AMG: SETUP phase
 *
 *
 *  \note Setup A, P, PT and levels using the unsmoothed aggregation algorithm;
 *        Refer to P. Vanek, J. Madel and M. Brezina
 *        "Algebraic Multigrid on Unstructured Meshes", 1994
 *
 *
 */

#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"
#include "aggregation_csr.inl"
#include "aggregation_bsr.inl"

static SHORT amg_setup_smoothP_smoothR (AMG_data *, AMG_param *);
static SHORT amg_setup_smoothP_unsmoothR (AMG_data *, AMG_param *);
static SHORT amg_setup_smoothP_smoothR_bsr (AMG_data_bsr *mgl, AMG_param *param);

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
 * \return       FASP_SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * Modified by Chensong Zhang on 04/06/2010.
 * Modified by Chensong Zhang on 05/09/2010.
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle.
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
SHORT fasp_amg_setup_sa (AMG_data *mgl,
                         AMG_param *param)
{
    SHORT status  = FASP_SUCCESS;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
#if TRUE
    status = amg_setup_smoothP_smoothR(mgl, param);
#else // smoothed P, unsmoothed R
    status = amg_setup_smoothP_unsmoothR(mgl, param);
#endif
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/**
 * \fn INT fasp_amg_setup_sa_bsr (AMG_data_bsr *mgl, AMG_param *param)
 *
 * \brief Set up phase of smoothed aggregation AMG (BSR format)
 *
 * \param mgl    Pointer to AMG data: AMG_data_bsr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       FASP_SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
SHORT fasp_amg_setup_sa_bsr (AMG_data_bsr *mgl,
                             AMG_param *param)
{
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    SHORT status = amg_setup_smoothP_smoothR_bsr(mgl, param);
    
#if DEBUG_MODE
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
static SHORT amg_setup_smoothP_smoothR (AMG_data *mgl,
                                        AMG_param *param)
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
    Schwarz_param swzparam;
    
#if DEBUG_MODE
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
    mgl->Schwarz_levels = param->Schwarz_levels;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
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
        
#if DEBUG_MODE
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
        if ( lvl < param->Schwarz_levels ) {
            mgl[lvl].Schwarz.A = fasp_dcsr_sympat(&mgl[lvl].A);
            fasp_dcsr_shift(&(mgl[lvl].Schwarz.A), 1);
            fasp_Schwarz_setup(&mgl[lvl].Schwarz, &swzparam);
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
            fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
            fasp_dcsr_sort(&Ac_tran);
            fasp_dcsr_cp(&Ac_tran, &mgl[lvl].A);
            fasp_dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = fasp_umfpack_factorize(&mgl[lvl].A, 0);
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
        mgl[lvl].Schwarz_levels = param->Schwarz_levels -lvl; // initialize Schwarz!
        
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = fasp_dvec_create(3*mm);
        else
            mgl[lvl].w = fasp_dvec_create(2*mm);
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        print_amgcomplexity(mgl,prtlvl);
        print_cputime("Smoothed aggregation setup", setup_end - setup_start);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    fasp_mem_free(tentp);
    
#if DEBUG_MODE
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
static SHORT amg_setup_smoothP_unsmoothR (AMG_data *mgl,
                                          AMG_param *param)
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
    Schwarz_param swzparam;
    
#if DEBUG_MODE
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
    mgl->Schwarz_levels = param->Schwarz_levels;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_mmsize = param->Schwarz_mmsize;
        swzparam.Schwarz_maxlvl = param->Schwarz_maxlvl;
        swzparam.Schwarz_type   = param->Schwarz_type;
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
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
        if ( lvl < param->Schwarz_levels ) {
            mgl[lvl].Schwarz.A = fasp_dcsr_sympat(&mgl[lvl].A);
            fasp_dcsr_shift(&(mgl[lvl].Schwarz.A), 1);
            fasp_Schwarz_setup(&mgl[lvl].Schwarz, &swzparam);
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
            fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
            fasp_dcsr_sort(&Ac_tran);
            fasp_dcsr_cp(&Ac_tran, &mgl[lvl].A);
            fasp_dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = fasp_umfpack_factorize(&mgl[lvl].A, 5);
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
        mgl[lvl].Schwarz_levels = param->Schwarz_levels -lvl; // initialize Schwarz!

        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = fasp_dvec_create(3*mm);
        else
            mgl[lvl].w = fasp_dvec_create(2*mm);
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        print_amgcomplexity(mgl,prtlvl);
        print_cputime("Smoothed aggregation 1/2 setup", setup_end - setup_start);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    fasp_mem_free(tentp);
    fasp_mem_free(tentr);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}


/**
 * \fn static SHORT amg_setup_smoothP_smoothR_bsr (AMG_data_bsr *mgl,
 *                                                 AMG_param *param)
 *
 * \brief Set up phase of smoothed aggregation AMG, using smoothed P and smoothed A
 *        in BSR format
 *
 * \param mgl    Pointer to AMG data: AMG_data_bsr
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       FASP_SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 *
 */
static SHORT amg_setup_smoothP_smoothR_bsr (AMG_data_bsr *mgl,
                                            AMG_param *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT csolver    = param->coarse_solver;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.ROW;
    const INT   nb         = mgl[0].A.nb;
    
    ILU_param iluparam;
    SHORT     max_levels=param->max_levels;
    SHORT     i, lvl=0, status=FASP_SUCCESS;
    REAL      setup_start, setup_end;
    
    dCSRmat temp1, temp2;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n",
           mgl[0].A.ROW, mgl[0].A.COL, mgl[0].A.NNZ);
#endif
    
    fasp_gettime(&setup_start);
    
    /*-----------------------*/
    /*--local working array--*/
    /*-----------------------*/
    
    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels, sizeof(ivector));
    
    // each level stores the information of the number of aggregations
    INT *num_aggs = (INT *)fasp_mem_calloc(max_levels, sizeof(INT));
    
    // each level stores the information of the strongly coupled neighbourhood
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels, sizeof(dCSRmat));
    
    // each level stores the information of the tentative prolongations
    dBSRmat *tentp = (dBSRmat *)fasp_mem_calloc(max_levels,sizeof(dBSRmat));
    
    for ( i=0; i<max_levels; ++i ) num_aggs[i] = 0;
    
    /*-----------------------*/
    /*-- setup null spaces --*/
    /*-----------------------*/
    
    // null space for whole Jacobian
    //mgl[0].near_kernel_dim   = 1;
    //mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim, sizeof(REAL*));
    
    //for ( i=0; i < mgl->near_kernel_dim; ++i ) mgl[0].near_kernel_basis[i] = NULL;
    
    /*-----------------------*/
    /*-- setup ILU param   --*/
    /*-----------------------*/
    
    // initialize ILU parameters
    mgl->ILU_levels = param->ILU_levels;
    if ( param->ILU_levels > 0 ) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }
    
    /*----------------------------*/
    /*--- checking aggregation ---*/
    /*----------------------------*/
    
    if (param->aggregation_type == PAIRWISE)
        param->pair_number = MIN(param->pair_number, max_levels);
    
    // Main AMG setup loop
    while ( (mgl[lvl].A.ROW > min_cdof) && (lvl < max_levels-1) ) {
        
        printf("level = %d\n", lvl);
        
        /*-- setup ILU decomposition if necessary */
        if ( lvl < param->ILU_levels ) {
            status = fasp_ilu_dbsr_setup(&mgl[lvl].A, &mgl[lvl].LU, &iluparam);
            if ( status < 0 ) {
                if ( prtlvl > PRINT_MIN ) {
                    printf("### WARNING: ILU setup on level-%d failed!\n", lvl);
                    printf("### WARNING: Disable ILU for level >= %d.\n", lvl);
                }
                param->ILU_levels = lvl;
            }
        }
        
        /*-- get the diagonal inverse --*/
        mgl[lvl].diaginv = fasp_dbsr_getdiaginv(&mgl[lvl].A);
        
        /*-- Aggregation --*/
        //mgl[lvl].PP =  fasp_dbsr_getblk_dcsr(&mgl[lvl].A);
        mgl[lvl].PP = fasp_dbsr_Linfinity_dcsr(&mgl[lvl].A);
        
        switch ( param->aggregation_type ) {
                
            case VMB: // VMB aggregation
                // Same as default
                
            default: // only one aggregation is tested!!! --Chensong
                
                status = aggregation_vmb(&mgl[lvl].PP, &vertices[lvl], param, lvl+1,
                                         &Neighbor[lvl], &num_aggs[lvl]);
                
                /*-- Choose strength threshold adaptively --*/
                if ( num_aggs[lvl]*4 > mgl[lvl].PP.row )
                    param->strong_coupled /= 4;
                else if ( num_aggs[lvl]*1.25 < mgl[lvl].PP.row )
                    param->strong_coupled *= 1.5;
                
                break;
        }
        
        if ( status < 0 ) {
            // When error happens, force solver to use the current multigrid levels!
            if ( prtlvl > PRINT_MIN ) {
                printf("### WARNING: Aggregation on level-%d failed!\n", lvl);
            }
            status = FASP_SUCCESS; break;
        }
        
        /* -- Form Tentative prolongation --*/
        printf("before form tentative P\n");
        if (lvl == 0 && mgl[0].near_kernel_dim >0 ){
            form_tentative_p_bsr1(&vertices[lvl], &tentp[lvl], &mgl[0], lvl+1,
                                  num_aggs[lvl], mgl[0].near_kernel_dim,
                                  mgl[0].near_kernel_basis);
        }
        else{
            form_boolean_p_bsr(&vertices[lvl], &tentp[lvl], &mgl[0], lvl+1, num_aggs[lvl]);
        }
        
        /* -- Smoothing -- */
        printf("Smoothing P\n");
        smooth_agg_bsr(&mgl[lvl].A, &tentp[lvl], &mgl[lvl].P, param, lvl+1,
                       &Neighbor[lvl]);
        
        /*-- Form restriction --*/
        printf("Form Restriction\n");
        fasp_dbsr_trans(&mgl[lvl].P, &mgl[lvl].R);
        
        /*-- Form coarse level stiffness matrix --*/
        printf("RAP\n");
        fasp_blas_dbsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);
        
        /*-- Form extra near kernel space if needed --*/
        if (mgl[lvl].A_nk != NULL){
            
            mgl[lvl+1].A_nk = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            mgl[lvl+1].P_nk = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            mgl[lvl+1].R_nk = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            
            temp1 = fasp_format_dbsr_dcsr(&mgl[lvl].R);
            fasp_blas_dcsr_mxm(&temp1, mgl[lvl].P_nk, mgl[lvl+1].P_nk);
            fasp_dcsr_trans(mgl[lvl+1].P_nk, mgl[lvl+1].R_nk);
            temp2 = fasp_format_dbsr_dcsr(&mgl[lvl+1].A);
            fasp_blas_dcsr_rap(mgl[lvl+1].R_nk, &temp2, mgl[lvl+1].P_nk, mgl[lvl+1].A_nk);
            fasp_dcsr_free(&temp1);
            fasp_dcsr_free(&temp2);
            
        }
        
        printf("clean\n");
        fasp_dcsr_free(&Neighbor[lvl]);
        fasp_ivec_free(&vertices[lvl]);
        fasp_dbsr_free(&tentp[lvl]);
        
        ++lvl;
    }
    
    // Setup coarse level systems for direct solvers (BSR version)
    switch (csolver) {
            
#if WITH_MUMPS
        case SOLVER_MUMPS: {
            // Setup MUMPS direct solver on the coarsest level
            mgl[lvl].mumps.job = 1;
            mgl[lvl].Ac = fasp_format_dbsr_dcsr(&mgl[lvl].A);
            fasp_solver_mumps_steps(&mgl[lvl].Ac, &mgl[lvl].b, &mgl[lvl].x, &mgl[lvl].mumps);
            break;
        }
#endif
            
#if WITH_UMFPACK
        case SOLVER_UMFPACK: {
            // Need to sort the matrix A for UMFPACK to work
            mgl[lvl].Ac = fasp_format_dbsr_dcsr(&mgl[lvl].A);
            dCSRmat Ac_tran;
            fasp_dcsr_trans(&mgl[lvl].Ac, &Ac_tran);
            fasp_dcsr_sort(&Ac_tran);
            fasp_dcsr_cp(&Ac_tran, &mgl[lvl].Ac);
            fasp_dcsr_free(&Ac_tran);
            mgl[lvl].Numeric = fasp_umfpack_factorize(&mgl[lvl].Ac, 0);
            break;
        }
#endif
            
#if WITH_SuperLU
        case SOLVER_SUPERLU: {
            /* Setup SuperLU direct solver on the coarsest level */
            mgl[lvl].Ac = fasp_format_dbsr_dcsr(&mgl[lvl].A);
        }
#endif
            
        default:
            // Do nothing!
            break;
    }
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = lvl+1;
    mgl[0].w = fasp_dvec_create(3*m*nb);
    
    if (mgl[0].A_nk != NULL){
        
#if WITH_UMFPACK
        // Need to sort the matrix A_nk for UMFPACK
        fasp_dcsr_trans(mgl[0].A_nk, &temp1);
        fasp_dcsr_sort(&temp1);
        fasp_dcsr_cp(&temp1, mgl[0].A_nk);
        fasp_dcsr_free(&temp1);
        mgl[0].Numeric = fasp_umfpack_factorize(mgl[0].A_nk, 0);
#endif
        
    }
    
    for ( lvl = 1; lvl < max_levels; lvl++ ) {
        const INT mm = mgl[lvl].A.ROW*nb;
        mgl[lvl].num_levels = max_levels;
        mgl[lvl].b          = fasp_dvec_create(mm);
        mgl[lvl].x          = fasp_dvec_create(mm);
        mgl[lvl].w          = fasp_dvec_create(3*mm);
        mgl[lvl].ILU_levels = param->ILU_levels - lvl; // initialize ILU levels!
        
        if (mgl[lvl].A_nk != NULL){
            
#if WITH_UMFPACK
            // Need to sort the matrix A_nk for UMFPACK
            fasp_dcsr_trans(mgl[lvl].A_nk, &temp1);
            fasp_dcsr_sort(&temp1);
            fasp_dcsr_cp(&temp1, mgl[lvl].A_nk);
            fasp_dcsr_free(&temp1);
            mgl[lvl].Numeric = fasp_umfpack_factorize(mgl[lvl].A_nk, 0);
#endif
            
        }
        
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        print_amgcomplexity_bsr(mgl,prtlvl);
        print_cputime("Smoothed aggregation (BSR) setup", setup_end - setup_start);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggs);
    fasp_mem_free(Neighbor);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
