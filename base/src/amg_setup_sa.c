/*! \file amg_setup_sa.c
 *
 *  \brief Smoothed aggregation AMG: SETUP phase
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
static SHORT amg_setup_smoothP_smoothR_bsr (AMG_data_bsr *mgl,
                                            AMG_param *param);

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
 */
static SHORT amg_setup_smoothP_smoothR (AMG_data *mgl,
                                        AMG_param *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;
    
    // local variables
    SHORT       max_levels = param->max_levels, level=0, status=FASP_SUCCESS;
    INT         i, j;
    REAL        setup_start, setup_end;
    ILU_param   iluparam;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    fasp_gettime(&setup_start);
    
    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    // each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));
    
    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    // each level stores the information of the tentative prolongations
    dCSRmat *tentp = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    // Initialzie level information
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
    while ( (mgl[level].A.row > min_cdof) && (level < max_levels-1) ) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %d, row = %d, nnz = %d\n",
               level,mgl[level].A.row,mgl[level].A.nnz);
#endif
        
        /*-- setup ILU decomposition if necessary */
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A, &mgl[level].LU, &iluparam);
            if ( status < 0 ) {
                param->ILU_levels = level;
                printf("### WARNING: ILU setup on level%d failed!\n", level);
            }
        }
        
        /* -- setup Schwarz smoother if necessary */
		if ( level < param->schwarz_levels ) {
            const INT smmsize = param->schwarz_mmsize;
            const INT smaxlvl = param->schwarz_maxlvl;
            const INT schtype = param->schwarz_type;
            
            mgl->schwarz_levels  = param->schwarz_levels;
            mgl[level].schwarz.A = fasp_dcsr_sympat(&mgl[level].A);
            fasp_dcsr_shift(&(mgl[level].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[level].schwarz, smmsize, smaxlvl, schtype);
		}
        
        /*-- Aggregation --*/
        status = aggregation_vmb(&mgl[level].A, &vertices[level], param, level+1,
                                 &Neighbor[level], &num_aggregations[level]);
		
        if ( status < 0 ) {
            // When error happens, force solver to use the current multigrid levels!
            printf("### WARNING: Coarsening on level%d is not successful!\n", level);
            status = FASP_SUCCESS; break;
        }
        
        /* -- Form Tentative prolongation --*/
        form_tentative_p(&vertices[level], &tentp[level], mgl[0].near_kernel_basis,
                         level+1, num_aggregations[level]);
        
        /* -- Smoothing -- */
        smooth_agg(&mgl[level].A, &tentp[level], &mgl[level].P, param,
                   level+1, &Neighbor[level]);
        
        /*-- Perform coarsening only up to the specified level --*/
        if ( mgl[level].P.col < MIN_CDOF ) break;
        
        /*-- Form restriction --*/
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        
        /*-- Form coarse level stiffness matrix --*/
        fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        
        fasp_dcsr_free(&Neighbor[level]);
        fasp_dcsr_free(&tentp[level]);
        fasp_ivec_free(&vertices[level]);
        
        ++level;
        
#if DIAGONAL_PREF
        // reorder each row to make diagonal appear first
        fasp_dcsr_diagpref(&mgl[level].A);
#endif
        
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK format
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[level].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[level].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
#if WITH_MUMPS
    /* Setup MUMPS direct solver on the coarsest level */
    fasp_solver_mumps_steps(&mgl[level].A, &mgl[level].b, &mgl[level].x, 1);
#endif
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w          = fasp_dvec_create(m);
    
    for ( level = 1; level < max_levels; ++level) {
        INT mm = mgl[level].A.row;
        mgl[level].num_levels = max_levels;
        mgl[level].b          = fasp_dvec_create(mm);
        mgl[level].x          = fasp_dvec_create(mm);
        
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[level].w = fasp_dvec_create(3*mm);
        else
            mgl[level].w = fasp_dvec_create(2*mm);
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
 * \brief Set up phase of plain aggregation AMG, using smoothed P and unsmoothed R
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
static SHORT amg_setup_smoothP_unsmoothR (AMG_data *mgl,
                                          AMG_param *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;
    
    // local variables
    SHORT       max_levels = param->max_levels, level=0, status=FASP_SUCCESS;
    INT         i, j;
    REAL        setup_start, setup_end;
    ILU_param   iluparam;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    fasp_gettime(&setup_start);
    
    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    // each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));
    
    // each level stores the information of the strongly coupled neighborhoods
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
    
    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)fasp_mem_calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
    }
    
    // Main AMG setup loop
    while ( (mgl[level].A.row > min_cdof) && (level < max_levels-1) ) {
        
        /*-- setup ILU decomposition if necessary */
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A, &mgl[level].LU, &iluparam);
            if ( status < 0 ) {
                param->ILU_levels = level;
                printf("### WARNING: ILU setup on level%d failed!\n", level);
            }
        }
        
        /* -- setup Schwarz smoother if necessary */
		if ( level < param->schwarz_levels ) {
            const INT smmsize = param->schwarz_mmsize;
            const INT smaxlvl = param->schwarz_maxlvl;
            const INT schtype = param->schwarz_type;
            
            mgl->schwarz_levels  = param->schwarz_levels;
            mgl[level].schwarz.A = fasp_dcsr_sympat(&mgl[level].A);
            fasp_dcsr_shift(&(mgl[level].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[level].schwarz, smmsize, smaxlvl, schtype);
		}
        
        /*-- Aggregation --*/
        status = aggregation_vmb(&mgl[level].A, &vertices[level], param, level+1,
                                 &Neighbor[level], &num_aggregations[level]);
        
        if ( status < 0 ) {
            // When error happens, force solver to use the current multigrid levels!
            printf("### WARNING: Coarsening on level%d is not successful!\n", level);
            status = FASP_SUCCESS; break;
        }

        /* -- Form Tentative prolongation --*/
        form_tentative_p(&vertices[level], &tentp[level], mgl[0].near_kernel_basis,
                         level+1, num_aggregations[level]);
        
        /* -- Smoothing -- */
        smooth_agg(&mgl[level].A, &tentp[level], &mgl[level].P, param,
                   level+1, &Neighbor[level]);
        
        /*-- Perform coarsening only up to the specified level --*/
        if ( mgl[level].P.col < MIN_CDOF ) break;
        
        /*-- Form resitriction --*/
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        fasp_dcsr_trans(&tentp[level], &tentr[level]);
        
        /*-- Form coarse level stiffness matrix --*/
        fasp_blas_dcsr_rap_agg(&tentr[level], &mgl[level].A, &tentp[level], &mgl[level+1].A);
        
        fasp_dcsr_free(&Neighbor[level]);
        fasp_dcsr_free(&tentp[level]);
        fasp_ivec_free(&vertices[level]);
        
        ++level;
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK format
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[level].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[level].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
#if WITH_MUMPS
    /* Setup MUMPS direct solver on the coarsest level */
    fasp_solver_mumps_steps(&mgl[level].A, &mgl[level].b, &mgl[level].x, 1);
#endif
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w          = fasp_dvec_create(m);
    
    for ( level = 1; level < max_levels; ++level) {
        INT mm = mgl[level].A.row;
        mgl[level].num_levels = max_levels;
        mgl[level].b          = fasp_dvec_create(mm);
        mgl[level].x          = fasp_dvec_create(mm);
        
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[level].w = fasp_dvec_create(3*mm);
        else
            mgl[level].w = fasp_dvec_create(2*mm);
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
    const SHORT prtlvl   = param->print_level;
    const SHORT min_cdof = MAX(param->coarse_dof,50);
    const INT   m        = mgl[0].A.ROW;
    const INT   nb       = mgl[0].A.nb;
    
    ILU_param iluparam;
    SHORT     max_levels=param->max_levels;
    SHORT     i, level=0, status=FASP_SUCCESS;
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
    
    //each elvel stores the information of the number of aggregations
    INT *num_aggs = (INT *)fasp_mem_calloc(max_levels, sizeof(INT));
    
    // each level stores the information of the strongly coupled neighborhoods
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

    if (param->aggregation_type == PAIRWISE) param->pair_number = MIN(param->pair_number, max_levels);
    
    // Main AMG setup loop
    while ( (mgl[level].A.ROW > min_cdof) && (level < max_levels-1) ) {
        
        printf("level = %d\n", level);
        
        /*-- setup ILU decomposition if necessary */
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dbsr_setup(&mgl[level].A, &mgl[level].LU, &iluparam);
            if ( status < 0 ) {
                param->ILU_levels = level;
                printf("### WARNING: ILU setup on level%d failed!\n", level);
            }
        }
        
        /*-- get the diagonal inverse --*/
        mgl[level].diaginv = fasp_dbsr_getdiaginv(&mgl[level].A);
        
        /*-- Aggregation --*/
        //mgl[level].PP =  fasp_dbsr_getblk_dcsr(&mgl[level].A);
        mgl[level].PP = fasp_dbsr_Linfinity_dcsr(&mgl[level].A);
        
        switch ( param->aggregation_type ) {
                
            case VMB: // VMB aggregation
                
                status = aggregation_vmb(&mgl[level].PP, &vertices[level], param, level+1,
                            &Neighbor[level], &num_aggs[level]);
                
                /*-- Choose strenth threshold adaptively --*/
                if ( num_aggs[level]*4 > mgl[level].PP.row )
                    param->strong_coupled /= 4;
                else if ( num_aggs[level]*1.25 < mgl[level].PP.row )
                    param->strong_coupled *= 1.5;
                
                break;
                
            default: // pairwise matching aggregation
                
                //aggregation_coarsening(&mgl[level].PP, param, level, vertices, &num_aggs[level]);
                status = aggregation_vmb(&mgl[level].PP, &vertices[level], param, level+1,
                            &Neighbor[level], &num_aggs[level]);
                
                break;
        }

        if ( status < 0 ) {
            // When error happens, force solver to use the current multigrid levels!
            printf("### WARNING: Coarsening on level%d is not successful!\n", level);
            status = FASP_SUCCESS; break;
        }
        
        /* -- Form Tentative prolongation --*/
        printf("before form tentative P\n");
        if (level == 0 && mgl[0].near_kernel_dim >0 ){
            form_tentative_p_bsr1(&vertices[level], &tentp[level], &mgl[0], level+1, num_aggs[level], mgl[0].near_kernel_dim, mgl[0].near_kernel_basis);
        }
        else{
            form_boolean_p_bsr(&vertices[level], &tentp[level], &mgl[0], level+1, num_aggs[level]);
        }
        
        /* -- Smoothing -- */
        printf("Smoothing P\n");
        smooth_agg_bsr(&mgl[level].A, &tentp[level], &mgl[level].P, param,
                   level+1, &Neighbor[level]);
        
        /*-- Form resitriction --*/
        printf("Form Restriction\n");
        fasp_dbsr_trans(&mgl[level].P, &mgl[level].R);
        
        /*-- Form coarse level stiffness matrix --*/
        printf("RAP\n");
        fasp_blas_dbsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        
        /* -- Form extra near kernal space if needed --*/
        printf("Form extra near kernel space\n");
        if (mgl[level].A_nk != NULL){
            
            mgl[level+1].A_nk = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            mgl[level+1].P_nk = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            mgl[level+1].R_nk = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
            
            temp1 = fasp_format_dbsr_dcsr(&mgl[level].R);
            fasp_blas_dcsr_mxm(&temp1, mgl[level].P_nk, mgl[level+1].P_nk);
            fasp_dcsr_trans(mgl[level+1].P_nk, mgl[level+1].R_nk);
            temp2 = fasp_format_dbsr_dcsr(&mgl[level+1].A);
            fasp_blas_dcsr_rap(mgl[level+1].R_nk, &temp2, mgl[level+1].P_nk, mgl[level+1].A_nk);
            fasp_dcsr_free(&temp1);
            fasp_dcsr_free(&temp2);
            
        }
        
        printf("clean\n");
        fasp_dcsr_free(&Neighbor[level]);
        fasp_ivec_free(&vertices[level]);
        fasp_dbsr_free(&tentp[level]);
        
        ++level;
    }
    
#if WITH_SuperLU
    /* Setup SuperLU direct solver on the coarsest level */
    mgl[level].Ac = fasp_format_dbsr_dcsr(&mgl[level].A);
#endif
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK format
    mgl[level].Ac = fasp_format_dbsr_dcsr(&mgl[level].A);
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[level].Ac, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran, &mgl[level].Ac);
    fasp_dcsr_free(&Ac_tran);
    
    mgl[level].Numeric = fasp_umfpack_factorize(&mgl[level].Ac, 5);
#endif
    
#if WITH_MUMPS
    /* Setup MUMPS direct solver on the coarsest level */
    mgl[level].Ac = fasp_format_dbsr_dcsr(&mgl[level].A);
    fasp_solver_mumps_steps(&mgl[level].Ac, &mgl[level].b, &mgl[level].x, 1);
#endif
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(3*m*nb);
    
    if (mgl[0].A_nk != NULL){
        
#if WITH_UMFPACK
        // Need to sort the matrix A_nk for UMFPACK
        fasp_dcsr_trans(mgl[0].A_nk, &temp1);
        fasp_dcsr_sort(&temp1);
        fasp_dcsr_cp(&temp1, mgl[0].A_nk);
        fasp_dcsr_free(&temp1);
#endif
        
    }

    
    for ( level = 1; level < max_levels; level++ ) {
        INT mm = mgl[level].A.ROW*nb;
        mgl[level].num_levels = max_levels;
        mgl[level].b = fasp_dvec_create(mm);
        mgl[level].x = fasp_dvec_create(mm);
        mgl[level].w = fasp_dvec_create(3*mm);
        
        if (mgl[level].A_nk != NULL){
         
#if WITH_UMFPACK
            // Need to sort the matrix A_nk for UMFPACK
            fasp_dcsr_trans(mgl[level].A_nk, &temp1);
            fasp_dcsr_sort(&temp1);
            fasp_dcsr_cp(&temp1, mgl[level].A_nk);
            fasp_dcsr_free(&temp1);
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
