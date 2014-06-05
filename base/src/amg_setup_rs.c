/*! \file amg_setup_rs.c
 *
 *  \brief Ruge-Stuben AMG: SETUP phase
 *
 *  \note Setup A, P, R, levels using classic AMG!
 *        Refter to "Multigrid" by Stuben
 *        in U. Trottenberg, C. W. Oosterlee and A. Schuller.
 *        Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *        Academic Press Inc., San Diego, CA, 2001.
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
 * \fn INT fasp_amg_setup_rs (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of Ruge and Stuben's classic AMG
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \return       FASP_SUCCESS if successed, otherwise, error information.
 *
 * \author Chensong Zhang
 * \date   05/09/2010
 *
 * Modified by Chensong Zhang on 04/04/2009.
 * Modified by Chensong Zhang on 05/09/2010.
 * Modified by Zhiyang Zhou on 11/17/2010.
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle.
 * Modified by Chensong zhang on 09/09/2011: add min dof.
 * Modified by Xiaozhe Hu on 04/24/2013: aggressive coarsening.
 * Modified by Chensong Zhang on 05/03/2013: add error handling in setup.
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
INT fasp_amg_setup_rs (AMG_data *mgl,
                       AMG_param *param)
{
    const SHORT  prtlvl     = param->print_level;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  min_cdof   = MAX(param->coarse_dof,MIN_CDOF);
    const INT    m          = mgl[0].A.row;
    
    // local variables
    INT        level = 0, status = FASP_SUCCESS;
    INT        max_levels = param->max_levels;
    REAL       setup_start, setup_end;
    ILU_param  iluparam;
    iCSRmat    S; // strong n-couplings
    
    // level info (fine: 0; coarse: 1)
    ivector    vertices = fasp_ivec_create(m);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_rs ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    fasp_gettime(&setup_start);
    
    // Make sure classical AMG will not call fasp_blas_dcsr_mxv_agg !!!
    param->tentative_smooth = 1.0;
    
    // If user want to use aggressive coarsening but did not specify number of
    // levels use aggressive coarsening, make sure apply aggresive coarsening
    // on the finest level only !!!
    if ( param->coarsening_type == COARSE_AC ) {
        param->aggressive_level = MAX(param->aggressive_level, 1);
    }
    
    // Initialize AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)fasp_mem_calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0, lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
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
    
#if DIAGONAL_PREF
    // Reorder each row to keep the diagonal entries appear first !!!
    fasp_dcsr_diagpref(&mgl[0].A);
#endif
    
    // Main AMG setup loop
    while ( (mgl[level].A.row > min_cdof) && (level < max_levels-1) ) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",
               level, mgl[level].A.row, mgl[level].A.nnz);
#endif
        
        /*-- Setup ILU decomposition if needed --*/
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A, &mgl[level].LU, &iluparam);
            if ( status < 0 ) {
                param->ILU_levels = level;
                printf("### WARNING: ILU setup on level %d failed!\n", level);
            }
        }
        
        /*-- Setup Schwarz smoother if needed --*/
        if ( level < param->schwarz_levels ) {
            const INT smmsize = param->schwarz_mmsize;
            const INT smaxlvl = param->schwarz_maxlvl;
            const INT schtype = param->schwarz_type;
            
            mgl->schwarz_levels  = param->schwarz_levels;
            mgl[level].schwarz.A = fasp_dcsr_sympat(&mgl[level].A);
            fasp_dcsr_shift(&(mgl[level].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[level].schwarz, smmsize, smaxlvl, schtype);
        }
        
        /*-- Coarseing and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[level].A, &vertices, &mgl[level].P, &S, param);
        if ( status < 0 ) {
            if ( prtlvl > PRINT_NONE ) {
                printf("### ERROR: Coarsening on level %d is not successful!\n", level);
            }
            fasp_mem_free(S.IA);
            fasp_mem_free(S.JA);
            return status;
        }
        
        /*-- Perform aggressive coarsening only up to the specified level --*/
        if ( mgl[level].P.col < MIN_CDOF ) break;
        if ( mgl[level].P.col*1.5 > mgl[level].A.row ) param->coarsening_type = COARSE_RS;
        if ( level == param->aggressive_level ) param->coarsening_type = COARSE_RS;
        
        /*-- Store the C/F marker --*/
        {
            INT size = mgl[level].A.row;
            mgl[level].cfmark = fasp_ivec_create(size);
            memcpy(mgl[level].cfmark.val, vertices.val, size*sizeof(INT));
        }
        
        /*-- Form interpolation --*/
        fasp_amg_interp(&mgl[level].A, &vertices, &mgl[level].P, &S, param);

        /*-- Form coarse level matrix: two RAP routines available! --*/
#if TRUE
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        fasp_blas_dcsr_rap (&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
#else
        fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        // TODO: Make a new ptap using (A,P) only. R is not needed as an input! --Chensong
#endif
        
        /*-- Clean up S! S is generated in coarsening --*/
        fasp_mem_free(S.IA);
        fasp_mem_free(S.JA);
        
        ++level;
        
#if DIAGONAL_PREF
        // reorder each row to make diagonal appear first
        fasp_dcsr_diagpref(&mgl[level].A);
#endif
        
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK to work
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[level].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran, &mgl[level].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
#if WITH_MUMPS
    // Setup MUMPS direct solver on the coarsest level
    fasp_solver_mumps_steps(&mgl[level].A, &mgl[level].b, &mgl[level].x, 1);
#endif
	
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w          = fasp_dvec_create(m);
    
#if FALSE
    INT groups;
    INT *results;
    for (level=0; level<max_levels && mgl[level].A.row>2000; level++) {
        results=(INT *)malloc(sizeof(INT)*mgl[level].A.row);
        dCSRmat_Division_Groups(mgl[level].A, results, &groups);
        printf("row = %d groups = %d \n", mgl[level].A.row, groups);
        free(results);
    }
#endif
    
    for ( level = 1; level < max_levels; ++level ) {
        INT mm = mgl[level].A.row;
        mgl[level].num_levels = max_levels;
        mgl[level].b          = fasp_dvec_create(mm);
        mgl[level].x          = fasp_dvec_create(mm);
        
        // allocate work arraies for the solve phase
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[level].w = fasp_dvec_create(3*mm);
        else
            mgl[level].w = fasp_dvec_create(2*mm);
    }
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        print_amgcomplexity(mgl, prtlvl);
        print_cputime("Classical AMG setup", setup_end - setup_start);
    }
    
    fasp_ivec_free(&vertices);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_rs ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn INT fasp_amg_setup_rs_omp (AMG_data *mgl, AMG_param *param)
 *
 * \brief  Setup of AMG based on R-S coarsening
 *
 * \param mgl    Pointer to AMG_data data
 * \param param  Pointer to AMG parameters
 *
 * \return       FASP_SUCCESS if successed, otherwise, error information.
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   03/11/2011
 */
INT fasp_amg_setup_rs_omp (AMG_data *mgl,
                           AMG_param *param)
{
    INT status = FASP_SUCCESS;
    
#ifdef _OPENMP
    const INT prtlvl      = param->print_level;
    const INT m           = mgl[0].A.row;
    const INT cycle_type  = param->cycle_type;
    const INT interp_type = param->interpolation_type;
    
	//local variables
    INT mm, size, nthreads;
    INT level = 0;
    INT max_levels = param->max_levels;
    REAL setup_start, setup_end, setup_duration;
    ILU_param iluparam;
    
    // set thread number
    nthreads = FASP_GET_NUM_THREADS();
    
    // stores level info (fine: 0; coarse: 1)
    ivector vertices = fasp_ivec_create(m);
    INT *icor_ysk    = (INT *)fasp_mem_calloc(5*nthreads+2, sizeof(INT));
    memset(icor_ysk, 0x0, sizeof(INT)*(5*nthreads+2));
    
    // strong n-couplings
    iCSRmat S;
    
#if DEBUG_MODE
    printf("fasp_amg_setup_rs ...... [Start]\n");
    printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    param->tentative_smooth = 1.0;
    
    fasp_gettime(&setup_start);
    
    //setup AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        const INT amlideg = param->amli_degree;
        param->amli_coef = (REAL *)fasp_mem_calloc(amlideg+1,sizeof(REAL));
        REAL lambda_max = 2.0;
        REAL lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, amlideg, param->amli_coef);
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
	mgl->schwarz_levels = param->schwarz_levels;
	INT schwarz_mmsize  = param->schwarz_mmsize;
    INT schwarz_maxlvl  = param->schwarz_maxlvl;
	INT schwarz_type    = param->schwarz_type;
    
#if DIAGONAL_PREF
	fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
	// main AMG setup loop
    while ( (mgl[level].A.row>MAX(param->coarse_dof,50)) && (level<max_levels-1) ) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
#endif
        
        /*-- setup ILU decomposition if necessary --*/
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
            if ( status < 0 ) goto FINISHED;
        }
        
        /*-- setup Schwarz smoother if necessary --*/
        if ( level < param->schwarz_levels ) {
            mgl[level].schwarz.A=fasp_dcsr_sympat(&mgl[level].A);
            fasp_dcsr_shift (&(mgl[level].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[level].schwarz, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
        }
        
        /*-- Coarseing and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[level].A, &vertices, &mgl[level].P, &S, param);
        if ( status < 0 ) goto FINISHED;
        
        size = mgl[level].A.row;
        mgl[level].cfmark = fasp_ivec_create(size);
        fasp_iarray_cp(size, vertices.val, mgl[level].cfmark.val);
        
        if (mgl[level].P.col == 0) break;        
        
        //status = fasp_amg_interp1(&mgl[level].A, &vertices, &mgl[level].P, param, &S, icor_ysk);
        fasp_amg_interp1(&mgl[level].A, &vertices, &mgl[level].P, param, &S, icor_ysk);

        //if ( status < 0 ) goto FINISHED;
        
        /*-- Form coarse level stiffness matrix --*/
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        
        if ( interp_type == INTERP_DIR ) {
            fasp_blas_dcsr_rap4(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, icor_ysk);
        }
        else {
            fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        }
        
        /*-- clean up! --*/
        fasp_mem_free(S.IA);
        fasp_mem_free(S.JA);
        
        ++level;
        
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif
        
    }
    
    fasp_mem_free(icor_ysk);
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w          = fasp_dvec_create(m);
    
    for ( level = 1; level < max_levels; ++level ) {
        mm                    = mgl[level].A.row;
        mgl[level].num_levels = max_levels;
        mgl[level].b          = fasp_dvec_create(mm);
        mgl[level].x          = fasp_dvec_create(mm);
        
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[level].w = fasp_dvec_create(3*mm);
        else
            mgl[level].w = fasp_dvec_create(2*mm);
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK to work
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
#if WITH_MUMPS
    // Setup MUMPS direct solver on the coarsest level
    fasp_solver_mumps_steps(&mgl[max_levels-1].A, &mgl[max_levels-1].b, &mgl[max_levels-1].x, 1);
#endif
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_amgcomplexity(mgl,prtlvl);
        print_cputime("Classical AMG setup", setup_duration);
    }
    
    status = FASP_SUCCESS;
    
FINISHED:
    fasp_ivec_free(&vertices);
    
#if DEBUG_MODE
	printf("### DEBUG: fasp_amg_setup_rs ...... [Finish]\n");
#endif
    
#endif  // end of OMP
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
