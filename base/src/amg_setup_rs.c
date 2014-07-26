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
 * \fn SHORT fasp_amg_setup_rs (AMG_data *mgl, AMG_param *param)
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
 * Modified by Chensong Zhang on 07/26/2014: handle coarsening errors.
 */
SHORT fasp_amg_setup_rs (AMG_data *mgl,
                         AMG_param *param)
{
    const SHORT  prtlvl     = param->print_level;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  min_cdof   = MAX(param->coarse_dof,MIN_CDOF);
    const INT    m          = mgl[0].A.row;
    
    // local variables
    SHORT      status = FASP_SUCCESS;
    INT        lvl = 0, max_lvls = param->max_levels;
    REAL       setup_start, setup_end;
    ILU_param  iluparam;
    iCSRmat    Scouple; // strong n-couplings
    
    // level info (fine: 0; coarse: 1)
    ivector    vertices = fasp_ivec_create(m);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nrow = %d, ncol = %d, nnz = %d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    fasp_gettime(&setup_start);
    
    // Make sure classical AMG will not call fasp_blas_dcsr_mxv_agg!
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
    
#if DIAGONAL_PREF
    // Reorder each row to keep the diagonal entries appear first !!!
    fasp_dcsr_diagpref(&mgl[0].A);
#endif
    
    // Main AMG setup loop
    while ( (mgl[lvl].A.row > min_cdof) && (lvl < max_lvls-1) ) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %d, row = %d, nnz = %d\n",
               lvl, mgl[lvl].A.row, mgl[lvl].A.nnz);
#endif
        
        /*-- Setup ILU decomposition if needed --*/
        if ( lvl < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[lvl].A, &mgl[lvl].LU, &iluparam);
            if ( status < 0 ) {
                param->ILU_levels = lvl;
                printf("### WARNING: ILU setup on level=%d failed!\n", lvl);
            }
        }
        
        /*-- Setup Schwarz smoother if needed --*/
        if ( lvl < param->schwarz_levels ) {
            const INT smmsize = param->schwarz_mmsize;
            const INT smaxlvl = param->schwarz_maxlvl;
            const INT schtype = param->schwarz_type;
            
            mgl->schwarz_levels  = param->schwarz_levels;
            mgl[lvl].schwarz.A = fasp_dcsr_sympat(&mgl[lvl].A);
            fasp_dcsr_shift(&(mgl[lvl].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[lvl].schwarz, smmsize, smaxlvl, schtype);
        }
        
        /*-- Coarseing and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[lvl].A, &vertices, &mgl[lvl].P, &Scouple, param);
        if ( status < 0 ) {
            if ( lvl > 0 ) { // Cannot continue coarsening, use the current levels!
                status = FASP_SUCCESS; break;
            }
            else {
                printf("### WARNING: Coarsening on level=%d failed!\n", lvl);
                fasp_mem_free(Scouple.IA);
                fasp_mem_free(Scouple.JA);
                return status;
            }
        }
        
        /*-- Perform aggressive coarsening only up to the specified level --*/
        if ( mgl[lvl].P.col < MIN_CDOF ) break;
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
        fasp_blas_dcsr_rap (&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);
        
        // TODO: Make a new ptap using (A,P) only. R is not needed as an input! --Chensong
        // fasp_blas_dcsr_ptap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);

        /*-- Clean up Scouple generated in coarsening --*/
        fasp_mem_free(Scouple.IA);
        fasp_mem_free(Scouple.JA);
        
        ++lvl;
        
#if DIAGONAL_PREF
        // reorder each row to make diagonal appear first
        fasp_dcsr_diagpref(&mgl[lvl].A);
#endif
        
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK to work
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[lvl].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran, &mgl[lvl].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
#if WITH_MUMPS
    // Setup MUMPS direct solver on the coarsest level
    fasp_solver_mumps_steps(&mgl[lvl].A, &mgl[lvl].b, &mgl[lvl].x, 1);
#endif
	
    // setup total level number and current level
    mgl[0].num_levels = max_lvls = lvl+1;
    mgl[0].w          = fasp_dvec_create(m);
    
/*    
    INT groups;
    INT *results;
    for (lvl=0; lvl<max_lvls && mgl[lvl].A.row>2000; lvl++) {
        results=(INT *)malloc(sizeof(INT)*mgl[lvl].A.row);
        dCSRmat_Division_Groups(mgl[lvl].A, results, &groups);
        printf("row = %d groups = %d \n", mgl[lvl].A.row, groups);
        free(results);
    }
*/
    
    for ( lvl = 1; lvl < max_lvls; ++lvl ) {
        INT mm = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_lvls;
        mgl[lvl].b          = fasp_dvec_create(mm);
        mgl[lvl].x          = fasp_dvec_create(mm);
        
        // allocate work arraies for the solve phase
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
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/**
 * \fn SHORT fasp_amg_setup_rs_omp (AMG_data *mgl, AMG_param *param)
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
 *
 * Modified by Chensong Zhang on 07/26/2014: handle coarsening errors.
 */
SHORT fasp_amg_setup_rs_omp (AMG_data *mgl,
                             AMG_param *param)
{
    SHORT status = FASP_SUCCESS;
    
#ifdef _OPENMP
    const INT prtlvl      = param->print_level;
    const INT m           = mgl[0].A.row;
    const INT cycle_type  = param->cycle_type;
    const INT interp_type = param->interpolation_type;
    
	//local variables
    INT mm, size, nthreads;
    INT lvl = 0;
    INT max_lvls = param->max_levels;
    REAL setup_start, setup_end, setup_duration;
    ILU_param iluparam;
    
    // set thread number
    nthreads = FASP_GET_NUM_THREADS();
    
    // strong n-couplings
    iCSRmat Scouple;
    
    // stores level info (fine: 0; coarse: 1)
    ivector vertices = fasp_ivec_create(m);
    INT *icor_ysk    = (INT *)fasp_mem_calloc(5*nthreads+2, sizeof(INT));
    memset(icor_ysk, 0x0, sizeof(INT)*(5*nthreads+2));
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nrow = %d, ncol = %d, nnz = %d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    param->tentative_smooth = 1.0;
    
    fasp_gettime(&setup_start);
    
    //setup AMLI coefficients
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
	mgl->schwarz_levels = param->schwarz_levels;
	INT schwarz_mmsize  = param->schwarz_mmsize;
    INT schwarz_maxlvl  = param->schwarz_maxlvl;
	INT schwarz_type    = param->schwarz_type;
    
#if DIAGONAL_PREF
	fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
	// main AMG setup loop
    while ( (mgl[lvl].A.row>MAX(param->coarse_dof,50)) && (lvl<max_lvls-1) ) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %d, row = %d, nnz = %d\n",
               lvl, mgl[lvl].A.row, mgl[lvl].A.nnz);
#endif
        
        /*-- setup ILU decomposition if necessary --*/
        if ( lvl < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[lvl].A,&mgl[lvl].LU,&iluparam);
            if ( status < 0 ) {
                param->ILU_levels = lvl;
                printf("### WARNING: ILU setup on level %d failed!\n", lvl);
            }
        }
        
        /*-- setup Schwarz smoother if necessary --*/
        if ( lvl < param->schwarz_levels ) {
            mgl[lvl].schwarz.A=fasp_dcsr_sympat(&mgl[lvl].A);
            fasp_dcsr_shift (&(mgl[lvl].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[lvl].schwarz, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
        }
        
        /*-- Coarseing and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[lvl].A, &vertices, &mgl[lvl].P, &Scouple, param);
        if ( status < 0 ) {
            if ( lvl > 0 ) { // Cannot continue coarsening. Stop right here!
                status = FASP_SUCCESS; break;
            }
            else {
                printf("### WARNING: Coarsening on level=%d failed!\n", lvl);
                fasp_mem_free(Scouple.IA);
                fasp_mem_free(Scouple.JA);
                return status;
            }
        }
        
        size = mgl[lvl].A.row;
        mgl[lvl].cfmark = fasp_ivec_create(size);
        fasp_iarray_cp(size, vertices.val, mgl[lvl].cfmark.val);
        
        if (mgl[lvl].P.col == 0) break;        
        
        fasp_amg_interp1(&mgl[lvl].A, &vertices, &mgl[lvl].P, param, &Scouple, icor_ysk);
        
        /*-- Form coarse level stiffness matrix --*/
        fasp_dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);
        
        if ( interp_type == INTERP_DIR ) {
            fasp_blas_dcsr_rap4(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A, icor_ysk);
        }
        else {
            fasp_blas_dcsr_rap(&mgl[lvl].R, &mgl[lvl].A, &mgl[lvl].P, &mgl[lvl+1].A);
        }
        
        /*-- clean up! --*/
        fasp_mem_free(Scouple.IA);
        fasp_mem_free(Scouple.JA);
        
        ++lvl;
        
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[lvl].A); // reorder each row to make diagonal appear first
#endif
        
    }
    
    fasp_mem_free(icor_ysk);
    
    // setup total level number and current level
    mgl[0].num_levels = max_lvls = lvl+1;
    mgl[0].w          = fasp_dvec_create(m);
    
    for ( lvl = 1; lvl < max_lvls; ++lvl ) {
        mm                  = mgl[lvl].A.row;
        mgl[lvl].num_levels = max_lvls;
        mgl[lvl].b          = fasp_dvec_create(mm);
        mgl[lvl].x          = fasp_dvec_create(mm);
        
        if ( cycle_type == NL_AMLI_CYCLE )
            mgl[lvl].w = fasp_dvec_create(3*mm);
        else
            mgl[lvl].w = fasp_dvec_create(2*mm);
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK to work
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_lvls-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[max_lvls-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
#if WITH_MUMPS
    // Setup MUMPS direct solver on the coarsest level
    fasp_solver_mumps_steps(&mgl[max_lvls-1].A, &mgl[max_lvls-1].b, &mgl[max_lvls-1].x, 1);
#endif
    
    if ( prtlvl > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_amgcomplexity(mgl,prtlvl);
        print_cputime("Classical AMG setup", setup_duration);
    }
    
    fasp_ivec_free(&vertices);
    
#if DEBUG_MODE
	printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
#endif  // end of OMP
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
