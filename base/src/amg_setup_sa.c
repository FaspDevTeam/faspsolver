/*! \file amg_setup_sa.c
 *
 *  \brief Smoothed aggregation AMG: SETUP phase
 *
 *  \note Setup A, P, PT and levels using the unsmoothed aggregation algorithm;
 *        Refer to P. Vanek, J. Madel and M. Brezina
 *        "Algebraic Multigrid on Unstructured Meshes", 1994
 */

#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"
#include "aggregation_csr.inl"

static SHORT amg_setup_smoothP_smoothA(AMG_data *, AMG_param *);
static SHORT amg_setup_smoothP_unsmoothA(AMG_data *, AMG_param *);
static void smooth_agg(dCSRmat *, dCSRmat *, dCSRmat *, AMG_param *, INT, dCSRmat *);

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
 * \return       SUCCESS if succeed, error otherwise
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * Modified by Chensong Zhang on 04/06/2010.
 * Modified by Chensong Zhang on 05/09/2010.
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
SHORT fasp_amg_setup_sa (AMG_data *mgl,
                         AMG_param *param)
{
    SHORT status  = SUCCESS;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_sa ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n",
           mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
#if TRUE
    status = amg_setup_smoothP_smoothA(mgl, param);
#else // smoothed P, unsmoothed A
    status = amg_setup_smoothP_unsmoothA(mgl, param);
#endif
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_sa ...... [Finish]\n");
#endif
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static SHORT amg_setup_smoothP_smoothA (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of smoothed aggregation AMG, using smoothed P and smoothed A
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
static SHORT amg_setup_smoothP_smoothA (AMG_data *mgl,
                                        AMG_param *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;
    
    // local variables
    SHORT       max_levels = param->max_levels, level=0, status=SUCCESS;
    INT         i, j;
    REAL        setup_start, setup_end;
    ILU_param   iluparam;
    
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
    
    // Main AMG setup loop
    while ( (mgl[level].A.row > min_cdof) && (level < max_levels-1) ) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",
               level,mgl[level].A.row,mgl[level].A.nnz);
#endif
        
        /*-- setup ILU decomposition if necessary */
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A, &mgl[level].LU, &iluparam);
            if ( status < 0 ) {
                param->ILU_levels = level;
                printf("### WARNING: ILU setup on level %d failed!\n", level);
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
        aggregation(&mgl[level].A, &vertices[level], param, level+1,
                    &Neighbor[level], &num_aggregations[level]);
        
        /* -- Form Tentative prolongation --*/
        form_tentative_p(&vertices[level], &tentp[level], &mgl[0],
                         level+1, num_aggregations[level]);
        
        /* -- Smoothing -- */
        smooth_agg(&mgl[level].A, &tentp[level], &mgl[level].P, param,
                   level+1, &Neighbor[level]);
        
        /*-- Perform aggressive coarsening only up to the specified level --*/
        if ( mgl[level].P.col < 50 ) break; // If coarse < 50, stop!!!
        
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
    
    return status;
}

/**
 * \fn static SHORT amg_setup_smoothP_unsmoothA (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of plain aggregation AMG, using smoothed P and unsmoothed A
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 * Modified by Chensong Zhang on 05/10/2013: adjust the structure.
 */
static SHORT amg_setup_smoothP_unsmoothA (AMG_data *mgl,
                                          AMG_param *param)
{
    const SHORT prtlvl     = param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const SHORT min_cdof   = MAX(param->coarse_dof,50);
    const INT   m          = mgl[0].A.row;
    
    // local variables
    SHORT       max_levels = param->max_levels, level=0, status=SUCCESS;
    INT         i, j;
    REAL        setup_start, setup_end;
    ILU_param   iluparam;
    
    fasp_gettime(&setup_start);
    
    // level info (fine: 0; coarse: 1)
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    // each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));
    
    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    // each level stores the information of the tentative prolongations
    dCSRmat *tentp    = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    dCSRmat *tentpt   = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
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
                printf("### WARNING: ILU setup on level %d failed!\n", level);
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
        aggregation(&mgl[level].A, &vertices[level], param, level+1,
                    &Neighbor[level], &num_aggregations[level]);
        
        /* -- Form Tentative prolongation --*/
        form_tentative_p(&vertices[level], &tentp[level], &mgl[0],
                         level+1, num_aggregations[level]);
        
        /* -- Smoothing -- */
        smooth_agg(&mgl[level].A, &tentp[level], &mgl[level].P, param,
                   level+1, &Neighbor[level]);
        
        /*-- Perform aggressive coarsening only up to the specified level --*/
        if ( mgl[level].P.col < 50 ) break; // If coarse < 50, stop!!!
        
        /*-- Form resitriction --*/
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        fasp_dcsr_trans(&tentp[level], &tentpt[level]);
        
        /*-- Form coarse level stiffness matrix --*/
        fasp_blas_dcsr_rap_agg(&tentpt[level], &mgl[level].A, &tentp[level], &mgl[level+1].A);
        
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
    fasp_mem_free(tentpt);
    
    return status;
}

/**
 * \fn static void smooth_agg (dCSRmat *A, dCSRmat *tentp, dCSRmat *P,
 *                             AMG_param *param, INT levelNum, dCSRmat *N)
 *
 * \brief Smooth the tentative prolongation
 *
 * \param A         pointer to the coefficient matrices
 * \param tentp     pointer to the tentative prolongation operators
 * \param P         pointer to the prolongation operators
 * \param param     pointer to AMG parameters
 * \param levelNum  Current level number
 * \param N         pointer to strongly coupled neighborhoods
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * Modified by Chunsheng Feng, Zheng Li on 10/12/2012
 */
static void smooth_agg (dCSRmat *A,
                        dCSRmat *tentp,
                        dCSRmat *P,
                        AMG_param *param,
                        INT levelNum,
                        dCSRmat *N)
{
    const REAL smooth_factor = param->tentative_smooth;
    const SHORT filter = param->smooth_filter;
    
    INT row = A->row, col= A->col;
    INT i,j;
    
    dCSRmat S;
    dvector diag;  // diaganoal entries
    
    REAL row_sum_A, row_sum_N;
    
    // local variables
#ifdef _OPENMP
    INT myid, mybegin, myend;
    INT nthreads = FASP_GET_NUM_THREADS();
#endif
    
    /* Step 1. Form smoother */
    
    /* Using A for damped Jacobian smoother */
    if (filter == 0){
        
        S = fasp_dcsr_create(row, col, A->IA[row]); // copy structure from A
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
        for (i=0; i<row+1; ++i) S.IA[i] = A->IA[i];
        for (i=0; i<S.IA[S.row]; ++i) S.JA[i] = A->JA[i];
        
        fasp_dcsr_getdiag(0, A, &diag);  // get the diaganol entries of A
        
        // check the diaganol entries.
        // if it is too small, use Richardson smoother for the corresponding row
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
        for (i=0;i<row;++i){
            if (ABS(diag.val[i]) < 1e-6){
                diag.val[i] = 1.0;
            }
        }
        
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS) private(j)
#endif
        for (i=0;i<row;++i){
            for (j=S.IA[i]; j<S.IA[i+1]; ++j){
                if (S.JA[j] == i) {
                    S.val[j] = 1 -  smooth_factor * A->val[j] / diag.val[i];
                }
                else {
                    S.val[j] = - smooth_factor * A->val[j] / diag.val[i];
                }
            }
        }
    }
    
    /* Using filtered A for damped Jacobian smoother */
    else {
        /* Form filtered A and store in N */
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, j, row_sum_A, row_sum_N) if (row>OPENMP_HOLDS)
        for (myid=0; myid<nthreads; myid++) {
            FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) {
#else
            for (i=0; i<row; ++i) {
#endif
                row_sum_A = 0.0;
                row_sum_N = 0.0;
                    
                for (j=A->IA[i]; j<A->IA[i+1]; ++j){
                    if (A->JA[j] != i){
                        row_sum_A+=A->val[j];
                    }
                }
                for (j=N->IA[i]; j<N->IA[i+1]; ++j){
                    if (N->JA[j] != i){
                        row_sum_N+=N->val[j];
                    }
                }
                    
                for (j=N->IA[i]; j<N->IA[i+1]; ++j){
                    if (N->JA[j] == i){
                        N->val[j] = N->val[j] - row_sum_A + row_sum_N;
                    }
                }
            }
#ifdef _OPENMP
        }
#endif
        S = fasp_dcsr_create(row, col, N->IA[row]); // copy structure from N (filtered A)
            
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
        for (i=0; i<row+1; ++i) S.IA[i] = N->IA[i];
        for (i=0; i<S.IA[S.row]; ++i) S.JA[i] = N->JA[i];
            
        fasp_dcsr_getdiag(0, N, &diag);  // get the diaganol entries of N (filtered A)
            
        // check the diaganol entries.
        // if it is too small, use Richardson smoother for the corresponding row
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
        for (i=0;i<row;++i) {
            if (ABS(diag.val[i]) < 1e-6) diag.val[i] = 1.0;
        }
            
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS) private(i,j)
#endif
        for (i=0;i<row;++i) {
            for (j=S.IA[i]; j<S.IA[i+1]; ++j) {
                if (S.JA[j] == i) {
                    S.val[j] = 1 -  smooth_factor * N->val[j] / diag.val[i];
                }
                else {
                    S.val[j] = - smooth_factor * N->val[j] / diag.val[i];
                }
            }
        }
            
    }
        
    fasp_dvec_free(&diag);
        
    /* Step 2. Smooth the tentative prolongation */
    /* P = S*tenp */
    fasp_blas_dcsr_mxm(&S, tentp, P); // Note: think twice about this.
        
    P->nnz = P->IA[P->row];
        
    fasp_dcsr_free(&S);
}
    
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
