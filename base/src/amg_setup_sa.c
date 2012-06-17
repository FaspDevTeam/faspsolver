/*! \file amg_setup_sa.c
 *  \brief Smoothed Aggregation AMG: SETUP phase
 *
 */

#include <math.h>
#include <time.h>
#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "amg_setup_aggregation.inl"

static SHORT amg_setup_smoothP_smoothA(AMG_data *mgl, AMG_param *param);
static SHORT amg_setup_smoothP_unsmoothA(AMG_data *mgl, AMG_param *param);
static void smooth_agg(dCSRmat *A, dCSRmat *tentp, dCSRmat *P, AMG_param *param, INT levelNum, dCSRmat *N);
static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl,AMG_param *param);
/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_amg_setup_sa (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of smoothed aggregation AMG
 * 
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \return        SUCCESS if succeed, error otherwise
 * 
 * \author Xiaozhe Hu
 * \date   09/29/2009 
 *
 * \note Setup A, P, PT, levels using smoothed aggregation concrete algorithm;
 *       Refer to Peter Vanek, Jan Madel and Marin Brezina, 
 *       "Algebraic Multigrid on Unstructured Meshes", 1994
 *
 * Modified by Chensong Zhang on 04/06/2010.
 * Modified by Chensong Zhang on 05/09/2010.
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle
 */
SHORT fasp_amg_setup_sa (AMG_data *mgl, 
                         AMG_param *param)
{    
	// only for test smoothed P and unsmoothed A, not used in general. 
    SHORT type = 0;  
    SHORT status=SUCCESS;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_sa ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    if (param->tentative_smooth > SMALLREAL) type = 2;
    switch (type) {
    case 1:
        status = amg_setup_smoothP_unsmoothA(mgl, param);
        break;
    case 2:
        status = amg_setup_smoothP_smoothA(mgl, param);
        break;    
	default:
		status = amg_setup_unsmoothP_unsmoothA(mgl, param);
		break;
    }
    
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
 * \brief Set up phase of smoothed aggregation AMG, using smoothed P and smoothed A
 * 
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011 
 */
static SHORT amg_setup_smoothP_smoothA (AMG_data *mgl, 
                                        AMG_param *param)
{
    const SHORT print_level=param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const INT   m=mgl[0].A.row;
    
    SHORT       max_levels=param->max_levels, level=0, status=SUCCESS;
    INT         i, j;
    clock_t     setup_start, setup_end;
    REAL        setupduration;
    
    setup_start=clock();
    
    if (cycle_type == AMLI_CYCLE) {
        param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
        REAL lambda_max = 2.0;
        REAL lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
    }
    
    //each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));
    
    for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
    
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    // each level stores the information of the tentative prolongations
    dCSRmat *tentp    = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(REAL*));
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        mgl[0].near_kernel_basis[i] = (REAL *)fasp_mem_calloc(m,sizeof(REAL));
        for (j=0;j<m;++j) mgl[0].near_kernel_basis[i][j] = 1.0;
    }
    
    // initialize ILU parameters
    ILU_param iluparam;
    if (param->ILU_levels>0) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }
    
    // initialize Schwarz parameters
	mgl->schwarz_levels = param->schwarz_levels;
	INT schwarz_mmsize = param->schwarz_mmsize;
	INT schwarz_maxlvl = param->schwarz_maxlvl;
	INT schwarz_type = param->schwarz_type;

    
#if DIAGONAL_PREF
    fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
    while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1)) {

#if DEBUG_MODE
        printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",
               level,mgl[level].A.row,mgl[level].A.nnz);
#endif
    
        /*-- setup ILU decomposition if necessary */
        if (level<param->ILU_levels) fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
        
        /* -- setup Schwarz smoother if necessary */
		if (level<param->schwarz_levels){
			mgl[level].schwarz.A=fasp_dcsr_sympat(&mgl[level].A);
			fasp_dcsr_shift (&(mgl[level].schwarz.A), 1);
			fasp_schwarz_setup(&mgl[level].schwarz, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
		}

    
        /*-- Aggregation --*/
        aggregation(&mgl[level].A, &vertices[level], param, level+1, &Neighbor[level], &num_aggregations[level]);
    
        /* -- Form Tentative prolongation --*/      
        form_tentative_p(&vertices[level], &tentp[level], &mgl[0], level+1, num_aggregations[level]);
    
        /* -- Smoothing -- */
        smooth_agg(&mgl[level].A, &tentp[level], &mgl[level].P, param, level+1, &Neighbor[level]);
    
        /*-- Form restriction --*/    
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
    
        /*-- Form coarse level stiffness matrix --*/    
        fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
    
        fasp_dcsr_free(&Neighbor[level]);
        fasp_dcsr_free(&tentp[level]);
        fasp_ivec_free(&vertices[level]);
    
        ++level;
        
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif      
    }
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(m);    
    
    for (level=1; level<max_levels; ++level) {
        INT    m = mgl[level].A.row;
        mgl[level].num_levels = max_levels;     
        mgl[level].b = fasp_dvec_create(m);
        mgl[level].x = fasp_dvec_create(m);
    
        if (cycle_type == NL_AMLI_CYCLE)  mgl[level].w = fasp_dvec_create(3*m);    
        else mgl[level].w = fasp_dvec_create(2*m);
    }
    
#if With_UMFPACK    
    // Need to sort the matrix A for UMFPACK format
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
    if (print_level>PRINT_NONE) {
        setup_end=clock();
        setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
        print_amgcomplexity(mgl,print_level);
        print_cputime("Smoothed Aggregation AMG setup",setupduration);
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
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011 
 */
static SHORT amg_setup_smoothP_unsmoothA (AMG_data *mgl,
                                          AMG_param *param)
{
    const SHORT print_level=param->print_level;
    const SHORT cycle_type = param->cycle_type;
    const INT   m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
    
    // local variables
    SHORT       max_levels=param->max_levels, level=0, status=SUCCESS;
    INT         i, j;
    clock_t     setup_start, setup_end;
    REAL        setupduration;
    
#if DEBUG_MODE
    printf("### DEBUG: amg_setup_sa ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if (print_level>8)    printf("amg_setup: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
    
    setup_start=clock();
    
    if (cycle_type == AMLI_CYCLE) 
        {
            param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
            REAL lambda_max = 2.0;
            REAL lambda_min = lambda_max/4;
            fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
        }
    
    //each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));
    
    for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
    
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    // each level stores the information of the tentative prolongations
    dCSRmat *tentp    = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    dCSRmat *tentpt    = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat));
    
    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(REAL*));
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        mgl[0].near_kernel_basis[i] = (REAL *)fasp_mem_calloc(m,sizeof(REAL));
        for (j=0;j<m;++j) mgl[0].near_kernel_basis[i][j] = 1.0;
    }
    
    // initialize ILU parameters
    ILU_param iluparam;
    if (param->ILU_levels>0) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }
    
    // initialize Schwarz parameters
	mgl->schwarz_levels = param->schwarz_levels;
	INT schwarz_mmsize = param->schwarz_mmsize;
	INT schwarz_maxlvl = param->schwarz_maxlvl;
	INT schwarz_type = param->schwarz_type;

    
    while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1)) {
        /*-- setup ILU decomposition if necessary */
        if (level<param->ILU_levels) fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
        
        /* -- setup Schwarz smoother if necessary */
		if (level<param->schwarz_levels){
			mgl[level].schwarz.A=fasp_dcsr_sympat(&mgl[level].A);
			fasp_dcsr_shift (&(mgl[level].schwarz.A), 1);
			fasp_schwarz_setup(&mgl[level].schwarz, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
		}

    
        /*-- Aggregation --*/
        aggregation(&mgl[level].A, &vertices[level], param, level+1, &Neighbor[level], &num_aggregations[level]);
    
        /* -- Form Tentative prolongation --*/      
        form_tentative_p(&vertices[level], &tentp[level], &mgl[0], level+1, num_aggregations[level]);
    
        /* -- Smoothing -- */
        smooth_agg(&mgl[level].A, &tentp[level], &mgl[level].P, param, level+1, &Neighbor[level]);
    
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
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(m);    
    
    for (level=1; level<max_levels; ++level) {
        INT    m = mgl[level].A.row;
        mgl[level].num_levels = max_levels;     
        mgl[level].b = fasp_dvec_create(m);
        mgl[level].x = fasp_dvec_create(m);
    
        if (cycle_type == NL_AMLI_CYCLE)  mgl[level].w = fasp_dvec_create(3*m);    
        else mgl[level].w = fasp_dvec_create(2*m);
    }
    
#if With_UMFPACK    
    // Need to sort the matrix A for UMFPACK format
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
    if (print_level>PRINT_MIN) {
        REAL gridcom=0.0, opcom=0.0;
    
        printf("-----------------------------------------------\n");
        printf("  Level     Num of rows      Num of nonzeros\n");
        printf("-----------------------------------------------\n");
        for (level=0;level<max_levels;++level) {
            printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
            gridcom += mgl[level].A.row;
            opcom += mgl[level].A.nnz;
        }
        printf("-----------------------------------------------\n");
    
        gridcom /= mgl[0].A.row;
        opcom /= mgl[0].A.nnz;
        printf("Half smoothed Aggregation AMG grid complexity = %f\n", gridcom);
        printf("Half smoothed Aggregation AMG operator complexity = %f\n", opcom);
    }
    
    if (print_level>PRINT_NONE) {
        setup_end=clock();
        setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
        print_cputime("Half Smoothed Aggregation AMG setup",setupduration);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    fasp_mem_free(tentp);
    
#if DEBUG_MODE
    printf("### DEBUG: amg_setup_sa ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn static void smooth_agg (dCSRmat *A, dCSRmat *tentp, dCSRmat *P, AMG_param *param, 
 *                             INT levelNum, dCSRmat *N)
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
    
    /* Step 1. Form smoother */
    
    /* Using A for damped Jacobian smoother */
    if (filter == 0){ 
        S = fasp_dcsr_create(row, col, A->IA[row]); // copy structure from A
        for (i=0; i<row+1; ++i) S.IA[i] = A->IA[i];
        for (i=0; i<S.IA[S.row]; ++i) S.JA[i] = A->JA[i];
    
        fasp_dcsr_getdiag(0, A, &diag);  // get the diaganol entries of A
    
        // check the diaganol entries. 
        // if it is too small, use Richardson smoother for the corresponding row 
        for (i=0;i<row;++i){
            if (ABS(diag.val[i]) < 1e-6){
                diag.val[i] = 1.0;
            }    
        }
    
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
        for (i=0; i<row; ++i){
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
    
        S = fasp_dcsr_create(row, col, N->IA[row]); // copy structure from N (filtered A)
        for (i=0; i<row+1; ++i) S.IA[i] = N->IA[i];
        for (i=0; i<S.IA[S.row]; ++i) S.JA[i] = N->JA[i];
    
        fasp_dcsr_getdiag(0, N, &diag);  // get the diaganol entries of N (filtered A)
    
        // check the diaganol entries. 
        // if it is too small, use Richardson smoother for the corresponding row 
        for (i=0;i<row;++i){
            if (ABS(diag.val[i]) < 1e-6){
                diag.val[i] = 1.0;
            }    
        }
    
        for (i=0;i<row;++i){
            for (j=S.IA[i]; j<S.IA[i+1]; ++j){
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


/**
 * \fn static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 *
 * \author        Chunsheng Feng, Xiaoqiang Yue
 * \date 03/15/2011
 *
 * \note Setup A, P, PT, levels using smoothed aggregation
 *       concrete algorithm see paper
 *       Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 */
static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl,
                                         AMG_param *param)
{
    INT status=SUCCESS;

    const INT print_level=param->print_level;
    const INT m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
    
    INT max_levels=param->max_levels;
    INT i, level=0;
    REAL setup_start, setup_end;
    REAL setupduration;
    
#if DEBUG_MODE
    printf("fasp_amg_setup_sa ...... [Start]\n");
    printf("fasp_amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if (print_level>8) printf("fasp_amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
   
#if FASP_USE_OPENMP
    setup_start=omp_get_wtime();
#endif
    
    if (param->cycle_type == AMLI_CYCLE) 
        {
            param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
            REAL lambda_max = 2.0;
            REAL lambda_min = lambda_max/4;
            fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
        }
   
	 //each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT)); 
    
    for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
    
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
	// each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); 
    
    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (REAL **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(REAL*));
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        mgl[0].near_kernel_basis[i] = (REAL *)fasp_mem_calloc(m,sizeof(REAL));
        fasp_array_set(m, mgl[0].near_kernel_basis[i], 1.0);
    }
    
    // initialize ILU parameters
    ILU_param iluparam;
    if (param->ILU_levels>0) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }
    
    while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
        {
            /*-- setup ILU decomposition if necessary */
            if (level<param->ILU_levels) fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
    
            /*-- Aggregation --*/
            aggregation(&mgl[level].A, &vertices[level], param, level+1, &Neighbor[level], &num_aggregations[level]);
            if (num_aggregations[level] * 4 > mgl[level].A.row) param->strong_coupled /=2.0; 
    
            /* -- Form Prolongation --*/      
            form_tentative_p(&vertices[level], &mgl[level].P, &mgl[0], level+1, num_aggregations[level]);
    
            /*-- Form resitriction --*/    
            fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
    
            /*-- Form coarse level stiffness matrix --*/    
            fasp_blas_dcsr_rap_agg(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);

            fasp_dcsr_free(&Neighbor[level]);
            fasp_ivec_free(&vertices[level]);
    
            ++level;
        }
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(m);    
    
    for (level=1; level<max_levels; ++level) {
        INT    m = mgl[level].A.row;
        mgl[level].num_levels = max_levels;     
        mgl[level].b = fasp_dvec_create(m);
        mgl[level].x = fasp_dvec_create(m);
        mgl[level].w = fasp_dvec_create(2*m);    
    }
    
#if With_UMFPACK
    // Need to sort the matrix A for UMFPACK format
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
    if (print_level>1) {
        REAL gridcom=0.0, opcom=0.0;
    
        printf("-----------------------------------------------\n");
        printf("  Level     Num of rows      Num of nonzeros\n");
        printf("-----------------------------------------------\n");
        for (level=0;level<max_levels;++level) {
            printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
            gridcom += mgl[level].A.row;
            opcom += mgl[level].A.nnz;
        }
        printf("-----------------------------------------------\n");
    
        gridcom /= mgl[0].A.row;
        opcom /= mgl[0].A.nnz;
        printf("Unsmoothed Aggregation AMG grid complexity = %f\n", gridcom);
        printf("Unsmoothed Aggregation AMG operator complexity = %f\n", opcom);
    }
    
    if (print_level>0) {
#if FASP_USE_OPENMP
        setup_end=omp_get_wtime();
#endif
        setupduration = setup_end - setup_start;
        print_cputime("Unsmoothed Aggregation AMG setup",setupduration);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    
#if DEBUG_MODE
    printf("amg_setup_sa ...... [Finish]\n");
#endif

    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
