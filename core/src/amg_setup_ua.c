/*! \file amg_setup_ua.c
 *  \brief Unsmoothed Aggregation AMG: SETUP phase
 *
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "amg_setup_aggregation.inl"

static SHORT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param);
static SHORT amg_setup_unsmoothP_unsmoothA_bsr(AMG_data_bsr *mgl, AMG_param * param);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_amg_setup_ua(AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of unsmoothed aggregation AMG
 * 
 * \param mgl     Pointer to AMG_data data
 * \param param   Pointer to AMG parameters
 *
 * \return        SUCCESS if succeed, error otherwise
 *
 * \note Setup A, P, PT, levels using unsmoothed aggregation aggregation algorithm;
 *       Refer to Peter Vanek, Jan Madel and Marin Brezina
 *       "Algebraic Multigrid on Unstructured Meshes", 1994
 * 
 * \author Xiaozhe Hu
 * \date   12/28/2011 
 */
SHORT fasp_amg_setup_ua (AMG_data *mgl, 
                         AMG_param *param)
{
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_ua ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    SHORT status = amg_setup_unsmoothP_unsmoothA(mgl, param);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_ua ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn INT fasp_amg_setup_ua_bsr(AMG_data_bsr *mgl, AMG_param *param)
 * \brief Set up phase of unsmoothed aggregation AMG (BSR format) 
 * 
 * \param *mgl     pointer to AMG_data_bsr data
 * \param *param   pointer to AMG parameters
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 03/16/2012 
 *
 */
SHORT fasp_amg_setup_ua_bsr (AMG_data_bsr *mgl, 
                             AMG_param *param)
{
    SHORT status=SUCCESS;
    
    status = amg_setup_unsmoothP_unsmoothA_bsr(mgl, param);
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static SHORT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param mgl     Pointer to AMG_data data
 * \param param   Pointer to AMG parameters
 * 
 * \author Xiaozhe Hu
 * \date   02/21/2011 
 */
static SHORT amg_setup_unsmoothP_unsmoothA (AMG_data *mgl, 
                                            AMG_param *param)
{
    const SHORT   cycle_type=param->cycle_type;
    const SHORT   print_level=param->print_level;
    const INT     m=mgl[0].A.row;
    
    // local variables
    clock_t       setup_start, setup_end;
    REAL          setupduration;
    SHORT         level=0, status=SUCCESS, max_levels=param->max_levels;
    INT           i, j;
    
    setup_start = clock();
    
    if (cycle_type == AMLI_CYCLE) {
        param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
        REAL lambda_max = 2.0;
        REAL lambda_min = lambda_max/8;
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
        for (j=0;j<m;++j) mgl[0].near_kernel_basis[i][j] = 1.0;
    }
    
    // initialize ILU parameters
    mgl->ILU_levels = param->ILU_levels;
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
        print_cputime("Unsmoothed Aggregation AMG setup",setupduration);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    
    return status;
}

/**
 * \fn static SHORT amg_setup_unsmoothP_unsmoothA_bsr(AMG_data_bsr *mgl, AMG_param *param)
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A (BSR format)
 * 
 * \param *mgl     pointer to AMG_data_bsr data
 * \param *param   pointer to AMG parameters
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 03/16/2012 
 */
static SHORT amg_setup_unsmoothP_unsmoothA_bsr(AMG_data_bsr *mgl, AMG_param *param)
{
    const SHORT print_level=param->print_level;
    const INT m=mgl[0].A.ROW, n=mgl[0].A.COL, nnz=mgl[0].A.NNZ, nb = mgl[0].A.nb;
    
    SHORT max_levels=param->max_levels;
    SHORT i, level=0, status=SUCCESS;
    clock_t setup_start, setup_end;
    REAL setupduration;
    
#if DEBUG_MODE
    printf("### DEBUG: amg_setup_unsmoothP_unsmoothA_bsr ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if (print_level>PRINT_MOST)    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
    
    setup_start=clock();
    
    //each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(INT));
    
    for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
    
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    // each level stores the information of the strongly coupled neighborhoods
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); 
    
    mgl[0].near_kernel_dim   = 1;
    //mgl[0].near_kernel_basis = (double **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(double*));
    mgl[0].near_kernel_basis = NULL;
    
    // initialize ILU parameters
    mgl->ILU_levels = param->ILU_levels;
    ILU_param iluparam;
    if (param->ILU_levels>0) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }
    
    while ((mgl[level].A.ROW>param->coarse_dof) && (level<max_levels-1)) {

        /*-- setup ILU decomposition if necessary */
        if (level<param->ILU_levels) fasp_ilu_dbsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
    
        /*-- Aggregation --*/
        mgl[level].PP =  fasp_dbsr_getblk_dcsr(&mgl[level].A);  // use first block now, need to bechanged
        aggregation(&mgl[level].PP, &vertices[level], param, level+1, &Neighbor[level], &num_aggregations[level]);
        if (num_aggregations[level] * 4 > mgl[level].A.ROW) param->strong_coupled /=8.0; 
    
        /* -- Form Prolongation --*/      
        form_tentative_p_bsr(&vertices[level], &mgl[level].P, &mgl[0], level+1, num_aggregations[level]);
    
        /*-- Form resitriction --*/    
        fasp_dbsr_trans(&mgl[level].P, &mgl[level].R);
        
        /*-- Form coarse level stiffness matrix --*/    
        fasp_blas_dbsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
    
        fasp_dcsr_free(&Neighbor[level]);
        fasp_ivec_free(&vertices[level]);
    
        ++level;
    }
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(3*m*nb);    
    
    for (level=1; level<max_levels; ++level) {
        INT    m = mgl[level].A.ROW;
        mgl[level].num_levels = max_levels;     
        mgl[level].b = fasp_dvec_create(m*nb);
        mgl[level].x = fasp_dvec_create(m*nb);
        mgl[level].w = fasp_dvec_create(3*m*nb);    
    }
    
    mgl[max_levels-1].Ac = fasp_format_dbsr_dcsr(&mgl[max_levels-1].A);
    
#if With_UMFPACK    
    // Need to sort the matrix A for UMFPACK format
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].Ac, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran, &mgl[max_levels-1].Ac);
    fasp_dcsr_free(&Ac_tran);
#endif
    
    if (print_level>1) {
        double gridcom=0.0, opcom=0.0;
    
        printf("-----------------------------------------------\n");
        printf("  Level     Num of rows      Num of nonzeros\n");
        printf("-----------------------------------------------\n");
        for (level=0;level<max_levels;++level) {
            printf("%5d  %14d  %16d\n",level,mgl[level].A.ROW,mgl[level].A.NNZ);
            gridcom += mgl[level].A.ROW;
            opcom += mgl[level].A.NNZ;
        }
        printf("-----------------------------------------------\n");
    
        gridcom /= mgl[0].A.ROW;
        opcom /= mgl[0].A.NNZ;
        printf("(BSR) Unsmoothed Aggregation AMG grid complexity = %f\n", gridcom);
        printf("(BSR) Unsmoothed Aggregation AMG operator complexity = %f\n", opcom);
    }
    
    if (print_level>PRINT_NONE) {
        setup_end=clock();
        setupduration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
        print_cputime("(BSR) Unsmoothed Aggregation AMG setup",setupduration);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    
#if DEBUG_MODE
    printf("### DEBUG: amg_setup_unsmoothP_unsmoothA_bsr ...... [Finish]\n");
#endif
    
    return status;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
