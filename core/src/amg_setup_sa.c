/*! \file amg_setup_sa.c
 *  \brief Smoothed Aggregation AMG: SETUP phase
 *
 *  TODO: Make unsmoothed aggregation standard alone and controlled with input? --Chensong
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

static INT amg_setup_smoothP_smoothA(AMG_data *mgl, AMG_param *param);
static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param);
static INT amg_setup_smoothP_unsmoothA(AMG_data *mgl, AMG_param *param);

static void aggregation(dCSRmat *A, ivector *vertices, AMG_param *param, INT levelNum, dCSRmat *N, INT *num_aggregations);
static void form_tentative_p(ivector *vertices, dCSRmat *tentp, AMG_data *mgl, INT levelNum, INT num_aggregations);
static void smooth_agg(dCSRmat *A, dCSRmat *tentp, dCSRmat *P, AMG_param *param, INT levelNum, dCSRmat *N);

#if FASP_USE_OPENMP
static INT amg_setup_unsmoothP_unsmoothA_omp(AMG_data *mgl,
                                             AMG_param *param,
                                             INT nthreads,
                                             INT openmp_holds);
static void aggregation_omp(dCSRmat *A,
                            ivector *vertices,
                            AMG_param *param,
                            INT levelNum,
                            dCSRmat *Neigh,
                            INT *num_aggregations,
                            INT nthreads,
                            INT openmp_holds);
static void form_tentative_p_omp(ivector *vertices,
                                 dCSRmat *tentp,
                                 AMG_data *mgl,
                                 INT levelNum,
                                 INT num_aggregations,
                                 INT nthreads,
                                 INT openmp_holds);
#endif // OMP

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_setup_sa(AMG_data *mgl, AMG_param *param)
 * \brief Set up phase of smoothed aggregation AMG
 * 
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 09/29/2009 
 *
 *  Modified by Chensong Zhang on 04/06/2010.
 *  Modified by Chensong Zhang on 05/09/2010.
 *  Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle
 */
INT fasp_amg_setup_sa (AMG_data *mgl, 
                       AMG_param *param)
{
	INT status=SUCCESS;
	
	INT type = 0; 
	if (param->tentative_smooth > SMALLREAL) type = 1;
	
	switch (type)
	{
		case 1: // Todo: Need to define types later -- Chensong
			status = amg_setup_smoothP_smoothA(mgl, param);
			break;
		case 2:
			status = amg_setup_smoothP_unsmoothA(mgl, param);
			break;
		default:
			status = amg_setup_unsmoothP_unsmoothA(mgl, param);
			break;	
	}
	
	return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 02/21/2011 
 */
static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
{
	const INT print_level=param->print_level;
	const INT m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
	
	INT max_levels=param->max_levels;
	INT i, j, level=0, status=SUCCESS;
	clock_t setup_start, setup_end;
	REAL setupduration;
	
#if DEBUG_MODE
	printf("fasp_amg_setup_sa ...... [Start]\n");
	printf("fasp_amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8)	printf("fasp_amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	setup_start=clock();
	
	if (param->cycle_type == AMLI_CYCLE) 
	{
		param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
		REAL lambda_max = 2.0;
		REAL lambda_min = lambda_max/4;
		fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
	}
	
	//each elvel stores the information of the number of aggregations
    INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(int));
	
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
		int	m = mgl[level].A.row;
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
		setup_end=clock();
		setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Unsmoothed Aggregation AMG setup costs %f seconds.\n", setupduration);	
	}
	
	fasp_mem_free(vertices);
	fasp_mem_free(num_aggregations);
	fasp_mem_free(Neighbor);
	
#if DEBUG_MODE
	printf("amg_setup_sa ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn static INT amg_setup_smoothP_smoothA(AMG_data *mgl, AMG_param *param)
 * \brief Set up phase of smoothed aggregation AMG, using smoothed P and smoothed A
 * 
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 02/21/2011 
 */
static INT amg_setup_smoothP_smoothA(AMG_data *mgl, AMG_param *param)
{
	const INT print_level=param->print_level;
	const INT m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
	
	INT max_levels=param->max_levels;
	INT i, j, level=0, status=SUCCESS;
	clock_t setup_start, setup_end;
	REAL setupduration;
	
#if DEBUG_MODE
	printf("amg_setup_sa ...... [Start]\n");
	printf("amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8)	printf("amg_setup: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	setup_start=clock();
	
	if (param->cycle_type == AMLI_CYCLE) 
	{
		param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
		REAL lambda_max = 2.0;
		REAL lambda_min = lambda_max/4;
		fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
	}
	
	INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(int)); //each elvel stores the information of the number of aggregations
	
	for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
	
	ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
	
	dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); // each level stores the information of the strongly coupled neighborhoods
	
	dCSRmat *tentp    = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); // each level stores the information of the tentative prolongations
	
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
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
		
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
	}
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);	
	
	for (level=1; level<max_levels; ++level) {
		int	m = mgl[level].A.row;
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
		printf("Smoothed Aggregation AMG grid complexity = %f\n", gridcom);
		printf("Smoothed Aggregation AMG operator complexity = %f\n", opcom);
	}
	
	if (print_level>0) {
		setup_end=clock();
		setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Smoothed Aggregation AMG setup costs %f seconds.\n", setupduration);	
	}
	
	fasp_mem_free(vertices);
	fasp_mem_free(num_aggregations);
	fasp_mem_free(Neighbor);
	fasp_mem_free(tentp);
	
#if DEBUG_MODE
	printf("amg_setup_sa ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn static INT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 02/21/2011 
 */
static INT amg_setup_smoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
{
	const INT print_level=param->print_level;
	const INT m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
	
	INT max_levels=param->max_levels;
	INT i, j, level=0, status=SUCCESS;
	clock_t setup_start, setup_end;
	REAL setupduration;
	
#if DEBUG_MODE
	printf("amg_setup_sa ...... [Start]\n");
	printf("amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8)	printf("amg_setup: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	setup_start=clock();
	
	if (param->cycle_type == AMLI_CYCLE) 
	{
		param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
		REAL lambda_max = 2.0;
		REAL lambda_min = lambda_max/4;
		fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
	}
	
	INT *num_aggregations = (INT *)fasp_mem_calloc(max_levels,sizeof(int)); //each elvel stores the information of the number of aggregations
	
	for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
	
	ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
	
	dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); // each level stores the information of the strongly coupled neighborhoods
	
	dCSRmat *tentp    = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); // each level stores the information of the tentative prolongations
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
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
		
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
		int	m = mgl[level].A.row;
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
		printf("Half smoothed Aggregation AMG grid complexity = %f\n", gridcom);
		printf("Half smoothed Aggregation AMG operator complexity = %f\n", opcom);
	}
	
	if (print_level>0) {
		setup_end=clock();
		setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Half smoothed Aggregation AMG setup costs %f seconds.\n", setupduration);	
	}
	
	fasp_mem_free(vertices);
	fasp_mem_free(num_aggregations);
	fasp_mem_free(Neighbor);
	fasp_mem_free(tentp);
	
#if DEBUG_MODE
	printf("amg_setup_sa ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn static void aggregation(dCSRmat *A, ivector *vertices, AMG_param *param, INT levelNum, dCSRmar *Neigh, INT *num_aggregations)
 * \brief Form aggregation based on strong coupled neighborhoods 
 * \param *A pointer to the coefficient matrices
 * \param *vertices pointer to the aggregation of vertics
 * \param *param pointer to AMG parameters
 * \param levelNum level number
 * \param *Neigh pointer to strongly coupled neighborhoods
 * \param *num_aggregations pointer to number of aggregations 
 * 
 * \author Xiaozhe Hu
 * \date 09/29/2009
 */
static void aggregation(dCSRmat *A, ivector *vertices, AMG_param *param, INT levelNum, dCSRmat *Neigh, INT *num_aggregations)
{
	// member of A
	INT row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
	INT *AIA =  A->IA;
	INT *AJA = A->JA;
	REAL *Aval = A->val;
	
	// local variable
	REAL strongly_coupled; 
	if (GE(param->tentative_smooth, SMALLREAL))
	{
		strongly_coupled= param->strong_coupled * pow(0.5, levelNum-1);
	}
	else
	{
		strongly_coupled= param->strong_coupled;
	}
	REAL strongly_coupled2 = pow(strongly_coupled,2);
	
	INT i,j,index, row_start, row_end;
	INT *NIA = NULL;
	INT *NJA = NULL;
	REAL *Nval = NULL;
	
	/*------------------------------------------*/
	/* Form strongly coupled neighborhood */
	/*------------------------------------------*/
	dvector diag; fasp_dcsr_getdiag(0, A, &diag);  // get the diagonal entries
	
	fasp_dcsr_alloc(row,col,nnz, Neigh);
	NIA =  Neigh->IA;
	NJA = Neigh->JA;
	Nval = Neigh->val;
    
	for(i=row+1; i--;) NIA[i] = AIA[i];
	
	index = 0;
	for (i=0; i<row; ++i) {
		NIA[i] = index;
		row_start = AIA[i]; row_end = AIA[i+1];
		for (j = row_start; j<row_end; ++j) {
			if ((AJA[j] == i) || (pow(Aval[j],2) >= strongly_coupled2 * fabs(diag.val[i]*diag.val[AJA[j]])))
			{
				NJA[index] = AJA[j];
				Nval[index] = Aval[j];
				index++;
			}
		}
	}
	
	NIA[row] = index;
	Neigh->nnz = index;
	
	Neigh->JA = (int*)fasp_mem_realloc(Neigh->JA, (Neigh->IA[row])*sizeof(int));
	Neigh->val = (REAL*)fasp_mem_realloc(Neigh->val, (Neigh->IA[row])*sizeof(REAL));
	
	NIA =  Neigh->IA;
	NJA = Neigh->JA;
	Nval = Neigh->val;
	
	fasp_dvec_free(&diag);
	
	/*------------------------------------------*/
	/* Initialization */
	/*------------------------------------------*/
	fasp_ivec_alloc(row, vertices);
	for (i=row;i--;) vertices->val[i] = -2;
	
	INT num_left = row;
	INT subset;
	INT max_aggregation = param->max_aggregation;
	INT *num_each_aggregation;
	INT count;
	
	*num_aggregations = 0;
	
	/*------------------------------------------*/
	/* Step 1. */
	/*------------------------------------------*/
	for (i=0; i<row; ++i){
		if ((AIA[i+1] - AIA[i] - 1) == 0){
			vertices->val[i] = -1;
			num_left--;
		}
		else{
			subset = 1;
			row_start = NIA[i]; row_end = NIA[i+1];
			for (j=row_start; j<row_end; ++j){
				if (vertices->val[NJA[j]] >= -1){
					subset = 0;
					break;
				}
			}
			
			if (subset == 1){
				count = 0;
				vertices->val[i] = *num_aggregations;
				num_left--;
				count++;
				row_start = NIA[i]; row_end = NIA[i+1];
				for (j=row_start; j<row_end;++j){
					if ((NJA[j]!=i) && (count < max_aggregation)){
						vertices->val[NJA[j]] = *num_aggregations;
						num_left--;
						count ++;
					}
				}
				(*num_aggregations)++;
			}
		}
	}
	
	/*------------------------------------------*/
	/* Step 2. */
	/*------------------------------------------*/
	INT *temp_C = (int*)fasp_mem_calloc(row,sizeof(int));
	
	num_each_aggregation = (int*)fasp_mem_calloc(*num_aggregations,sizeof(int));
	
	for (i=row;i--;){
		temp_C[i] = vertices->val[i];  
		if (vertices->val[i] >= 0){
			num_each_aggregation[vertices->val[i]] ++;
		}
	}
	
	for(i=0; i<row; ++i){
		if (vertices->val[i] < -1){
			row_start = NIA[i]; row_end = NIA[i+1];
			for (j=row_start;j<row_end;++j){
				if(temp_C[NJA[j]] >= -1 && num_each_aggregation[temp_C[NJA[j]]] < max_aggregation){
					vertices->val[i] = temp_C[NJA[j]];
					num_left--;
					num_each_aggregation[temp_C[NJA[j]]] ++ ;
					break;
				}
			}
		}
	}
	
	/*------------------------------------------*/
	/* Step 3. */
	/*------------------------------------------*/
	while (num_left > 0){
		for (i=0; i<row; ++i){
			if (vertices->val[i] < -1){
				count = 0;
				vertices->val[i] = *num_aggregations;
				num_left--;
				count++;
				row_start = NIA[i]; row_end = NIA[i+1];
				for (j=row_start; j<row_end;++j){
					if ((NJA[j]!=i) && (vertices->val[NJA[j]] < -1) && (count<max_aggregation)){
						vertices->val[NJA[j]] = *num_aggregations;
						num_left--;
						count++;
					}
				}
				(*num_aggregations)++;
			}
		}
	}
	
	fasp_mem_free(temp_C);
	fasp_mem_free(num_each_aggregation);
}

/**
 * \fn static void form_tentative_p(ivectors *vertices, dCSRmat *tentp, AMG_data *mgl, INT levelNum, INT num_aggregations)
 * \brief Form aggregation based on strong coupled neighborhoods 
 * \param *A pointer to the coefficient matrices
 * \param *vertices pointer to the aggregation of vertices
 * \param *P pointer to the prolongation operators 
 * \param *mgl pointer to AMG levele data
 * \param levelNum level number
 * \param num_aggregations number of aggregations
 *
 * \author Xiaozhe Hu
 * \date 09/29/2009
 */
static void form_tentative_p(ivector *vertices, dCSRmat *tentp, AMG_data *mgl, INT levelNum, INT num_aggregations)
{
	INT i, j;
	REAL **basis = mgl->near_kernel_basis;
	
	/* Form tentative prolongation */
	tentp->row = vertices->row;
	tentp->col = num_aggregations;
	tentp->nnz = vertices->row;
	
	tentp->IA  = (int*)fasp_mem_calloc(tentp->row+1,sizeof(int));	
	
	// local variables
	INT * IA = tentp->IA;
	INT *JA; 
	REAL *val; 
	INT *vval = vertices->val;
	
	const INT row = tentp->row;
	
	// first run
	for (i = 0, j = 0; i < row; i ++)
	{
		IA[i] = j;
		if (vval[i] > -1)
		{
			j ++;
		}
	}
	IA[row] = j;
	
	// allocate
	tentp->nnz = j;
	tentp->JA = (int*)fasp_mem_calloc(tentp->nnz, sizeof(int));
	tentp->val = (REAL*)fasp_mem_calloc(tentp->nnz, sizeof(REAL));
	
	JA = tentp->JA;
	val = tentp->val;
	
	// second run
	for (i = 0, j = 0; i < row; i ++)
	{
		IA[i] = j;
		if (vval[i] > -1)
		{
			JA[j] = vval[i];
			val[j] = basis[0][i];
			j ++;
		}
	}
}

/**
 * \fn static void smooth_agg(dCSRmat *A, dCSRmat *tentp, dCSRmat *P, AMG_param *param, INT levelNum, dCSRmat *N)
 * \brief Smooth the tentative prolongation
 * \param *A pointer to the coefficient matrices
 * \param *tentp pointer to the tentative prolongation operators
 * \param *P pointer to the prolongation operators 
 * \param *param pointer to AMG parameters
 * \param levelNum level number
 * \param *N pointer to strongly coupled neighborhoods
 *
 * \author Xiaozhe Hu
 * \date 09/29/2009
 */
static void smooth_agg(dCSRmat *A, dCSRmat *tentp, dCSRmat *P, AMG_param *param, INT levelNum, dCSRmat *N)
{
	INT i,j;
	
	REAL smooth_factor = param->tentative_smooth;
	
	INT row = A->row, col= A->col;
	
	INT filter = param->smooth_filter;
	
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
