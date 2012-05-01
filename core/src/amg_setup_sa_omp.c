/*! \file amg_setup_sa_omp.c
 *  \brief Smoothed Aggregation AMG: SETUP phase
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"


#if FASP_USE_OPENMP
static int amg_setup_unsmoothP_unsmoothA_omp(AMG_data *mgl,
                                             AMG_param *param,
                                             int nthreads,
                                             int openmp_holds);
static void aggregation_omp(dCSRmat *A,
                            ivector *vertices,
                            AMG_param *param,
                            int levelNum,
                            dCSRmat *Neigh,
                            int *num_aggregations,
                            int nthreads,
                            int openmp_holds);
static void form_tentative_p_omp(ivector *vertices,
                                 dCSRmat *tentp,
                                 AMG_data *mgl,
                                 int levelNum,
                                 int num_aggregations,
                                 int nthreads,
                                 int openmp_holds);
#endif // OMP

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*---------------------------------omp----------------------------------------*/

/**
 * \fn int fasp_amg_setup_sa_omp(AMG_data *mgl, AMG_param *param, int nthreads, int openmp_holds)
 * \brief Set up phase of smoothed aggregation AMG
 * 
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
int fasp_amg_setup_sa_omp (AMG_data *mgl, AMG_param *param, int nthreads, int openmp_holds)
{
    int status=SUCCESS;
#if FASP_USE_OPENMP
    int type = 0; 
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
    status = amg_setup_unsmoothP_unsmoothA_omp(mgl, param, nthreads, openmp_holds);
    break;
    }
#endif
    return status;
}

#if FASP_USE_OPENMP

/**
 * \fn static int amg_setup_unsmoothP_unsmoothA_omp(AMG_data *mgl, AMG_param *param, int nthreads, int openmp_holds)
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param mgl     pointer to AMG_data data
 * \param param   pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * Setup A, P, PT, levels using smoothed aggregation
 * concrete algorithm see paper
 * Peter Vanek, Jan Madel and Marin Brezina, Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/15/2011
 */
static int amg_setup_unsmoothP_unsmoothA_omp(AMG_data *mgl,
                                             AMG_param *param,
                                             int nthreads,
                                             int openmp_holds)
{
    int status=SUCCESS;
#if FASP_USE_OPENMP

    const int print_level=param->print_level;
    const int m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
    
    int max_levels=param->max_levels;
    int i, level=0;
    double setup_start, setup_end;
    double setupduration;
    
#if DEBUG_MODE
    printf("fasp_amg_setup_sa ...... [Start]\n");
    printf("fasp_amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if (print_level>8) printf("fasp_amg_setup_sa: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
    
    setup_start=omp_get_wtime();
    
    if (param->cycle_type == AMLI_CYCLE) 
    {
    param->amli_coef = (double *)fasp_mem_calloc(param->amli_degree+1,sizeof(double));
    double lambda_max = 2.0;
    double lambda_min = lambda_max/4;
    fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
    }
    
    int *num_aggregations = (int *)fasp_mem_calloc(max_levels,sizeof(int)); //each elvel stores the information of the number of aggregations
    
    for (i=0; i<max_levels; ++i) num_aggregations[i] = 0;
    
    ivector *vertices = (ivector *)fasp_mem_calloc(max_levels,sizeof(ivector));
    
    dCSRmat *Neighbor = (dCSRmat *)fasp_mem_calloc(max_levels,sizeof(dCSRmat)); // each level stores the information of the strongly coupled neighborhoods
    
    mgl[0].near_kernel_dim   = 1;
    mgl[0].near_kernel_basis = (double **)fasp_mem_calloc(mgl->near_kernel_dim,sizeof(double*));
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
    mgl[0].near_kernel_basis[i] = (double *)fasp_mem_calloc(m,sizeof(double));
    fasp_array_set_omp (m, mgl[0].near_kernel_basis[i], 1.0, nthreads, openmp_holds);
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
    aggregation_omp(&mgl[level].A, &vertices[level], param, level+1, &Neighbor[level], &num_aggregations[level], nthreads, openmp_holds);
    if (num_aggregations[level] * 4 > mgl[level].A.row) param->strong_coupled /=2.0; 
    
    /* -- Form Prolongation --*/      
    form_tentative_p_omp(&vertices[level], &mgl[level].P, &mgl[0], level+1, num_aggregations[level], nthreads, openmp_holds);
    
    /*-- Form resitriction --*/    
    fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
    
    /*-- Form coarse level stiffness matrix --*/    
#if 0
    fasp_blas_dcsr_rap_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A,nthreads,openmp_holds);
#else    
    fasp_blas_dcsr_rap_agg_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A,nthreads,openmp_holds);
#endif    
    fasp_dcsr_free(&Neighbor[level]);
    fasp_ivec_free(&vertices[level]);
    
    ++level;
    }
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(m);    
    
    for (level=1; level<max_levels; ++level) {
    int    m = mgl[level].A.row;
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
    double gridcom=0.0, opcom=0.0;
    
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
    setup_end=omp_get_wtime();
    setupduration = setup_end - setup_start;
    print_cputime("Unsmoothed Aggregation AMG setup",setupduration);
    }
    
    fasp_mem_free(vertices);
    fasp_mem_free(num_aggregations);
    fasp_mem_free(Neighbor);
    
#if DEBUG_MODE
    printf("amg_setup_sa ...... [Finish]\n");
#endif

#endif // OMP
    return status;
}

/**
 * \fn static void aggregation_omp(dCSRmat *A, ivector *vertices, AMG_param *param, int levelNum, dCSRmar *Neigh, int *num_aggregations, int nthreads, int openmp_holds)
 * \brief Form aggregation based on strong coupled neighborhoods 
 * \param A pointer to the coefficient matrices
 * \param vertices pointer to the aggregation of vertics
 * \param param pointer to AMG parameters
 * \param levelNum level number
 * \param Neigh pointer to strongly coupled neighborhoods
 * \param num_aggregations pointer to number of aggregations 
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
static void aggregation_omp(dCSRmat *A,
                            ivector *vertices,
                            AMG_param *param,
                            int levelNum,
                            dCSRmat *Neigh,
                            int *num_aggregations,
                            int nthreads,
                            int openmp_holds)
{
#if FASP_USE_OPENMP
    double strongly_coupled;
    
    //strongly_coupled= param->strong_coupled * pow(0.5, levelNum-1);
    
    if (param->tentative_smooth > SMALLREAL)
    {
    strongly_coupled= param->strong_coupled * pow(0.5, levelNum-1);
    }
    else
    {
    strongly_coupled= param->strong_coupled;
    }
    
    int i,j,index;
    int row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
    
    /* Form strongly coupled neighborhood */
    dvector diag; fasp_dcsr_getdiag_omp(0, A, &diag, nthreads, openmp_holds);  // get the diagonal entries
    
    fasp_dcsr_alloc(row,col,nnz, Neigh);
    
    //for(i=0; i<=row; ++i) Neigh->IA[i] = A->IA[i];
    fasp_iarray_cp_omp (row+1, A->IA, Neigh->IA, nthreads, openmp_holds);
    
    index = 0;
    for (i=0; i<row; ++i) {
    Neigh->IA[i] = index;
    for (j = A->IA[i]; j<A->IA[i+1]; ++j) {
    if (A->JA[j] == i) {
    Neigh->JA[index] = i;
    Neigh->val[index] = A->val[j];
    index++;
    }
    else if (ABS(A->val[j]) >= strongly_coupled * sqrt(fabs(diag.val[i]*diag.val[A->JA[j]]))){
    Neigh->JA[index] = A->JA[j];
    Neigh->val[index] = A->val[j];
    index++;
    }
    }
    }
    
    Neigh->IA[row] = index;
    Neigh->nnz = index;
    
    Neigh->JA = (int*)fasp_mem_realloc(Neigh->JA, (Neigh->IA[row])*sizeof(int));
    Neigh->val = (double*)fasp_mem_realloc(Neigh->val, (Neigh->IA[row])*sizeof(double));
    
    fasp_dvec_free(&diag);
    
    /* Initialization */
    fasp_ivec_alloc(row, vertices);
    fasp_iarray_set_omp(row, vertices->val, -2, nthreads, openmp_holds);
    
    int num_left = row;
    int subset;
    int max_aggregation = param->max_aggregation;
    int *num_each_aggregation;
    int count;
    
    *num_aggregations = 0;
    
    /* Step 1. */
    for (i=0; i<row; ++i){
    if (A->IA[i+1] - A->IA[i] == 1){
    vertices->val[i] = -1;
    num_left--;
    }
    else{
    subset = 1;
    for (j=Neigh->IA[i]; j<Neigh->IA[i+1]; ++j){
    if (vertices->val[Neigh->JA[j]] >= -1){
    subset = 0;
    break;
    }
    }
    
    if (subset == 1){
    count = 0;
    vertices->val[i] = *num_aggregations;
    num_left--;
    count++;
    for (j=Neigh->IA[i]; j<Neigh->IA[i+1];++j){
    if ((Neigh->JA[j]!=i) && (count < max_aggregation)){
    vertices->val[Neigh->JA[j]] = *num_aggregations;
    num_left--;
    count ++;
    }
    }
    (*num_aggregations)++;
    }
    }
    }
    
    /* Step 2. */
    int *temp_C = (int*)fasp_mem_calloc(row,sizeof(int));
    
    num_each_aggregation = (int*)fasp_mem_calloc(*num_aggregations,sizeof(int));
    
    for (i=0;i<row;++i){
    temp_C[i] = vertices->val[i];
    if (vertices->val[i] >= 0){
    num_each_aggregation[vertices->val[i]] ++;
    }
    }
    
    for(i=0; i<row; ++i){
    if (vertices->val[i] < -1){
    for (j=Neigh->IA[i];j<Neigh->IA[i+1];++j){
    if(temp_C[Neigh->JA[j]] >= -1 && num_each_aggregation[temp_C[Neigh->JA[j]]] < max_aggregation){
    vertices->val[i] = temp_C[Neigh->JA[j]];
    num_left--;
    num_each_aggregation[temp_C[Neigh->JA[j]]] ++ ;
    break;
    }
    }
    }
    }
    
    /* Step 3. */
    while (num_left > 0){
    for (i=0; i<row; ++i){
    if (vertices->val[i] < -1){
    count = 0;
    vertices->val[i] = *num_aggregations;
    num_left--;
    count++;
    for (j=Neigh->IA[i]; j<Neigh->IA[i+1];++j){
    if ((Neigh->JA[j]!=i) && (vertices->val[Neigh->JA[j]] < -1) && (count<max_aggregation)){
    vertices->val[Neigh->JA[j]] = *num_aggregations;
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
#endif // OMP

}

/**
 * \fn static void form_tentative_p(ivectors *vertices, dCSRmat *tentp, AMG_data *mgl, int levelNum, int num_aggregations, int nthreads, int openmp_holds)
 * \brief Form aggregation based on strong coupled neighborhoods 
 * \param A pointer to the coefficient matrices
 * \param vertices pointer to the aggregation of vertices
 * \param P pointer to the prolongation operators 
 * \param mgl pointer to AMG levele data
 * \param levelNum level number
 * \param num_aggregations number of aggregations
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
static void form_tentative_p_omp(ivector *vertices,
                                 dCSRmat *tentp,
                                 AMG_data *mgl,
                                 int levelNum,
                                 int num_aggregations,
                                 int nthreads,
                                 int openmp_holds)
{
#if FASP_USE_OPENMP
    int i, myid, mybegin, myend, row_plus_one;
    double **basis = mgl->near_kernel_basis;
    
    /* Form tentative prolongation */
    tentp->row = vertices->row;
    tentp->col = num_aggregations;
    tentp->nnz = vertices->row;
    row_plus_one = tentp->row + 1;
    
    tentp->IA  = (int*)fasp_mem_calloc(row_plus_one,sizeof(int));
    tentp->JA  = (int*)fasp_mem_calloc(tentp->nnz,sizeof(int));
    tentp->val = (double*)fasp_mem_calloc(tentp->nnz,sizeof(double));
    
    if (row_plus_one > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, row_plus_one, mybegin, myend);
    for (i=mybegin; i<myend; ++i) {
    tentp->IA[i] = i;
    }
    }
    }
    else {
    for(i=0; i<row_plus_one;++i) {
    tentp->IA[i] = i;
    }
    }
    
    if (tentp->nnz > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i) ////num_threads(nthreads)
    {
    myid = omp_get_thread_num();
    FASP_GET_START_END(myid, nthreads, tentp->nnz, mybegin, myend);
    for (i=mybegin; i<myend; ++i) {
    if (vertices->val[i] == -1) {
    tentp->JA[i] = 0;
    tentp->val[i] = 0;
    }
    else {
    tentp->JA[i] = vertices->val[i];
    tentp->val[i] = basis[0][i];
    }
    }
    }
    }
    else {
    for(i=0;i<tentp->nnz;++i) {
    if (vertices->val[i] == -1){
    tentp->JA[i] = 0;
    tentp->val[i] = 0;
    }
    else{
    tentp->JA[i] = vertices->val[i];
    tentp->val[i] = basis[0][i];
    }
    }
    }

#endif
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
