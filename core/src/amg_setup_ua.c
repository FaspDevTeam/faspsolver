/*! \file amg_setup_ua.c
 *  \brief Unsmoothed Aggregation AMG: SETUP phase
 *
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

static SHORT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param);

static void aggregation(dCSRmat *A, ivector *vertices, AMG_param *param, INT levelNum, dCSRmat *N, INT *num_aggregations);
static void form_tentative_p(ivector *vertices, dCSRmat *tentp, AMG_data *mgl, INT levelNum, INT num_aggregations);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_amg_setup_ua(AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of unsmoothed aggregation AMG
 * 
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * \note Setup A, P, PT, levels using unsmoothed aggregation aggregation algorithm;
 *       Refer to Peter Vanek, Jan Madel and Marin Brezina
 *       Algebraic Multigrid on Unstructured Meshes, 1994
 * 
 * \author Xiaozhe Hu
 * \date 12/28/2011 
 *
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

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static SHORT amg_setup_unsmoothP_unsmoothA(AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of plain aggregation AMG, using unsmoothed P and unsmoothed A
 * 
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 * 
 * \author Xiaozhe Hu
 * \date 02/21/2011 
 */
static SHORT amg_setup_unsmoothP_unsmoothA (AMG_data *mgl, 
                                            AMG_param *param)
{
    const SHORT   cycle_type=param->cycle_type;
	const SHORT   print_level=param->print_level;
	const INT     m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;
	
    // local variables
	clock_t       setup_start, setup_end;
	REAL          setupduration;
    SHORT         level=0, status=SUCCESS, max_levels=param->max_levels;
	INT           i, j;
    
	setup_start = clock();
	
	if (cycle_type == AMLI_CYCLE) 
	{
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
	
#if DIAGONAL_PREF
    fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
#if DEBUG_MODE
		printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",
               level,mgl[level].A.row,mgl[level].A.nnz);
#endif
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
        
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif      
	}
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);	
	
	for (level=1; level<max_levels; ++level) {
		int	m = mgl[level].A.row;
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
        printf("Unsmoothed aggregation setup costs %f seconds.\n\n", setupduration);	
	}
	
	fasp_mem_free(vertices);
	fasp_mem_free(num_aggregations);
	fasp_mem_free(Neighbor);
    
	return status;
}

/**
 * \fn static void aggregation (dCSRmat *A, ivector *vertices, AMG_param *param, 
 *                              INT levelNum, dCSRmar *Neigh, INT *num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighborhoods 
 *
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
static void aggregation (dCSRmat *A, 
                         ivector *vertices, 
                         AMG_param *param, 
                         INT levelNum, 
                         dCSRmat *Neigh, 
                         INT *num_aggregations)
{
	// member of A
	INT row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
	INT *AIA = A->IA;
	INT *AJA = A->JA;
	REAL *Aval = A->val;
	
	// local variable
	REAL strongly_coupled  = param->strong_coupled;
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
	
	Neigh->JA = (int*)fasp_mem_realloc(Neigh->JA, (Neigh->IA[row])*sizeof(INT));
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
	INT *temp_C = (int*)fasp_mem_calloc(row,sizeof(INT));
	
	num_each_aggregation = (int*)fasp_mem_calloc(*num_aggregations,sizeof(INT));
	
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
 * \fn static void form_tentative_p (ivectors *vertices, dCSRmat *tentp, AMG_data *mgl, 
 *                                   INT levelNum, INT num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighborhoods 
 *
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
static void form_tentative_p (ivector *vertices, 
                              dCSRmat *tentp, 
                              AMG_data *mgl, 
                              INT levelNum, 
                              INT num_aggregations)
{
	INT i, j;
	REAL **basis = mgl->near_kernel_basis;
	
	/* Form tentative prolongation */
	tentp->row = vertices->row;
	tentp->col = num_aggregations;
	tentp->nnz = vertices->row;
	
	tentp->IA  = (int*)fasp_mem_calloc(tentp->row+1,sizeof(INT));	
	
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
	tentp->JA = (int*)fasp_mem_calloc(tentp->nnz, sizeof(INT));
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
