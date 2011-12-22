/*! \file init.c
 *  \brief Initialize important data structures
 *
 *  \note Every structures should be initialized before usage.
 *
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn AMG_data * fasp_amg_data_create (INT max_levels)
 *
 * \brief Create and initialize AMG_data for classical and SA AMG
 *
 * \param   max_levels   max number of levels allowed
 *
 * \return *mgl          pointer to the AMG_data data structure
 *
 * \author Chensong Zhang
 * \date 2010/04/06
 */
AMG_data * fasp_amg_data_create (INT max_levels)
{			
	AMG_data *mgl = (AMG_data *)fasp_mem_calloc(max_levels,sizeof(AMG_data));
	
	INT i;
	for (i=0; i<max_levels; ++i) {
		mgl[i].max_levels = max_levels;
		mgl[i].num_levels = 0;
		mgl[i].near_kernel_dim = 0;
		mgl[i].near_kernel_basis = NULL;
	}
	
	return(mgl);
}

/**
 * \fn void fasp_ilu_data_alloc (INT iwk, INT nwork, ILU_data *iludata)
 *
 * \brief Allocate workspace for ILU factorization
 *
 * \param iwk       integer, length of the index array
 * \param nwork     integer, length of the work array
 * \param *iludata  pointer, the ILU facotrization
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
void fasp_ilu_data_alloc (INT iwk, 
                          INT nwork, 
                          ILU_data *iludata)
{	
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_data_alloc: iwk=%d, nwork=%d\n", iwk, nwork);
#endif
	
	iludata->ijlu=(int*)fasp_mem_calloc(iwk, sizeof(int));
	
#if CHMEM_MODE		
	total_alloc_mem += iwk*sizeof(int);
#endif
	
	iludata->luval=(REAL*)fasp_mem_calloc(iwk, sizeof(REAL)); 
	
#if CHMEM_MODE		
	total_alloc_mem += iwk*sizeof(REAL);
#endif
	
	iludata->work=(REAL*)fasp_mem_calloc(nwork, sizeof(REAL)); 
	
#if CHMEM_MODE		
	total_alloc_mem += nwork*sizeof(REAL);
#endif
	
	return;
}

/**
 * \fn void fasp_amg_data_free (AMG_data *mgl)
 * \brief Free AMG_data data memeory space
 *
 * \param *mgl  pointer to the AMG_data data
 *
 * \author Chensong Zhang
 * \date 2010/04/06 
 */
void fasp_amg_data_free (AMG_data *mgl)
{		
	const INT max_levels = mgl[0].max_levels;
	unsigned INT i;
	
	for (i=0; i<max_levels; ++i) {
		if (&mgl[i].A) { fasp_dcsr_free(&mgl[i].A); }
		if (&mgl[i].P) { fasp_dcsr_free(&mgl[i].P); }
		if (&mgl[i].R) { fasp_dcsr_free(&mgl[i].R); }
		if (&mgl[i].b) { fasp_dvec_free(&mgl[i].b); }
		if (&mgl[i].x) { fasp_dvec_free(&mgl[i].x); }
		if (&mgl[i].w) { fasp_dvec_free(&mgl[i].w); }
		if (&mgl[i].cfmark) { fasp_ivec_free(&mgl[i].cfmark); }
		if (&mgl[i].LU) { fasp_ilu_data_free(&mgl[i].LU); }
	}
	
	for (i=0; i<mgl->near_kernel_dim; ++i) {
		fasp_mem_free(mgl->near_kernel_basis[i]); 
		mgl->near_kernel_basis[i]=NULL;
	}
	fasp_mem_free(mgl->near_kernel_basis); mgl->near_kernel_basis = NULL;
	fasp_mem_free(mgl); mgl = NULL;
}

/**
 * \fn void fasp_ilu_data_free (ILU_data *ILUdata)
 *
 * \brief Create ILU_data sturcture
 *
 * \param nnz       number of nonzeros
 * \param ILUdata   pointer to ILU_data
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_ilu_data_free(ILU_data *ILUdata)
{		
	if (ILUdata==NULL) return;
	
	fasp_mem_free(ILUdata->ijlu);  ILUdata->ijlu  = NULL;
	fasp_mem_free(ILUdata->luval); ILUdata->luval = NULL;
	fasp_mem_free(ILUdata->work);  ILUdata->work  = NULL;
	
	ILUdata->row = ILUdata->col = ILUdata->nzlu = ILUdata ->nwork = ILUdata->nb= 0;
}

/**
 * \fn void fasp_ilu_data_init (ILU_data *ILUdata)
 *
 * \brief Initialize ILU data
 *
 * \param Input   pointer to ILU_data
 *
 * \author Chensong Zhang
 * \date 2010/03/23 
 */
void fasp_ilu_data_init (ILU_data *ILUdata)
{
	ILUdata->row = ILUdata->col = ILUdata->nzlu = 0;
	ILUdata->ijlu = NULL; ILUdata->luval = NULL;
}

/**
 * \fn void fasp_precond_init (precond *pdata)
 *
 * \brief Initialize precond data
 *
 * \param Input   pointer to precond
 *
 * \author Chensong Zhang
 * \date 2010/03/23 
 */
void fasp_precond_init (precond *pdata)
{
	pdata->data = NULL;
	pdata->fct  = NULL;
}

/**
 * \fn void fasp_precond_data_init (precond_data *pdata)
 *
 * \brief Initialize precond_data
 *
 * \param Input   pointer to precond_data
 *
 * \author Chensong Zhang
 * \date 2010/03/23 
 */
void fasp_precond_data_init (precond_data *pdata)
{
	pdata->print_level     = 0;
	pdata->max_iter        = 500;
	pdata->max_levels      = 12;
	pdata->tol             = 1e-8;
	pdata->cycle_type      = V_CYCLE;
	pdata->smoother        = GS;
	pdata->coarse_scaling  = OFF;
	pdata->presmooth_iter  = 2;
	pdata->postsmooth_iter = 2;
	pdata->coarsening_type = 1;
	pdata->relaxation      = 0.9;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
