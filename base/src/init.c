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
 * \fn void fasp_precond_data_null (precond_data *pcdata)
 *
 * \brief Initialize precond_data
 *
 * \param pcdata   Preconditioning data structure
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_precond_data_null (precond_data *pcdata)
{
    pcdata->AMG_type            = CLASSIC_AMG;
    pcdata->print_level         = PRINT_NONE;
    pcdata->maxit               = 500;
    pcdata->max_levels          = 20;
    pcdata->tol                 = 1e-8;
    pcdata->cycle_type          = V_CYCLE;
    pcdata->smoother            = SMOOTHER_GS;
    pcdata->smooth_order        = CF_ORDER;
    pcdata->presmooth_iter      = 1;
    pcdata->postsmooth_iter     = 1;
    pcdata->relaxation          = 1.1;
    pcdata->coarsening_type     = 1;
    pcdata->coarse_scaling      = ON;
    pcdata->amli_degree         = 1;
    pcdata->nl_amli_krylov_type = SOLVER_GCG;
}

/**
 * \fn AMG_data * fasp_amg_data_create (SHORT max_levels)
 *
 * \brief Create and initialize AMG_data for classical and SA AMG
 *
 * \param max_levels   Max number of levels allowed
 *
 * \return Pointer to the AMG_data data structure
 *
 * \author Chensong Zhang
 * \date   2010/04/06
 */
AMG_data * fasp_amg_data_create (SHORT max_levels)
{
    max_levels = MAX(1, max_levels); // at least allocate one level
    
    AMG_data *mgl = (AMG_data *)fasp_mem_calloc(max_levels,sizeof(AMG_data));
    
    INT i;
    for ( i=0; i<max_levels; ++i ) {
        mgl[i].max_levels = max_levels;
        mgl[i].num_levels = 0;
        mgl[i].near_kernel_dim = 0;
        mgl[i].near_kernel_basis = NULL;
        mgl[i].cycle_type = 0;
    }
    
    return(mgl);
}

/**
 * \fn AMG_data_bsr * fasp_amg_data_bsr_create (SHORT max_levels)
 *
 * \brief Create and initialize AMG_data data sturcture for AMG/SAMG (BSR format)
 *
 * \param max_levels   Max number of levels allowed
 *
 * \return Pointer to the AMG_data data structure
 *
 * \author Xiaozhe Hu
 * \date   08/07/2011
 */
AMG_data_bsr * fasp_amg_data_bsr_create (SHORT max_levels)
{
    SHORT i;
    
    AMG_data_bsr *mgl = (AMG_data_bsr *)fasp_mem_calloc(max_levels,sizeof(AMG_data_bsr));
    
    for (i=0; i<max_levels; ++i) {
        mgl[i].max_levels = max_levels;
        mgl[i].num_levels = 0;
        mgl[i].near_kernel_dim = 0;
        mgl[i].near_kernel_basis = NULL;
        mgl[i].A_nk = NULL;
        mgl[i].P_nk = NULL;
        mgl[i].R_nk = NULL;
    }
    
    return(mgl);
}

/**
 * \fn void fasp_ilu_data_alloc (INT iwk, INT nwork, ILU_data *iludata)
 *
 * \brief Allocate workspace for ILU factorization
 *
 * \param iwk       Size of the index array
 * \param nwork     Size of the work array
 * \param iludata   Pointer to the ILU_data
 *
 * \author Chensong Zhang
 * \date   2010/04/06
 */
void fasp_ilu_data_alloc (INT iwk,
                          INT nwork,
                          ILU_data *iludata)
{
#if DEBUG_MODE
    printf("### DEBUG: iwk=%d, nwork=%d\n", iwk, nwork);
#endif
    
    iludata->ijlu=(INT*)fasp_mem_calloc(iwk, sizeof(INT));
    
    iludata->luval=(REAL*)fasp_mem_calloc(iwk, sizeof(REAL));
    
    iludata->work=(REAL*)fasp_mem_calloc(nwork, sizeof(REAL));
    
    return;
}

/**
 * \fn void fasp_schwarz_data_free (Schwarz_data *schwarz)
 * \brief Free Schwarz_data data memeory space
 *
 * \param *schwarz  pointer to the AMG_data data
 *
 * \author Xiaozhe Hu
 * \date 2010/04/06
 */
void fasp_schwarz_data_free (Schwarz_data *schwarz)
{
    INT i;
	fasp_dcsr_free(&schwarz->A);
	
	for (i=0; i<schwarz->nblk; ++i) fasp_dcsr_free (&((schwarz->blk_data)[i]));

	schwarz->nblk = 0;
	fasp_mem_free (schwarz->iblock);
	fasp_mem_free (schwarz->jblock);
	//fasp_mem_free (schwarz->rhsloc);
	fasp_dvec_free (&schwarz->rhsloc1);
	fasp_dvec_free (&schwarz->xloc1);
	//fasp_mem_free (schwarz->au);
	//fasp_mem_free (schwarz->al);
	
	schwarz->memt = 0;
	fasp_mem_free (schwarz->mask);
	fasp_mem_free (schwarz->maxa);

#if WITH_MUMPS
	if (schwarz->mumps == NULL) return;
	else 
	for (i=0; i<schwarz->nblk; ++i) fasp_mumps_free (&((schwarz->mumps)[i]));
#endif
}

/**
 * \fn void fasp_amg_data_free (AMG_data *mgl, AMG_param *param)
 *
 * \brief Free AMG_data data memeory space
 *
 * \param mgl    Pointer to the AMG_data
 * \param param  Pointer to AMG parameters
 *
 * \author Chensong Zhang
 * \date   2010/04/06
 *
 * Modified by Chensong Zhang on 05/05/2013: Clean up param as well!
 */
void fasp_amg_data_free (AMG_data *mgl,
                         AMG_param *param)
{
    const INT max_levels = MAX(1,mgl[0].max_levels);
    
    INT i;
    
    for (i=0; i<max_levels; ++i) {
        if (&mgl[i].A) { fasp_dcsr_free(&mgl[i].A); }
        if (&mgl[i].P) { fasp_dcsr_free(&mgl[i].P); }
        if (&mgl[i].R) { fasp_dcsr_free(&mgl[i].R); }
        if (&mgl[i].b) { fasp_dvec_free(&mgl[i].b); }
        if (&mgl[i].x) { fasp_dvec_free(&mgl[i].x); }
        if (&mgl[i].w) { fasp_dvec_free(&mgl[i].w); }
        if (&mgl[i].cfmark) { fasp_ivec_free(&mgl[i].cfmark); }
        if (&mgl[i].LU) { fasp_ilu_data_free(&mgl[i].LU); }
        if (&mgl[i].schwarz) {fasp_schwarz_data_free (&mgl[i].schwarz);}
    }
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        fasp_mem_free(mgl->near_kernel_basis[i]);
        mgl->near_kernel_basis[i]=NULL;
    }
    
	// Clean direct solver data in necessary
	switch (param->coarse_solver) {

#if WITH_MUMPS   /* Destroy MUMPS direct solver on the coarsest level */
		case SOLVER_MUMPS: {
            mgl[max_levels-1].mumps.job = 3;
			fasp_solver_mumps_steps(&mgl[max_levels-1].A, &mgl[max_levels-1].b, &mgl[max_levels-1].x, &mgl[max_levels-1].mumps);
			break;
		}
#endif

		default: // Do nothing!
			break;
	}
    
    fasp_mem_free(mgl->near_kernel_basis); mgl->near_kernel_basis = NULL;
    
    fasp_mem_free(mgl); mgl = NULL;
    
    if (param != NULL) {
        if ( param->cycle_type == AMLI_CYCLE ) fasp_mem_free(param->amli_coef);
    }
}

/**
 * \fn void fasp_amg_data_bsr_free (AMG_data_bsr *mgl)
 *
 * \brief Free AMG_data_bsr data memeory space
 *
 * \param mgl  Pointer to the AMG_data_bsr
 *
 * \author Xiaozhe Hu
 * \date   2013/02/13
 */
void fasp_amg_data_bsr_free (AMG_data_bsr *mgl)
{
    const INT max_levels = mgl[0].max_levels;
    INT i;
    
    for (i=0; i<max_levels; ++i) {
        if (&mgl[i].A) { fasp_dbsr_free(&mgl[i].A); }
        if (&mgl[i].P) { fasp_dbsr_free(&mgl[i].P); }
        if (&mgl[i].R) { fasp_dbsr_free(&mgl[i].R); }
        if (&mgl[i].b) { fasp_dvec_free(&mgl[i].b); }
        if (&mgl[i].x) { fasp_dvec_free(&mgl[i].x); }
        if (&mgl[i].diaginv) { fasp_dvec_free(&mgl[i].diaginv); }
        if (&mgl[i].Ac) { fasp_dcsr_free(&mgl[i].Ac); }
        if (&mgl[i].PP) { fasp_dcsr_free(&mgl[i].PP); }
        if (&mgl[i].pw) { fasp_mem_free(mgl[i].pw); }
        if (&mgl[i].SS) { fasp_dbsr_free(&mgl[i].SS); }
        if (&mgl[i].sw) { fasp_mem_free(mgl[i].sw); }
        if (&mgl[i].diaginv_SS) { fasp_dvec_free(&mgl[i].diaginv_SS); }
        if (&mgl[i].PP_LU) { fasp_ilu_data_free(&mgl[i].PP_LU); }
        if (&mgl[i].w) { fasp_dvec_free(&mgl[i].w); }
        if (&mgl[i].cfmark) { fasp_ivec_free(&mgl[i].cfmark); }
        if (&mgl[i].LU) { fasp_ilu_data_free(&mgl[i].LU); }
        //if (&mgl[i].schwarz) {fasp_schwarz_data_free (&mgl[i].schwarz);}
    }
    
    for (i=0; i<mgl->near_kernel_dim; ++i) {
        if (&mgl->near_kernel_basis[i]) fasp_mem_free(mgl->near_kernel_basis[i]);
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
 * \param ILUdata   Pointer to ILU_data
 *
 * \author Chensong Zhang
 * \date   2010/04/03
 */
void fasp_ilu_data_free (ILU_data *ILUdata)
{
    if (ILUdata==NULL) return;
    
    fasp_mem_free(ILUdata->ijlu);  ILUdata->ijlu  = NULL;
    fasp_mem_free(ILUdata->luval); ILUdata->luval = NULL;
    fasp_mem_free(ILUdata->work);  ILUdata->work  = NULL;
    
    ILUdata->row = ILUdata->col = ILUdata->nzlu = ILUdata ->nwork = ILUdata->nb= 0;
}

/**
 * \fn void fasp_ilu_data_null (ILU_data *ILUdata)
 *
 * \brief Initialize ILU data
 *
 * \param ILUdata   Pointer to ILU_data
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_ilu_data_null (ILU_data *ILUdata)
{
    ILUdata->row = ILUdata->col = ILUdata->nzlu = 0;
    ILUdata->ijlu = NULL; ILUdata->luval = NULL;
}

/**
 * \fn void fasp_precond_null (precond *pcdata)
 *
 * \brief Initialize precond data
 *
 * \param pcdata   Pointer to precond
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_precond_null (precond *pcdata)
{
    pcdata->data = NULL;
    pcdata->fct  = NULL;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
