/*! \file precond_omp.c
 *  \brief Preconditioners
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*---------------------------------omp----------------------------------------*/

/**
 * \fn void fasp_precond_amg_omp(double *r, double *z, void *data, int nthreads, int openmp_holds)
 * \brief get z from r by classic AMG
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
void fasp_precond_amg_omp(double *r, double *z, void *data, int nthreads, int openmp_holds)
{
#if FASP_USE_OPENMP
	precond_data *predata=(precond_data *)data;
	const int m=predata->mgl_data[0].A.row;
	const int maxit=predata->maxit;
	unsigned int i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
	amgparam.cycle_type = predata->cycle_type;
	amgparam.smoother   = predata->smoother;
	amgparam.smooth_order = predata->smooth_order;
	amgparam.presmooth_iter  = predata->presmooth_iter;
	amgparam.postsmooth_iter = predata->postsmooth_iter;
	amgparam.relaxation = predata->relaxation;
	amgparam.coarse_scaling = predata->coarse_scaling;
	amgparam.tentative_smooth = predata->tentative_smooth;
	amgparam.ILU_levels = predata->mgl_data->ILU_levels;
	
	AMG_data *mgl = predata->mgl_data;
	mgl->b.row=m; fasp_array_cp_omp(m,r,mgl->b.val,nthreads,openmp_holds); // residual is an input 
	mgl->x.row=m; fasp_dvec_set_omp(m,&mgl->x,0.0,nthreads,openmp_holds);
	
	for (i=0;i<maxit;++i) fasp_solver_mgcycle_omp(mgl,&amgparam,nthreads,openmp_holds);
	//fasp_solver_mgrecurmgl,&amgparam,0); //
	
	fasp_array_cp_omp(m,mgl->x.val,z,nthreads,openmp_holds);
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
