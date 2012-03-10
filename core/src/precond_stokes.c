/*! \file precond_stokes.c
 *  \brief Preconditioners for Stokes-type problems
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_stokes_bdiag (REAL *r, REAL *z, void *data)
 * \brief block diagonal preconditioning for Stokes equation
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 */
void fasp_precond_stokes_bdiag (REAL *r, 
                                REAL *z, 
                                void *data)
{
	precond_Stokes_data *predata=(precond_Stokes_data *)data;
	
	const INT col = predata->col, colA = predata->colA, colB = predata->colB;
	const INT maxit = predata->maxit;
	REAL *diagptr=predata->diag_M->val;
	
	// local variables
	REAL	*tempr = predata->w;		
	INT i;
	
	//! prepare	 AMG preconditioner 
	AMG_data *mgl = predata->mgl_data;
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
	amgparam.cycle_type = predata->cycle_type;
	amgparam.smoother   = predata->smoother;
	amgparam.presmooth_iter  = predata->presmooth_iter;
	amgparam.postsmooth_iter = predata->postsmooth_iter;
	amgparam.relaxation      = predata->relaxation;
	amgparam.coarse_scaling  = predata->coarse_scaling;
	amgparam.ILU_levels      = predata->mgl_data->ILU_levels;		
	
	//! back up r, setup z;
	fasp_array_cp(col, r, tempr);
	fasp_array_set(col, z, 0.0);
	
	//! Solve A by AMG
	mgl->b.row=colA; fasp_array_cp(colA,r,mgl->b.val); // residual is an input 
	mgl->x.row=colA; fasp_dvec_set(colA,&mgl->x,0.0);	
	for (i=0;i<maxit;++i) fasp_solver_mgcycle(mgl,&amgparam); //fasp_solver_mgrecur(mgl,&amgparam,0); //
	
	//! Solve M by PCG
	for (i=0;i<colB;++i) {
		if (ABS(diagptr[i])>SMALLREAL) z[colA+i]/=diagptr[i];
	}		
	
	//! restore r
	fasp_array_cp(col, tempr, r);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
