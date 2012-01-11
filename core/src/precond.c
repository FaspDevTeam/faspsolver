/*! \file precond.c
 *  \brief Preconditioners
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_diag(double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Chensong Zhang
 * \date 04/06/2010
 */
void fasp_precond_diag (double *r, 
                        double *z, 
                        void *data)
{
	dvector *diag=(dvector *)data;
	double *diagptr=diag->val;
	unsigned int i, m=diag->row;	
	
	memcpy(z,r,m*sizeof(double));
	for (i=0;i<m;++i) {
		if (ABS(diag->val[i])>SMALLREAL) z[i]/=diagptr[i];
	}	
}

/**
 * \fn void fasp_precond_ilu(double *r, double *z, void *data)
 * \brief preconditioning using ILU decomposition
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Shiquan Zhang
 * \date 04/06/2010
 */
void fasp_precond_ilu (double *r, 
                       double *z, 
                       void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const int m=iludata->row, mm1=m-1, memneed=2*m;
	double *zz, *zr;
	
	if (iludata->nwork<memneed) goto MEMERR; // check this outside this subroutine!!
	
	zz = iludata->work; 
	zr = iludata->work+m;
	fasp_array_cp(m, r, zr); 	
	
	{
		int i, j, jj, begin_row, end_row, mm2=m-2;
		int *ijlu=iludata->ijlu;
		double *lu=iludata->luval;
		
		// forward sweep: solve unit lower matrix equation L*zz=zr
		zz[0]=zr[0];
		for (i=1;i<=mm1;++i) {
			begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
			for (j=begin_row;j<=end_row;++j) {
				jj=ijlu[j];
				if (jj<i) zr[i]-=lu[j]*zz[jj];
				else break;
			}
			zz[i]=zr[i];
		}
		
		// backward sweep: solve upper matrix equation U*z=zz
		z[mm1]=zz[mm1]*lu[mm1];
		for (i=mm2;i>=0;i--) {
			begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
			for (j=end_row;j>=begin_row;j--) {
				jj=ijlu[j];
				if (jj>i) zz[i]-=lu[j]*z[jj];
				else break;
			} 
			z[i]=zz[i]*lu[i];
		}
	}
	
	return;
	
MEMERR:
	printf("Error: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_ilu_forward(double *r, double *z, void *data)
 * \brief preconditioning using ILU decomposition: only forwear sweep
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu, Shiquang Zhang
 * \date 04/06/2010
 */
void fasp_precond_ilu_forward (double *r, 
                               double *z, 
                               void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const int m=iludata->row, mm1=m-1, memneed=2*m;
	double *zz, *zr;
	
	if (iludata->nwork<memneed) goto MEMERR; 
	
	zz = iludata->work; 
	zr = iludata->work+m;
	fasp_array_cp(m, r, zr); 	
	
	{
		int i, j, jj, begin_row, end_row;
		int *ijlu=iludata->ijlu;
		double *lu=iludata->luval;
		
		// forward sweep: solve unit lower matrix equation L*z=r
		zz[0]=zr[0];
		for (i=1;i<=mm1;++i) {
			begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
			for (j=begin_row;j<=end_row;++j) {
				jj=ijlu[j];
				if (jj<i) zr[i]-=lu[j]*zz[jj];
				else break;
			}
			zz[i]=zr[i];
		}
	}
	
	fasp_array_cp(m, zz, z); 
	
	return;
	
MEMERR:
	printf("Error: Need %d memory, only %d available!!!", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_ilu_backward (double *r, double *z, void *data)
 * \brief preconditioning using ILU decomposition: only backward sweep
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu, Shiquan  Zhang
 * \date 04/06/2010
 */
void fasp_precond_ilu_backward (double *r, 
                                double *z, 
                                void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const int m=iludata->row, mm1=m-1, memneed=2*m;
	double *zz;//, *zr;
	
	if (iludata->nwork<memneed) goto MEMERR; 
	
	zz = iludata->work; 
	//zr = iludata->work+m;
	fasp_array_cp(m, r, zz); 	
	
	{
		int i, j, jj, begin_row, end_row, mm2=m-2;
		int *ijlu=iludata->ijlu;
		double *lu=iludata->luval;
		
		// backward sweep: solve upper matrix equation U*z=zz
		z[mm1]=zz[mm1]*lu[mm1];
		for (i=mm2;i>=0;i--) {
			begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
			for (j=end_row;j>=begin_row;j--) {
				jj=ijlu[j];
				if (jj>i) zz[i]-=lu[j]*z[jj];
				else break;
			} 
			z[i]=zz[i]*lu[i];
		}
		
	}
	
	return;
	
MEMERR:
	printf("Error: Need %d memory, only %d available!!!", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_amg(double *r, double *z, void *data)
 * \brief get z from r by classic AMG
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Chensong Zhang
 * \date 04/06/2010
 */
void fasp_precond_amg (double *r, 
                       double *z, 
                       void *data)
{
	precond_data *predata=(precond_data *)data;
	const int m=predata->mgl_data[0].A.row;
	const int maxit=predata->maxit;
	unsigned int i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    amgparam.AMG_type = predata->AMG_type;
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
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_mgcycle(mgl,&amgparam); //fasp_solver_mgrecur(mgl,&amgparam,0);
	
	fasp_array_cp(m,mgl->x.val,z);	
}

/**
 * \fn void fasp_precond_famg(double *r, double *z, void *data)
 * \brief get z from r by Full AMG
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 02/27/2011
 */
void fasp_precond_famg (double *r, 
                        double *z, 
                        void *data)
{
	precond_data *predata=(precond_data *)data;
	const int m=predata->mgl_data[0].A.row;
	const int maxit=predata->maxit;
	unsigned int i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    amgparam.AMG_type = predata->AMG_type;
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
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_fmgcycle(mgl,&amgparam);
	
    fasp_array_cp(m,mgl->x.val,z);	
}

/**
 * \fn void fasp_precond_amli(double *r, double *z, void *data)
 * \brief get z from r by AMLI
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
void fasp_precond_amli (double *r, 
                        double *z, 
                        void *data)
{
	precond_data *predata=(precond_data *)data;
	const int m=predata->mgl_data[0].A.row;
	const int maxit=predata->maxit;
	unsigned int i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    amgparam.AMG_type = predata->AMG_type;
	amgparam.cycle_type = predata->cycle_type;
	amgparam.smoother   = predata->smoother;
	amgparam.presmooth_iter  = predata->presmooth_iter;
	amgparam.postsmooth_iter = predata->postsmooth_iter;
	amgparam.relaxation = predata->relaxation;
	amgparam.coarse_scaling = predata->coarse_scaling;
	amgparam.amli_degree = predata->amli_degree;
	amgparam.amli_coef = predata->amli_coef;
	amgparam.tentative_smooth = predata->tentative_smooth;
	amgparam.ILU_levels = predata->mgl_data->ILU_levels;
	
	AMG_data *mgl = predata->mgl_data;
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_amli(mgl,&amgparam,0); 
	
	fasp_array_cp(m,mgl->x.val,z);	
}

/**
 * \fn void fasp_precond_nl_amli(double *r, double *z, void *data)
 * \brief get z from r by nonliear AMLI-cycle
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 04/25/2011
 */
void fasp_precond_nl_amli (double *r, 
                           double *z, 
                           void *data)
{
	
	precond_data *predata=(precond_data *)data;
	const INT m=predata->mgl_data[0].A.row;
	const INT maxit=predata->maxit;
	const INT num_levels = predata->max_levels;
	unsigned INT i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
	amgparam.cycle_type = predata->cycle_type;
	amgparam.smoother   = predata->smoother;
	amgparam.presmooth_iter  = predata->presmooth_iter;
	amgparam.postsmooth_iter = predata->postsmooth_iter;
	amgparam.relaxation = predata->relaxation;
	amgparam.coarse_scaling = predata->coarse_scaling;
	amgparam.amli_degree = predata->amli_degree;
	amgparam.amli_coef = predata->amli_coef;
    amgparam.nl_amli_krylov_type = predata->nl_amli_krylov_type;
	amgparam.tentative_smooth = predata->tentative_smooth;
	amgparam.ILU_levels = predata->mgl_data->ILU_levels;
	
	AMG_data *mgl = predata->mgl_data;
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_nl_amli(mgl,&amgparam,0, num_levels); 
	
	fasp_array_cp(m,mgl->x.val,z);	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
