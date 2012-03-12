/*! \file precond.c
 *  \brief Preconditioners
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_diag (REAL *r, REAL *z, void *data)
 *
 * \brief Diagonal preconditioner z=inv(D)*r
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Chensong Zhang
 * \date 04/06/2010
 */
void fasp_precond_diag (REAL *r, 
                        REAL *z, 
                        void *data)
{
	dvector *diag=(dvector *)data;
	REAL *diagptr=diag->val;
	unsigned INT i, m=diag->row;	
	
	memcpy(z,r,m*sizeof(REAL));
	for (i=0;i<m;++i) {
		if (ABS(diag->val[i])>SMALLREAL) z[i]/=diagptr[i];
	}	
}

/**
 * \fn void fasp_precond_ilu (REAL *r, REAL *z, void *data)
 *
 * \brief preconditioning using ILU decomposition
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Shiquan Zhang
 * \date 04/06/2010
 */
void fasp_precond_ilu (REAL *r, 
                       REAL *z, 
                       void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const INT m=iludata->row, mm1=m-1, memneed=2*m;
	REAL *zz, *zr;
	
	if (iludata->nwork<memneed) goto MEMERR; // check this outside this subroutine!!
	
	zz = iludata->work; 
	zr = iludata->work+m;
	fasp_array_cp(m, r, zr); 	
	
	{
		INT i, j, jj, begin_row, end_row, mm2=m-2;
		INT *ijlu=iludata->ijlu;
		REAL *lu=iludata->luval;
		
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
	printf("### ERROR: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_ilu_forward (REAL *r, REAL *z, void *data)
 *
 * \brief preconditioning using ILU decomposition: only forwear sweep
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Xiaozhe Hu, Shiquang Zhang
 * \date 04/06/2010
 */
void fasp_precond_ilu_forward (REAL *r, 
                               REAL *z, 
                               void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const INT m=iludata->row, mm1=m-1, memneed=2*m;
	REAL *zz, *zr;
	
	if (iludata->nwork<memneed) goto MEMERR; 
	
	zz = iludata->work; 
	zr = iludata->work+m;
	fasp_array_cp(m, r, zr); 	
	
	{
		INT i, j, jj, begin_row, end_row;
		INT *ijlu=iludata->ijlu;
		REAL *lu=iludata->luval;
		
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
	printf("### ERROR: Need %d memory, only %d available!!!", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_ilu_backward (REAL *r, REAL *z, void *data)
 *
 * \brief preconditioning using ILU decomposition: only backward sweep
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Xiaozhe Hu, Shiquan  Zhang
 * \date 04/06/2010
 */
void fasp_precond_ilu_backward (REAL *r, 
                                REAL *z, 
                                void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const INT m=iludata->row, mm1=m-1, memneed=2*m;
	REAL *zz;//, *zr;
	
	if (iludata->nwork<memneed) goto MEMERR; 
	
	zz = iludata->work; 
	//zr = iludata->work+m;
	fasp_array_cp(m, r, zz); 	
	
	{
		INT i, j, jj, begin_row, end_row, mm2=m-2;
		INT *ijlu=iludata->ijlu;
		REAL *lu=iludata->luval;
		
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
	printf("### ERROR: Need %d memory, only %d available!!!", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_amg (REAL *r, REAL *z, void *data)
 *
 * \brief get z from r by classic AMG
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Chensong Zhang
 * \date 04/06/2010
 */
void fasp_precond_amg (REAL *r, 
                       REAL *z, 
                       void *data)
{
	precond_data *pcdata=(precond_data *)data;
	const INT m=pcdata->mgl_data[0].A.row;
	const INT maxit=pcdata->maxit;
	unsigned INT i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
	
	AMG_data *mgl = pcdata->mgl_data;
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_mgcycle(mgl,&amgparam); //fasp_solver_mgrecur(mgl,&amgparam,0);
	
	fasp_array_cp(m,mgl->x.val,z);	
}

/**
 * \fn void fasp_precond_famg (REAL *r, REAL *z, void *data)
 *
 * \brief get z from r by Full AMG
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 02/27/2011
 */
void fasp_precond_famg (REAL *r, 
                        REAL *z, 
                        void *data)
{
	precond_data *pcdata=(precond_data *)data;
	const INT m=pcdata->mgl_data[0].A.row;
	const INT maxit=pcdata->maxit;
	unsigned INT i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
	
	AMG_data *mgl = pcdata->mgl_data;
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_fmgcycle(mgl,&amgparam);
	
    fasp_array_cp(m,mgl->x.val,z);	
}

/**
 * \fn void fasp_precond_amli(REAL *r, REAL *z, void *data)
 *
 * \brief get z from r by AMLI
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
void fasp_precond_amli (REAL *r, 
                        REAL *z, 
                        void *data)
{
	precond_data *pcdata=(precond_data *)data;
	const INT m=pcdata->mgl_data[0].A.row;
	const INT maxit=pcdata->maxit;
	unsigned INT i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
	
	AMG_data *mgl = pcdata->mgl_data;
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_amli(mgl,&amgparam,0); 
	
	fasp_array_cp(m,mgl->x.val,z);	
}

/**
 * \fn void fasp_precond_nl_amli(REAL *r, REAL *z, void *data)
 *
 * \brief get z from r by nonliear AMLI-cycle
 *
 * \param r pointer to residual
 * \param z pointer to preconditioned residual
 * \param data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 04/25/2011
 */
void fasp_precond_nl_amli (REAL *r, 
                           REAL *z, 
                           void *data)
{
	
	precond_data *pcdata=(precond_data *)data;
	const INT m=pcdata->mgl_data[0].A.row;
	const INT maxit=pcdata->maxit;
	const SHORT num_levels = pcdata->max_levels;
	unsigned INT i;
	
	AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
	
	AMG_data *mgl = pcdata->mgl_data;
	mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
	mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
	
	for (i=0;i<maxit;++i) fasp_solver_nl_amli(mgl,&amgparam,0, num_levels); 
	
	fasp_array_cp(m,mgl->x.val,z);	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
