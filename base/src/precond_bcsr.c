/*! \file precond_bcsr.c
 *  \brief Preconditioners
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"
//#include "forts_ns.h"
//#include "fasp4melanie.h"
//#include "fasp4melanie_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_block_diag (double *r, double *z, void *data)
 * \brief block diagonal preconditioning reservoir-reservoir block: AMG for pressure-pressure block, Jacobi for saturation-saturation block
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Xiaozhe Hu
 */
void fasp_precond_block_diag (double *r,
                            double *z, 
                            void *data)
{
	
	precond_block_data_3 *precdata=(precond_block_data_3 *)data;
	block_dCSRmat *A = precdata->Abcsr;
	dvector *tempr = &(precdata->r);
	
	const int Nx = A->blocks[0]->row;
	const int Ny = A->blocks[4]->row;
    const int Nz = A->blocks[8]->row;
	const int N = Nx + Ny + Nz;
	
	int i;
	
	//! back up r, setup z;
	fasp_array_cp(N, r, tempr->val);
	fasp_array_set(N, z, 0.0);
	
	//! prepare
	dvector rx, ry, rz, zx, zy, zz;
	rx.row = Nx; zx.row = Nx; ry.row = Ny; zy.row = Ny; rz.row = Nz; zz.row = Nz;
	rx.val = r; ry.val = &(r[Nx]); rz.val = &(r[Nx+Ny]); zx.val = z; zy.val = &(z[Nx]); zz.val = &(z[Nx+Ny]);
	
	AMG_param *amgparam = precdata->amgparam;
	AMG_data *mglx = precdata->mgl1;
    AMG_data *mgly = precdata->mgl2;
	AMG_data *mglz = precdata->mgl3;
	
	//! Preconditioning Axx block
	mglx->b.row=Nx; fasp_array_cp(Nx, rx.val, mglx->b.val); // residual is an input 
	mglx->x.row=Nx; fasp_dvec_set(Nx, &mglx->x,0.0);
	
	for(i=0;i<1;++i) fasp_solver_mgcycle(mglx, amgparam);
	fasp_array_cp(Nx, mglx->x.val, zx.val);
    
    //! Preconditioning Ayy block
	mgly->b.row=Ny; fasp_array_cp(Ny, ry.val, mgly->b.val); // residual is an input 
	mgly->x.row=Ny; fasp_dvec_set(Ny, &mgly->x,0.0);
	
	for(i=0;i<1;++i) fasp_solver_mgcycle(mgly, amgparam);
	fasp_array_cp(Ny, mgly->x.val, zy.val);
    
    //! Preconditioning Azz block
	mglz->b.row=Nz; fasp_array_cp(Nz, rz.val, mglz->b.val); // residual is an input 
	mglz->x.row=Nz; fasp_dvec_set(Nz, &mglz->x,0.0);
	
	for(i=0;i<1;++i) fasp_solver_mgcycle(mglz, amgparam);
	fasp_array_cp(Nz, mglz->x.val, zz.val);
	
	//! restore r
	fasp_array_cp(N, tempr->val, r);
	
}

/**
 * \fn void fasp_precond_block_lower (double *r, double *z, void *data)
 * \brief block diagonal preconditioning reservoir-reservoir block: AMG for pressure-pressure block, Jacobi for saturation-saturation block
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 * \author Xiaozhe Hu
 */
void fasp_precond_block_lower (double *r,
                            double *z, 
                            void *data)
{
	
	precond_block_data_3 *precdata=(precond_block_data_3 *)data;
	block_dCSRmat *A = precdata->Abcsr;
	dvector *tempr = &(precdata->r);
	
	const int Nx = A->blocks[0]->row;
	const int Ny = A->blocks[4]->row;
    const int Nz = A->blocks[8]->row;
	const int N = Nx + Ny + Nz;
	
	int i;
	
	//! back up r, setup z;
	fasp_array_cp(N, r, tempr->val);
	fasp_array_set(N, z, 0.0);
	
	//! prepare
	dvector rx, ry, rz, zx, zy, zz;
	rx.row = Nx; zx.row = Nx; ry.row = Ny; zy.row = Ny; rz.row = Nz; zz.row = Nz;
	rx.val = r; ry.val = &(r[Nx]); rz.val = &(r[Nx+Ny]); zx.val = z; zy.val = &(z[Nx]); zz.val = &(z[Nx+Ny]);
	
	AMG_param *amgparam = precdata->amgparam;
	AMG_data *mglx = precdata->mgl1;
    AMG_data *mgly = precdata->mgl2;
	AMG_data *mglz = precdata->mgl3;
	
	//! Preconditioning Axx block
	mglx->b.row=Nx; fasp_array_cp(Nx, rx.val, mglx->b.val); // residual is an input 
	mglx->x.row=Nx; fasp_dvec_set(Nx, &mglx->x,0.0);
	
	for(i=0;i<1;++i) fasp_solver_mgcycle(mglx, amgparam);
	fasp_array_cp(Nx, mglx->x.val, zx.val);
    
    // ry = ry - Ayx*zx
	fasp_blas_dcsr_aAxpy(-1.0, A->blocks[3], zx.val, ry.val);
    
    //! Preconditioning Ayy block
	mgly->b.row=Ny; fasp_array_cp(Ny, ry.val, mgly->b.val); // residual is an input 
	mgly->x.row=Ny; fasp_dvec_set(Ny, &mgly->x,0.0);
	
	for(i=0;i<1;++i) fasp_solver_mgcycle(mgly, amgparam);
	fasp_array_cp(Ny, mgly->x.val, zy.val);
    
    // rz = rz - Azx*zx - Azy*zy
	fasp_blas_dcsr_aAxpy(-1.0, A->blocks[6], zx.val, rz.val);
    fasp_blas_dcsr_aAxpy(-1.0, A->blocks[7], zy.val, rz.val);
    
    //! Preconditioning Azz block
	mglz->b.row=Nz; fasp_array_cp(Nz, rz.val, mglz->b.val); // residual is an input 
	mglz->x.row=Nz; fasp_dvec_set(Nz, &mglz->x,0.0);
	
	for(i=0;i<1;++i) fasp_solver_mgcycle(mglz, amgparam);
	fasp_array_cp(Nz, mglz->x.val, zz.val);
	
	//! restore r
	fasp_array_cp(N, tempr->val, r);
	
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
