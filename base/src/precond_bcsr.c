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


/**
 * \fn void fasp_precond_sweeping (double *r, double *z, void *data)
 * \brief sweeping preconditioner for Maxwell equations
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date  05/01/2014
 */
void fasp_precond_sweeping (double *r,
                              double *z,
                              void *data)
{
	
	precond_sweeping_data *precdata=(precond_sweeping_data *)data;
    
    INT NumLayers = precdata->NumLayers;
	block_dCSRmat *A = precdata->A;
    block_dCSRmat *Ai = precdata->Ai;
    dCSRmat *local_A = precdata->local_A;
	ivector *local_index = precdata->local_index;
    void **local_LU = precdata->local_LU;

	dvector *r_backup = &(precdata->r);
	REAL *w = precdata->w;

    // local veriables
	INT i,l;
    dvector temp_r;
    dvector temp_e;
    
    dvector *local_r = (dvector *)fasp_mem_calloc(NumLayers, sizeof(dvector));
    dvector *local_e = (dvector *)fasp_mem_calloc(NumLayers, sizeof(dvector));

    // calculate the size and generate block local_r and local_z
    INT N=0;
    
    for (l=0;l<NumLayers; l++) {

        local_r[l].row = A->blocks[l*NumLayers+l]->row;
        local_r[l].val = r+N;
        
        local_e[l].row = A->blocks[l*NumLayers+l]->col;
        local_e[l].val = z+N;
        
        N = N+A->blocks[l*NumLayers+l]->col;
        
    }
    
    temp_r.val = w;
    temp_e.val = w+N;
	
	//! back up r, setup z;
	fasp_array_cp(N, r, r_backup->val);
    fasp_array_cp(N, r, z);
	//fasp_array_set(N, z, 0.0);
	
    // L^{-1}r
    for (l=0; l<NumLayers-1; l++){
        
        temp_r.row = local_A[l].row;
        temp_e.row = local_A[l].row;
        
        fasp_dvec_set(local_A[l].row, &temp_r, 0.0);
        
        for (i=0; i<local_e[l].row; i++){
            temp_r.val[local_index[l].val[i]] = local_e[l].val[i];
        }
        
#if  WITH_UMFPACK
        /* use UMFPACK direct solver */
        //fasp_solver_umfpack(&local_A[l], &temp_r, &temp_e, 0);
        fasp_umfpack_solve(&local_A[l], &temp_r, &temp_e, local_LU[l], 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(&local_A[l], &temp_r, &temp_e, 0);
#endif
      
        for (i=0; i<local_r[l].row; i++){
            local_r[l].val[i] = temp_e.val[local_index[l].val[i]];
        }
        
        fasp_blas_dcsr_aAxpy(-1.0, Ai->blocks[(l+1)*NumLayers+l], local_r[l].val, local_e[l+1].val);
        
    }

    // D^{-1}L^{-1}r
    for (l=0; l<NumLayers; l++){
        
        temp_r.row = local_A[l].row;
        temp_e.row = local_A[l].row;
        
        fasp_dvec_set(local_A[l].row, &temp_r, 0.0);
        
        for (i=0; i<local_e[l].row; i++){
            temp_r.val[local_index[l].val[i]] = local_e[l].val[i];
        }
        
#if  WITH_UMFPACK
        /* use UMFPACK direct solver */
        //fasp_solver_umfpack(&local_A[l], &temp_r, &temp_e, 0);
        fasp_umfpack_solve(&local_A[l], &temp_r, &temp_e, local_LU[l], 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(&local_A[l], &temp_r, &temp_e, 0);
#endif
        
        for (i=0; i<local_e[l].row; i++){
            local_e[l].val[i] = temp_e.val[local_index[l].val[i]];
        }
        
    }
    
    // L^{-t}D^{-1}L^{-1}u
    for (l=NumLayers-2; l>=0; l--){
        
        temp_r.row = local_A[l].row;
        temp_e.row = local_A[l].row;
        
        fasp_dvec_set(local_A[l].row, &temp_r, 0.0);
        
        fasp_blas_dcsr_mxv (Ai->blocks[l*NumLayers+l+1], local_e[l+1].val, local_r[l].val);
        
        for (i=0; i<local_r[l].row; i++){
            temp_r.val[local_index[l].val[i]] = local_r[l].val[i];
        }
        
#if  WITH_UMFPACK
        /* use UMFPACK direct solver */
        //fasp_solver_umfpack(&local_A[l], &temp_r, &temp_e, 0);
        fasp_umfpack_solve(&local_A[l], &temp_r, &temp_e, local_LU[l], 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(&local_A[l], &temp_r, &temp_e, 0);
#endif
       
        for (i=0; i<local_e[l].row; i++){
            local_e[l].val[i] = local_e[l].val[i] - temp_e.val[local_index[l].val[i]];
        }
        
    }
	
	//! restore r
	fasp_array_cp(N, r_backup->val, r);
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
