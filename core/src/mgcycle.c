/*! \file mgcycle.c
 *  \brief Abstract non-recursive multigrid cycle
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "mg_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_mgcycle (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle 
 *
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * \author Chensong Zhang
 * \date 10/06/2010
 *
 * \note
 *  Modified by Chensong Zhang on 12/13/2011
 */
void fasp_solver_mgcycle (AMG_data *mgl, 
                          AMG_param *param)
{	
    const SHORT  amg_type=param->AMG_type;
	const SHORT  print_level = param->print_level;
	const SHORT  smoother = param->smoother;
	const SHORT  smooth_order = param->smooth_order;
	const SHORT  cycle_type = param->cycle_type;
	const SHORT  nl = mgl[0].num_levels;
	const REAL   relax = param->relaxation;
	
    // local variables
	REAL alpha = 1.0;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_mgcycle ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    if (print_level >= PRINT_MOST) printf("AMG_level = %d, ILU_level = %d\n", nl, param->ILU_levels);
    
    INT num_lvl[MAX_AMG_LVL] = {0}, l = 0;
    
ForwardSweep:
	while (l<nl-1) { 
		num_lvl[l]++;
		
		// pre smoothing
		if (l<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
		}
		else {
            fasp_dcsr_presmoothing(smoother,&mgl[l].A,&mgl[l].b,&mgl[l].x,param->presmooth_iter,
                                   0,mgl[l].A.row-1,1,relax,smooth_order,mgl[l].cfmark.val);
		}
		
		// form residual r = b - A x
		fasp_array_cp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val); 
		fasp_blas_dcsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);
		
		// restriction r1 = R*r0
		switch (amg_type)
		{		
			case UA_AMG: 
				fasp_blas_dcsr_mxv_agg(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
				break;
			default:
				fasp_blas_dcsr_mxv(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
				break;
		}
		
		// prepare for the next level
		++l; fasp_dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);	
		
	}	
	
	// CoarseSpaceSolver:	
	{
#if With_DISOLVE 
        /* use Direct.lib in Windows */
		DIRECT_MUMPS(&mgl[nl-1].A.row, &mgl[nl-1].A.nnz, mgl[nl-1].A.IA, mgl[nl-1].A.JA, 
                     mgl[nl-1].A.val,  mgl[nl-1].b.val, mgl[nl-1].x.val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#else	
		/* use iterative solver on the coarest level */
        fasp_coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, param->tol, print_level);
#endif 
	}
	
	// BackwardSweep: 
	while (l>0) { 
		
		--l;
		
		// find the optimal scaling factor alpha
		if ( param->coarse_scaling == ON ) {
			alpha = fasp_blas_array_dotprod(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val)
			      / fasp_blas_dcsr_vmv(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val);
		}
		
		// prolongation u = u + alpha*P*e1
		switch (amg_type)
		{
			case UA_AMG:
				fasp_blas_dcsr_aAxpy_agg(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
				break;
			default:
				fasp_blas_dcsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
				break;
		}
		
		// post-smoothing
		if (l<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
		}
		else {
            fasp_dcsr_postsmoothing(smoother,&mgl[l].A,&mgl[l].b,&mgl[l].x,param->presmooth_iter,
                                    0,mgl[l].A.row-1,-1,relax,smooth_order,mgl[l].cfmark.val);
		}
		
		if (num_lvl[l]<cycle_type) break;
		else num_lvl[l] = 0;
	}
	
	if (l>0) goto ForwardSweep;
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_mgcycle ...... [Finish]\n");
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
