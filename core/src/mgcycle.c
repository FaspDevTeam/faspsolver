/*! \file mgcycle.c
 *  \brief Abstract non-recursive multigrid cycle
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if With_DISOLVE
extern "C" {
    void DIRECT_MUMPS(const int *n, const int *nnz, int *ia, int *ja, 
                      double *a, double *b, double *x);
}
#endif 

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_mgcycle (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid k-cycle 
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
	const REAL relax = param->relaxation;
	const INT nl = mgl[0].num_levels;
	const INT smoother = param->smoother;
	const INT smooth_order = param->smooth_order;
	const INT cycle_type = param->cycle_type;
	const INT print_level = param->print_level;
	const INT ndeg = 0;
	
    // local variables
	REAL alpha = 1.0;
    INT p_type = 1; // TODO: Why assign p_type this way? Input? --Chensong
	if (param->tentative_smooth < SMALLREAL) p_type = 0;
    
    if (print_level >= PRINT_MOST) {
        printf("AMG_level = %d, ILU_level = %d\n", nl, param->ILU_levels);
    }
    
    INT num_lvl[MAX_AMG_LVL] = {0}, l = 0;
    
ForwardSweep:
	while (l<nl-1) { 
		num_lvl[l]++;
		
		// pre smoothing
		if (l<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
		}
		else {
			unsigned INT steps = param->presmooth_iter;
			switch (smoother) {
				case GS:
					if (smooth_order == NO_ORDER || mgl[l].cfmark.val == NULL)
                        fasp_smoother_dcsr_gs(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps);
					else if (smooth_order == CF_ORDER)
                        fasp_smoother_dcsr_gs_cf(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps, mgl[l].cfmark.val, 1);
					break;
				case POLY:
					fasp_smoother_dcsr_poly(&mgl[l].A, &mgl[l].b, &mgl[l].x, mgl[l].A.row, ndeg, steps); 
					break;
				case JACOBI:
					fasp_smoother_dcsr_jacobi(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps);
					break;
				case SGS:
					fasp_smoother_dcsr_sgs(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps);
					break;
				case SOR:
					if (smooth_order == NO_ORDER || mgl[l].cfmark.val == NULL)				
						fasp_smoother_dcsr_sor(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps, relax);
					else if (smooth_order == CF_ORDER)
						fasp_smoother_dcsr_sor_cf(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps, relax, mgl[l].cfmark.val, 1);
					break;
				case SSOR:
					fasp_smoother_dcsr_sor(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps, relax);
					fasp_smoother_dcsr_sor(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps, relax);
					break;
				case GSOR:
					fasp_smoother_dcsr_gs(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps);
					fasp_smoother_dcsr_sor(&mgl[l].x, mgl[l].A.row-1, 0, -1, &mgl[l].A, &mgl[l].b, steps, relax);
					break;
				case SGSOR:
					fasp_smoother_dcsr_gs(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps);
					fasp_smoother_dcsr_gs(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps);
					fasp_smoother_dcsr_sor(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps, relax);
					fasp_smoother_dcsr_sor(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps, relax);
					break;
				case L1_DIAG:
					fasp_smoother_dcsr_L1diag(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps);
					break;
				case CG:
					//fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, pow(6,l), 1e-3, NULL, 0, 1);
					fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, 1, 1e-3, NULL, 0, 1);
					break;
				default:
					printf("### ERROR: wrong smoother type!\n"); 
                    exit(ERROR_INPUT_PAR);
			}
		}
		
		// form residual r = b - A x
		fasp_array_cp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val); 
		fasp_blas_dcsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);
		
		// restriction r1 = R*r0
		switch (p_type)
		{		
			case 0: 
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
                     mgl[nl-1].A.val, mgl[nl-1].b.val, mgl[nl-1].x.val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#else	
		/* use iterative solver on the coarest level */
		const INT csize = mgl[nl-1].A.row;
		const INT cmaxit = MAX(500,MIN(csize*csize, 10000)); // coarse level iteration number
		REAL ctol = param->tol; // coarse level tolerance
        
		INT flag = fasp_solver_dcsr_pcg (&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
		
        if (flag < 0) { // If PCG does not converge, use BiCGstab as a saft net.
            flag = fasp_solver_dcsr_pvgmres (&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1, 25);
        }
        
		if ( flag < 0 && print_level > PRINT_MIN ) {
			printf("### WARNING: coarse level solver does not converge in %d steps!\n", cmaxit);
		}
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
		switch (p_type)
		{
			case 0:
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
			unsigned INT steps = param->presmooth_iter;
			switch (smoother) {
				case GS:
					if (smooth_order == NO_ORDER || mgl[l].cfmark.val == NULL)
                        fasp_smoother_dcsr_gs(&mgl[l].x, mgl[l].A.row-1, 0, -1, &mgl[l].A, &mgl[l].b, steps);
					else if (smooth_order == CF_ORDER)
                        fasp_smoother_dcsr_gs_cf(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps, mgl[l].cfmark.val, -1);
					break;
				case POLY:
					fasp_smoother_dcsr_poly(&mgl[l].A, &mgl[l].b, &mgl[l].x, mgl[l].A.row, ndeg, steps); 
					break;
				case JACOBI:
					fasp_smoother_dcsr_jacobi(&mgl[l].x, mgl[l].A.row-1, 0, -1, &mgl[l].A, &mgl[l].b, steps);
					break;					
				case SGS:
					fasp_smoother_dcsr_sgs(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps);
					break;
				case SOR:
					if (smooth_order == NO_ORDER || mgl[l].cfmark.val == NULL)				
						fasp_smoother_dcsr_sor(&mgl[l].x, mgl[l].A.row-1, 0, -1, &mgl[l].A, &mgl[l].b, steps, relax);
					else if (smooth_order == CF_ORDER)
						fasp_smoother_dcsr_sor_cf(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps, relax, mgl[l].cfmark.val, -1);
					break;
				case SSOR:
					fasp_smoother_dcsr_sor(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps, relax);
					fasp_smoother_dcsr_sor(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps, relax);
					break;
				case GSOR:
					fasp_smoother_dcsr_sor(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps, relax);
					fasp_smoother_dcsr_gs(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps);
					break;
				case SGSOR:
					fasp_smoother_dcsr_sor(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps, relax);
					fasp_smoother_dcsr_sor(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps, relax);
					fasp_smoother_dcsr_gs(&mgl[l].x, 0, mgl[l].A.row-1, 1, &mgl[l].A, &mgl[l].b, steps);
					fasp_smoother_dcsr_gs(&mgl[l].x, mgl[l].A.row-1, 0,-1, &mgl[l].A, &mgl[l].b, steps);
					break;
				case L1_DIAG:
					fasp_smoother_dcsr_L1diag(&mgl[l].x, mgl[l].A.row-1, 0, -1, &mgl[l].A, &mgl[l].b, steps);
					break;
				case CG:
					//fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, pow(6,l), 1e-3, NULL, 0, 1);
					fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, 1, 1e-3, NULL, 0, 1);
					break;
				default:
					printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
			}
		}
		
		if (num_lvl[l]<cycle_type) break;
		else num_lvl[l] = 0;
	}
	
	if (l>0) goto ForwardSweep;
	
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
