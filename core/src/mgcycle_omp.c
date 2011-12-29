/*! \file mgcycle_omp.c
 *  \brief Abstract non-recursive multigrid cycle
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void mgcycle_omp(AMG_data *mgl, AMG_param *param, INT nthreads, INT openmp_holds)
 * \brief Solve Ax=b with non-recursive multigrid k-cycle 
 *
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_solver_mgcycle_omp (AMG_data *mgl, 
                              AMG_param *param, 
                              INT nthreads, 
                              INT openmp_holds)
{
#if FASP_USE_OPENMP
	const INT nl = mgl[0].num_levels;
	const INT smoother = param->smoother;
	const INT smooth_order = param->smooth_order;
	const INT cycle_type = param->cycle_type;
	const INT print_level = param->print_level;
	const REAL relax = param->relaxation;
	const INT ndeg = 3;
	
	INT num_lvl[MAX_AMG_LVL] = {0}, l = 0;
	REAL alpha = 1.0;
	
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
					{
						fasp_smoother_dcsr_gs_cf_omp(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps, mgl[l].cfmark.val, 1, nthreads,openmp_holds);
					}
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
				default:
					printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
			}
		}
		
		// form residual r = b - A x
		fasp_array_cp_omp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val, nthreads,openmp_holds); 
		fasp_blas_dcsr_aAxpy_omp(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val, nthreads,openmp_holds);
		
		// restriction r1 = R*r0
		fasp_blas_dcsr_mxv_omp(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val, nthreads,openmp_holds);
		
		// prepare for the next level
		++l; fasp_dvec_set_omp(mgl[l].A.row, &mgl[l].x, 0.0, nthreads,openmp_holds);
		
	}	
	
	// CoarseSpaceSolver:	
	{
#if With_DISOLVE /* use Direct.lib in Windows */
		DIRECT_MUMPS(&mgl[nl-1].A.row, &mgl[nl-1].A.nnz, mgl[nl-1].A.IA, mgl[nl-1].A.JA, 
                     mgl[nl-1].A.val, mgl[nl-1].b.val, mgl[nl-1].x.val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#else	
		/* use default iterative solver on the coarest level */
		const INT csize = mgl[nl-1].A.row;
		const INT cmaxit = MAX(5,MIN(csize*csize, 10000)); // coarse level iteration number
		REAL ctol = param->tol; // coarse level tolerance
        
		INT flag = fasp_solver_dcsr_pbcgs (&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
		
        if ( flag < 0 ) {
            flag = fasp_solver_dcsr_pcg (&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
		}
        
        if ( flag < 0 && print_level > 0 ) {
            printf("### ERROR: coarse level solver failed!\n"); exit;
        }
#endif 
	}
	
	// BackwardSweep: 
	while (l>0) { 
		
		--l;
		
		// find the optimal scaling factor alpha
		if ( param->coarse_scaling == ON ) {
			alpha = fasp_blas_array_dotprod_omp(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val,nthreads,openmp_holds)
			/ fasp_blas_dcsr_vmv_omp(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val, nthreads,openmp_holds);
		}
		
		// prolongation u = u + alpha*P*e1
		fasp_blas_dcsr_aAxpy_omp(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val, nthreads,openmp_holds);
		
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
					{
						fasp_smoother_dcsr_gs_cf_omp(&mgl[l].x, &mgl[l].A, &mgl[l].b, steps, mgl[l].cfmark.val, -1, nthreads,openmp_holds);
					}
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
				default:
					printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
			}
		}
		
		if (num_lvl[l]<cycle_type) break;
		else num_lvl[l] = 0;
	}
	
	if (l>0) goto ForwardSweep;
#endif	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
