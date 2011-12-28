/*! \file fmgcycle.c
 *  \brief Abstract non-recursive full multigrid cycle
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if With_DISOLVE
extern "C" {void DIRECT_MUMPS(const int *n, const int *nnz, int *ia, int *ja, double *a, double *b, double *x);}
#endif 
/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_fmgcycle(AMG_data *mgl, AMG_param *param)
 * \brief Solve Ax=b with non-recursive full multigrid k-cycle 
 *
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * \author Chensong Zhang
 * \date 02/27/2011
 */
void fasp_solver_fmgcycle (AMG_data *mgl, 
													 AMG_param *param)
{	
    const INT amg_type=param->AMG_type;
	const int nl = mgl[0].num_levels;
	const int smoother = param->smoother;
	const int smooth_order = param->smooth_order;
	const double relax = param->relaxation;
	const int ndeg = 0;
	
	// local variables
	//int p_type = 1, 
    int l = 0, i, lvl, num_cycle;
	double alpha = 1.0, relerr = BIGREAL;
	
	//if (param->tentative_smooth < SMALLREAL) p_type = 0;
	
	for ( l=0; l<nl-1; l++) { 
		// restriction r1 = R*r0
		switch (amg_type)
		{		
			case UA_AMG: 
				fasp_blas_dcsr_mxv_agg(&mgl[l].R, mgl[l].b.val, mgl[l+1].b.val);
				break;
			default:
				fasp_blas_dcsr_mxv(&mgl[l].R, mgl[l].b.val, mgl[l+1].b.val);
				break;
		}		
	}	
	
	fasp_dvec_set(mgl[l].A.row, &mgl[l].x, 0.0); // starting point
	
	// downward V-cycle
	for ( i=1; i<nl; i++ ) {
		
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
			const int csize = mgl[nl-1].A.row;
			unsigned int cmaxit = MIN(csize*csize, 1000); // coarse level iteration number
			double ctol = param->tol; // coarse level tolerance
			int flag = 0;		
			flag = fasp_solver_dcsr_pbcgs(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
			if (flag < 0)
			{
				printf("Warning: coarse level iterative solver does not converge !! (error message = %d)\n", flag);
			}
			//fasp_solver_dcsr_pcg(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
#endif 
		}
		
		// backslash 
		{
			--l; // go back to finer level
			
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
			
		}
		
		num_cycle = 0; relerr = BIGREAL;
		
		while ( relerr > 1e-8 && num_cycle < 4) {
			
			++num_cycle;
			
			// form residual r = b - A x
			fasp_array_cp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val); 
			fasp_blas_dcsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);
			relerr = fasp_blas_dvec_norm2(&mgl[l].w) / fasp_blas_dvec_norm2(&mgl[l].b);
			
			// Forward Sweep
			for ( lvl=0; lvl<i; lvl++ )  { 
				
				// pre smoothing
				if (l<param->ILU_levels) {
					fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
				}
				else {
					unsigned int steps = param->presmooth_iter;
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
							fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, steps, 1e-3, NULL, 0, 1);
							break;
						default:
							printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
					}
				} // end of smoother
				
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
				
				++l; 
				
				// prepare for the next level
				fasp_dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);	
				
			}	// end for lvl
			
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
				const int csize = mgl[nl-1].A.row;
				unsigned int cmaxit = csize*csize; // coarse level iteration number
				double ctol = param->tol; // coarse level tolerance		
				fasp_solver_dcsr_pbcgs(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
				//fasp_solver_dcsr_pcg(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, cmaxit, ctol, NULL, 0, 1);
#endif 
			}
			
			// BackwardSweep 
			for ( lvl=0; lvl<i; lvl++ ) {
				
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
					unsigned int steps = param->presmooth_iter;
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
							fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, steps, 1e-3, NULL, 0, 1);
							break;
						default:
							printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
					} // end of smoother
					
				} // end for lvl
				
			} // end while
			
		} //end while
		
	}	// end for
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
