/*! \file multigrid.c
 *  \brief Abstract recursive multigrid cycle
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_mgrecur (AMG_data *mgl, AMG_param *param, int level)
 * \brief Solve Ax=b with recursive multigrid k-cycle
 *
 * \param *mgl pointer to AMG_data data
 * \param *param pointer to AMG parameters
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date 04/06/2010
 */
void fasp_solver_mgrecur (AMG_data *mgl, AMG_param *param, int level)
{	

	const int print_level = param->print_level;
	const int smoother = param->smoother;
	const int cycle_type = param->cycle_type;
	const int ndeg = 3;
	const double relax = param->relaxation;
	
	dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
	dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
	
	dCSRmat *A_level0 = &mgl[level].A; // fine level matrix
	dCSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
	const int m0 = A_level0->row, m1 = A_level1->row;
	
	ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
	double *r = mgl[level].w.val;
	
	if (print_level>8) printf("MG level %d, pre-smoother %d.\n", level, smoother);
	
	if (level < mgl[level].num_levels-1) { 
		
		// pre smoothing
		if (level<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
			unsigned int steps = param->presmooth_iter;
			switch (smoother) {
				case GS:
					fasp_smoother_dcsr_gs(e0, 0, m0-1, 1, A_level0, b0, steps);
					break;
				case POLY:
					fasp_smoother_dcsr_poly(A_level0, b0, e0, m0, ndeg, steps); 
					break;
				case JACOBI:
					fasp_smoother_dcsr_jacobi(e0, 0, m0-1, 1, A_level0, b0, steps);
					break;
				case SGS:
					fasp_smoother_dcsr_sgs(e0, A_level0, b0, steps);
					break;
				case SOR:
					fasp_smoother_dcsr_sor(e0, 0, m0-1, 1, A_level0, b0, steps, relax);
					break;
				case SSOR:
					fasp_smoother_dcsr_sor(e0, 0, m0-1, 1, A_level0, b0, steps, relax);
					fasp_smoother_dcsr_sor(e0, m0-1, 0,-1, A_level0, b0, steps, relax);
					break;
				case GSOR:
					fasp_smoother_dcsr_gs(e0, 0, m0-1, 1, A_level0, b0, steps);
					fasp_smoother_dcsr_sor(e0, m0-1, 0, -1, A_level0, b0, steps, relax);
					break;
				case SGSOR:
					fasp_smoother_dcsr_gs(e0, 0, m0-1, 1, A_level0, b0, steps);
					fasp_smoother_dcsr_gs(e0, m0-1, 0,-1, A_level0, b0, steps);
					fasp_smoother_dcsr_sor(e0, 0, m0-1, 1, A_level0, b0, steps, relax);
					fasp_smoother_dcsr_sor(e0, m0-1, 0,-1, A_level0, b0, steps, relax);
					break;
				default:
					printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
			}
		}
		
		// form residual r = b - A x
		fasp_array_cp(m0,b0->val,r); 
		fasp_blas_dcsr_aAxpy(-1.0,A_level0,e0->val,r);
		
		// restriction r1 = R*r0
		fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val);
		
		{ // call MG recursively: type = 1 for V cycle, type = 2 for W cycle 
			unsigned int i;		
			fasp_dvec_set(m1,e1,0.0);	
			for (i=0; i<cycle_type; ++i) fasp_solver_mgrecur (mgl, param, level+1);
		}
		
		// prolongation e0 = e0 + P*e1
		fasp_blas_dcsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
		
		// post smoothing
		if (level < param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
			unsigned int steps = param->postsmooth_iter;
			switch (smoother) {
				case GS:
					fasp_smoother_dcsr_gs(e0, m0-1, 0, -1, A_level0, b0, steps); 
					break;
				case POLY:
					fasp_smoother_dcsr_poly(A_level0, b0, e0, m0, ndeg, steps); 
					break;
				case JACOBI:
					fasp_smoother_dcsr_jacobi(e0, m0-1, 0, -1, A_level0, b0, steps);
					break;					
				case SGS:
					fasp_smoother_dcsr_sgs(e0, A_level0, b0, steps);
					break;
				case SOR:
					fasp_smoother_dcsr_sor(e0, m0-1, 0, -1, A_level0, b0, steps, relax);
					break;
				case SSOR:
					fasp_smoother_dcsr_sor(e0, 0, m0-1, 1, A_level0, b0, steps, relax);
					fasp_smoother_dcsr_sor(e0, m0-1, 0,-1, A_level0, b0, steps, relax);
					break;
				case GSOR:
					fasp_smoother_dcsr_sor(e0, 0, m0-1, 1, A_level0, b0, steps, relax);
					fasp_smoother_dcsr_gs(e0, m0-1, 0, -1, A_level0, b0, steps);
					break;
				case SGSOR:
					fasp_smoother_dcsr_sor(e0, 0, m0-1, 1, A_level0, b0, steps, relax);
					fasp_smoother_dcsr_sor(e0, m0-1, 0,-1, A_level0, b0, steps, relax);
					fasp_smoother_dcsr_gs(e0, 0, m0-1, 1, A_level0, b0, steps);
					fasp_smoother_dcsr_gs(e0, m0-1, 0,-1, A_level0, b0, steps);
					break;
				default:
					printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
			}
		}
		
	}
	else // coarsest level solver
	{
#if With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(A_level0, b0, e0, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(A_level0, b0, e0, 0);
#else	
		/* use default iterative solver on the coarest level */
		unsigned int cmaxit = m0*m0; // coarse level iteration number
		double ctol = param->tol; // coarse level tolerance		
		fasp_solver_dcsr_pbcgs(A_level0, b0, e0, cmaxit, ctol, NULL, 0, 1);
		//fasp_solver_dcsr_pcg(A_level0, b0, e0, cmaxit, ctol, NULL, 0, 1);
#endif 
	}
	
	if (print_level>8) printf("MG level %d, post-smoother %d.\n", level, smoother);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
