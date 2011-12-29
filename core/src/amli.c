/*! \file amli.c
 *  \brief Abstract Algebraic Multilevel Iteration (recursive version)
 *
 *  TODO: Need a non-recursive version too. --Chensong 
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_amli (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive AMLI-cycle
 *
 * \param * mgl      pointer to AMG_data data
 * \param * param    pointer to AMG parameters
 * \param   level    integer of level indicator
 *
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
void fasp_solver_amli (AMG_data *mgl, 
                       AMG_param *param, 
                       INT level)
{	
    const INT amg_type=param->AMG_type;
	const REAL relax = param->relaxation;
    const INT  print_level = param->print_level;
	const INT  smoother = param->smoother;
	const INT  degree= param->amli_degree;
	const INT  ndeg = 3;
    
    // local variables
	//INT    p_type = 1;
	REAL   alpha  = 1.0;
	REAL * coef   = param->amli_coef;
		
	dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
	dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
	
	dCSRmat *A_level0 = &mgl[level].A; // fine level matrix
	dCSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
	const INT m0 = A_level0->row, m1 = A_level1->row;
	
	ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
	REAL *r = mgl[level].w.val, *r1 = mgl[level+1].w.val+m1;
	
	if (print_level>=PRINT_MOST) printf("AMLI level %d, pre-smoother %d.\n", level, smoother);
	
    //if (param->tentative_smooth < SMALLREAL) p_type = 0;

	if (level < mgl[level].num_levels-1) { 
		
		// pre smoothing
		if (level<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
			unsigned INT steps = param->presmooth_iter;
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
				case L1_DIAG:
					fasp_smoother_dcsr_L1diag(e0, 0, m0-1, 1, A_level0, b0, steps);
					break;
				case CG:
					//fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, pow(6,l), 1e-3, NULL, 0, 1);
					fasp_solver_dcsr_pcg(A_level0, b0, e0, 1, 1e-3, NULL, 0, 1);
					break;
				default:
					printf("### ERROR: Wrong smoother type %d!\n", smoother); 
                    exit(ERROR_INPUT_PAR);
			}
		}
		
		// form residual r = b - A x
		fasp_array_cp(m0,b0->val,r); 
		fasp_blas_dcsr_aAxpy(-1.0,A_level0,e0->val,r);
		
		// restriction r1 = R*r0
		switch (amg_type)
		{		
			case UA_AMG: 
				fasp_blas_dcsr_mxv_agg(&mgl[level].R, r, b1->val);
				break;
			default:
				fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val);
				break;
		}
		
		{ 
			fasp_array_cp(m1,b1->val,r1);
			
			unsigned INT i;
			for (i=1; i<=degree; i++)
			{
				fasp_dvec_set(m1,e1,0.0);
				fasp_solver_amli(mgl, param, level+1);
				
				// b1.val = (coef[degree-i]/coef[degree])*r1 + A_level1*e1;
                // First, compute b1.val = A_level1*e1
				fasp_blas_dcsr_mxv(A_level1, e1->val, b1->val); 
                // Then, compute b1.val = b1.val + (coef[degree-i]/coef[degree])*r1
				fasp_blas_array_axpy(m1, coef[degree-i]/coef[degree], r1, b1->val); 
				
			}
			
			fasp_dvec_set(m1,e1,0.0);	
			fasp_solver_amli(mgl, param, level+1);
		}
		
		// prolongation e0 = e0 + coef(ncoef) * P * e1
		fasp_blas_array_ax(m1, coef[degree], e1->val);
		if ( param->coarse_scaling == ON ) {
			alpha = fasp_blas_array_dotprod(m1, e1->val, r1) 
                  / fasp_blas_dcsr_vmv(A_level1, e1->val, e1->val);
		}
		
		switch (amg_type)
		{
			case UA_AMG:
				fasp_blas_dcsr_aAxpy_agg(alpha, &mgl[level].P, e1->val, e0->val);
				break;
			default:
				fasp_blas_dcsr_aAxpy(alpha, &mgl[level].P, e1->val, e0->val);
				break;
		}	
		
		// post smoothing
		if (level < param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
			unsigned INT steps = param->postsmooth_iter;
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
				case L1_DIAG:
					fasp_smoother_dcsr_L1diag(e0, m0-1, 0, -1, A_level0, b0, steps);
					break;
				case CG:
					//fasp_solver_dcsr_pcg(&mgl[l].A, &mgl[l].b, &mgl[l].x, pow(6,l), 1e-3, NULL, 0, 1);
					fasp_solver_dcsr_pcg(A_level0, b0, e0, 1, 1e-3, NULL, 0, 1);
					break;
				default:
					printf("### ERROR: Wrong smoother type %d!\n", smoother); 
                    exit(ERROR_INPUT_PAR);
			}
		}
		
	}
	else // coarsest level solver
	{
#if With_DISOLVE 
        /* use Direct.lib in Windows */
		DIRECT_MUMPS(A_level0->row, A_level0->nnz, A_level0->IA, A_level0->JA, A_level0->val, 
                     b0->val, e0->val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(A_level0, b0, e0, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(A_level0, b0, e0, 0);
#else	
		/* use iterative solver on the coarest level */
		const INT csize = A_level0->row;
		const INT cmaxit = MAX(500,MIN(csize*csize, 2000)); // coarse level iteration number
		REAL ctol = param->tol; // coarse level tolerance
        
		INT flag = fasp_solver_dcsr_pcg(A_level0, b0, e0, cmaxit, ctol, NULL, 0, 1);
		
        if (flag < 0) { // If PCG does not converge, use BiCGstab as a saft net.
            flag = fasp_solver_dcsr_pvgmres (A_level0, b0, e0, cmaxit, ctol, NULL, 0, 1, 25);
        }
        
		if ( flag < 0 && print_level > PRINT_MIN ) {
			printf("### WARNING: coarse level solver does not converge in %d steps!\n", cmaxit);
		}
#endif 
	}
	
	if (print_level>=PRINT_MOST) printf("AMLI level %d, post-smoother %d.\n", level, smoother);
}

/**
 * \fn void fasp_amg_amli_coef (REAL lambda_max, REAL lambda_min, 
 *                              INT degree, REAL *coef)
 *
 * \brief Compute the coefficients of the polynomial used by AMLI-cycle
 *
 * \param *A pointer to the coefficient matrices
 * \param degree degree of the polynomial
 * 
 * \author Xiaozhe Hu
 * \date 01/23/2011
 *
 * \note: This might be only a static function? --Chensong
 */
void fasp_amg_amli_coef (REAL lambda_max, 
                         REAL lambda_min, 
                         INT degree, 
                         REAL *coef)
{
	REAL mu0 = 1.0/lambda_max, mu1 = 1.0/lambda_min;
	REAL c = (sqrt(mu0)+sqrt(mu1))*(sqrt(mu0)+sqrt(mu1));
	REAL a = (4*mu0*mu1)/(c);
	
	REAL kappa = lambda_max/lambda_min; // condition number
	REAL delta = (sqrt(kappa) - 1.0)/(sqrt(kappa)+1.0);  
	REAL b = delta*delta;
	
	if (degree == 0)
	{
		coef[0] = 0.5*(mu0+mu1);
	}
	else if (degree == 1)
	{
		coef[0] = 0.5*c;
		coef[1] = -1.0*mu0*mu1;
	}
	else if (degree > 1)
	{
		INT i;
		
		// allocate memory 
		REAL *work = (REAL *)fasp_mem_calloc(2*degree-1, sizeof(REAL));
		REAL *coef_k, *coef_km1;
		coef_k = work; coef_km1 = work+degree;
		
		// get q_k
		fasp_amg_amli_coef(lambda_max, lambda_min, degree-1, coef_k);
		// get q_km1
		fasp_amg_amli_coef(lambda_max, lambda_min, degree-2, coef_km1);
		
		// get coef		
		coef[0] = a - b*coef_km1[0] + (1+b)*coef_k[0];
		
		for (i=1; i<degree-1; i++)
		{
			coef[i] = -b*coef_km1[i] + (1+b)*coef_k[i] - a*coef_k[i-1];
		}
		
		coef[degree-1] = (1+b)*coef_k[degree-1] - a*coef_k[degree-2];
		
		coef[degree] = -a*coef_k[degree-1];
		
		// clean memory
		if (work) fasp_mem_free(work);
	}
	else 
	{
		printf("### ERROR: Wrong degree number %d for AMLI polynomial!\n", degree);
		exit(ERROR_INPUT_PAR);
	}
	
	return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
