/*! \file amlirecur.c
 *  \brief Abstract AMLI Multilevel Iteration (recursive version)
 *
 *  \note Contains AMLI and nonlinear AMLI cycles
 *
 *  TODO: Need to add a non-recursive version. --Chensong 
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
 * \fn void fasp_solver_amli (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive AMLI-cycle
 *
 * \param  mgl       pointer to AMG_data data
 * \param  param     pointer to AMG parameters
 * \param  level     integer of level indicator
 *
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
void fasp_solver_amli (AMG_data *mgl, 
                       AMG_param *param, 
                       INT level)
{	
    const SHORT  amg_type=param->AMG_type;
    const SHORT  print_level = param->print_level;
	const SHORT  smoother = param->smoother;
	const SHORT  smooth_order = param->smooth_order;
	const SHORT  degree= param->amli_degree;
	const REAL   relax = param->relaxation;
    
    // local variables
	REAL   alpha  = 1.0;
	REAL * coef   = param->amli_coef;
    
	dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
	dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
	
	dCSRmat *A_level0 = &mgl[level].A; // fine level matrix
	dCSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
	const INT m0 = A_level0->row, m1 = A_level1->row;
	
	ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
	REAL *r = mgl[level].w.val, *r1 = mgl[level+1].w.val+m1; // for residual
    INT *ordering = mgl[level].cfmark.val; // for smoother ordering
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amli ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
	if (print_level>=PRINT_MOST) printf("AMLI level %d, pre-smoother %d.\n", level, smoother);
	
	if (level < mgl[level].num_levels-1) { 
		
		// pre smoothing
		if (level<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
            fasp_dcsr_presmoothing(smoother,A_level0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,smooth_order,ordering);
		}
		
		// form residual r = b - A x
		fasp_array_cp(m0,b0->val,r); 
		fasp_blas_dcsr_aAxpy(-1.0,A_level0,e0->val,r);
		
		// restriction r1 = R*r0
		switch (amg_type)
		{		
			case UA_AMG: 
				fasp_blas_dcsr_mxv_agg(&mgl[level].R, r, b1->val); break;
			default:
				fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val); break;
		}
		
        // coarse grid correction
		{ 
			fasp_array_cp(m1,b1->val,r1);
			
			unsigned INT i;
			for (i=1; i<=degree; i++)
			{
				fasp_dvec_set(m1,e1,0.0);
				fasp_solver_amli(mgl, param, level+1);
				
				// b1 = (coef[degree-i]/coef[degree])*r1 + A_level1*e1;
                // First, compute b1 = A_level1*e1
				fasp_blas_dcsr_mxv(A_level1, e1->val, b1->val); 
                // Then, compute b1 = b1 + (coef[degree-i]/coef[degree])*r1
				fasp_blas_array_axpy(m1, coef[degree-i]/coef[degree], r1, b1->val); 
			}
			
			fasp_dvec_set(m1,e1,0.0);	
			fasp_solver_amli(mgl, param, level+1);
		}
		
		// find the optimal scaling factor alpha
		fasp_blas_array_ax(m1, coef[degree], e1->val);
		if ( param->coarse_scaling == ON ) {
			alpha = fasp_blas_array_dotprod(m1, e1->val, r1) 
                  / fasp_blas_dcsr_vmv(A_level1, e1->val, e1->val);
		}
		
		// prolongation e0 = e0 + alpha * P * e1
		switch (amg_type)
		{
			case UA_AMG:
				fasp_blas_dcsr_aAxpy_agg(alpha, &mgl[level].P, e1->val, e0->val); break;
			default:
				fasp_blas_dcsr_aAxpy(alpha, &mgl[level].P, e1->val, e0->val); break;
		}	
		
		// post smoothing
		if (level < param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
            fasp_dcsr_postsmoothing(smoother,A_level0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,smooth_order,ordering);
		}
		
	}
    
    // coarsest level solver
	else 
	{
#if With_DISOLVE 
        /* use Direct.lib in Windows */
		DIRECT_MUMPS(&A_level0->row, &A_level0->nnz, A_level0->IA, A_level0->JA, A_level0->val, 
                     b0->val, e0->val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(A_level0, b0, e0, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(A_level0, b0, e0, 0);
#else	
		/* use iterative solver on the coarest level */
        fasp_coarse_itsolver(A_level0, b0, e0, param->tol, print_level);
#endif 
	}
	
	if (print_level>=PRINT_MOST) printf("AMLI level %d, post-smoother %d.\n", level, smoother);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amli ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_solver_nl_amli (AMG_data *mgl, AMG_param *param, INT level, INT num_levels)
 * \brief Solve Ax=b with recursive nonlinear AMLI-cycle
 *
 * \param mgl         pointer to AMG_data data
 * \param param       pointer to AMG parameters
 * \param level       current level
 * \param num_levels  total numebr of levels
 *
 * \author Xiaozhe Hu
 * \date 04/06/2010
 */
void fasp_solver_nl_amli (AMG_data *mgl, 
                          AMG_param *param, 
                          INT level, 
                          INT num_levels)
{	
    const SHORT  amg_type=param->AMG_type;
	const SHORT  print_level = param->print_level;
	const SHORT  smoother = param->smoother;
	const SHORT  smooth_order = param->smooth_order;
	const REAL   relax = param->relaxation;
	
	dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
	dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
	
	dCSRmat *A_level0 = &mgl[level].A; // fine level matrix
	dCSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
	const INT m0 = A_level0->row, m1 = A_level1->row;
	
	ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
	REAL *r = mgl[level].w.val; // for residual
    INT *ordering = mgl[level].cfmark.val; // for smoother ordering
    
	dvector uH, bH;  // for coarse level correction
	uH.row = m1; uH.val = mgl[level+1].w.val + m1;
	bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
	
	if (print_level>=PRINT_MOST) 
        printf("Nonlinear AMLI level %d, pre-smoother %d.\n", level, smoother);
	
	if (level < num_levels-1) { 
		
		// pre smoothing
		if (level<param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
            fasp_dcsr_presmoothing(smoother,A_level0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,smooth_order,ordering);
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
		
		// call nonlinear AMLI-cycle recursively
		{         
			fasp_dvec_set(m1,e1,0.0);	
			
            // The coarsest problem is solved exactly.
            // No need to call krylov method on second coarest level
			if (level == num_levels-2)  
			{
				fasp_solver_nl_amli(&mgl[level+1], param, 0, num_levels-1);
			}
			else{  // recursively call preconditioned Krylov method on coarse grid
				precond_data precdata;
                
                fasp_param_amg_to_prec(&precdata, param);
				precdata.maxit = 1;
				precdata.max_levels = num_levels-1;
				precdata.mgl_data = &mgl[level+1];
                
				precond prec;
				prec.data = &precdata; 
				prec.fct = fasp_precond_nl_amli;
                
				fasp_array_cp (m1, b1->val, bH.val);
				fasp_array_cp (m1, e1->val, uH.val);
                
				const INT maxit = param->amli_degree+1;
                const REAL tol = 1e-12;
                
                switch (param->nl_amli_krylov_type)
                {
                    case SOLVER_GCG: // Use GCG
                        fasp_solver_dcsr_pgcg(A_level1,&bH,&uH,maxit,tol,&prec,0,1);
                        break;
                    default: // Use FGMRES
                        fasp_solver_dcsr_pvfgmres(A_level1,&bH,&uH,maxit,tol,&prec,0,1,30);
                        break;
                }
                
				fasp_array_cp (m1, bH.val, b1->val);
				fasp_array_cp (m1, uH.val, e1->val);
			}
			
		}
		
		// prolongation e0 = e0 + P*e1
        switch (amg_type)
		{
			case UA_AMG:
				fasp_blas_dcsr_aAxpy_agg(1.0, &mgl[level].P, e1->val, e0->val);
				break;
			default:
				fasp_blas_dcsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
				break;
		}
		
		// post smoothing
		if (level < param->ILU_levels) {
			fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
            fasp_dcsr_postsmoothing(smoother,A_level0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,smooth_order,ordering);
		}
		
	}
    
	else // coarsest level solver
	{
#if With_DISOLVE 
        /* use Direct.lib in Windows */
		DIRECT_MUMPS(&A_level0->row, &A_level0->nnz, A_level0->IA, A_level0->JA, A_level0->val, 
                     b0->val, e0->val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(A_level0, b0, e0, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(A_level0, b0, e0, 0);
#else	
		/* use iterative solver on the coarest level */
        fasp_coarse_itsolver(A_level0, b0, e0, param->tol, print_level);
#endif 
	}
	
	if (print_level>=PRINT_MOST) printf("Nonlinear AMLI level %d, post-smoother %d.\n", level, smoother);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_solver_nl_amli_bsr (AMG_data_bsr *mgl, AMG_param *param, INT level, INT num_levels)
 * \brief Solve Ax=b with recursive nonlinear AMLI-cycle
 *
 * \param mgl         pointer to AMG_data_bsr data
 * \param param       pointer to AMG parameters
 * \param level       current level
 * \param num_levels  total numebr of levels
 *
 * \author Xiaozhe Hu
 * \date 04/06/2010
 */
void fasp_solver_nl_amli_bsr (AMG_data_bsr *mgl, 
                          AMG_param *param, 
                          INT level, 
                          INT num_levels)
{	
    const SHORT  print_level = param->print_level;
	const SHORT  smoother = param->smoother;
    INT i;
	
	dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
	dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
	
	dBSRmat *A_level0 = &mgl[level].A; // fine level matrix
	dBSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
	const INT m0 = A_level0->ROW*A_level0->nb, m1 = A_level1->ROW*A_level1->nb;
	
	ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
	REAL *r = mgl[level].w.val; // for residual
    
	dvector uH, bH;  // for coarse level correction
	uH.row = m1; uH.val = mgl[level+1].w.val + m1;
	bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
	
	if (print_level>=PRINT_MORE) printf("Nonlinear AMLI level %d, pre-smoother %d.\n", level, smoother);

	if (level < num_levels-1) { 
		
		// pre smoothing
		if (level<param->ILU_levels) {
			fasp_smoother_dbsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
            unsigned int steps = param->presmooth_iter;
			
			if (steps > 0){
				switch (smoother) {
                    case GS:
                        for (i=0; i<steps; i++) fasp_smoother_dbsr_gs (A_level0, b0, e0, ASCEND, NULL);
                        break;
                    default:
						printf("### ERROR: Unsupported smoother type %d!\n", smoother); 
						exit(ERROR_INPUT_PAR);
                }
			}
		}
        
		// form residual r = b - A x
		fasp_array_cp(m0,b0->val,r); 
		fasp_blas_dbsr_aAxpy(-1.0,A_level0,e0->val,r);
		
		// restriction r1 = R*r0
        //switch (amg_type)
		//{		
		//	case UA_AMG: 
		//		fasp_blas_dcsr_mxv_agg(&mgl[level].R, r, b1->val);
		//		break;
		//	default:
				fasp_blas_dbsr_mxv(&mgl[level].R, r, b1->val);
		//		break;
		//}
		
		// call nonlinear AMLI-cycle recursively
		{ 
			fasp_dvec_set(m1,e1,0.0);	
			
            // The coarsest problem is solved exactly.
            // No need to call krylov method on second coarest level
			if (level == num_levels-2)  
			{
				fasp_solver_nl_amli_bsr(&mgl[level+1], param, 0, num_levels-1);
			}
			else{  // recursively call preconditioned Krylov method on coarse grid
				precond_data_bsr precdata;
                
                fasp_param_amg_to_prec_bsr (&precdata, param);
				precdata.maxit = 1;
				precdata.max_levels = num_levels-1;
				precdata.mgl_data = &mgl[level+1];
                
				precond prec;
				prec.data = &precdata; 
				prec.fct = fasp_precond_dbsr_nl_amli;
                
				fasp_array_cp (m1, b1->val, bH.val);
				fasp_array_cp (m1, e1->val, uH.val);
                
				const INT maxit = param->amli_degree+1;
                const REAL tol = 1e-12;
                
               switch (param->nl_amli_krylov_type)
                {
                    //case SOLVER_GCG: // Use GCG
                     //   fasp_solver_dcsr_pgcg(A_level1,&bH,&uH,maxit,tol,&prec,0,1);
                      //  break;
                    default: // Use FGMRES
                        fasp_solver_dbsr_pvfgmres(A_level1,&bH,&uH, maxit,tol,&prec,0,1, MIN(maxit,30));
                        break;
                }
                
				fasp_array_cp (m1, bH.val, b1->val);
				fasp_array_cp (m1, uH.val, e1->val);
			}
			
		}
		
		// prolongation e0 = e0 + P*e1
      //  switch (amg_type)
		//{
		//	case UA_AMG:
		//		fasp_blas_dbsr_aAxpy_agg(1.0, &mgl[level].P, e1->val, e0->val);
		//		break;
		//	default:
				fasp_blas_dbsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
		//		break;
		//}
		
		// post smoothing
		if (level < param->ILU_levels) {
			fasp_smoother_dbsr_ilu(A_level0, b0, e0, LU_level);
		}
		else {
            unsigned int steps = param->postsmooth_iter;
			
			if (steps > 0){
				switch (smoother) {
                    case GS:
                        for (i=0; i<steps; i++) fasp_smoother_dbsr_gs (A_level0, b0, e0, ASCEND, NULL);
                        break;
                    default:
                        printf("### ERROR: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
                }
			}
		}
		
	}
    
	else // coarsest level solver
	{
#if With_DISOLVE 
        /* use Direct.lib in Windows */
		DIRECT_MUMPS(&(mgl[level].Ac.row), &(mgl[level].Ac.nnz), mgl[level].Ac.IA, mgl[level].Ac.JA, mgl[level].Ac.val, 
                     b0->val, e0->val);
#elif With_UMFPACK
		/* use UMFPACK direct solver on the coarsest level */
		umfpack(&mgl[level].Ac, b0, e0, 0);
#elif With_SuperLU
		/* use SuperLU direct solver on the coarsest level */
		superlu(&mgl[level].Ac, b0, e0, 0);
#else	
		/* use iterative solver on the coarest level */
        fasp_coarse_itsolver(&mgl[level].Ac, b0, e0, param->tol, print_level);  
#endif 
	}
	
	//if (print_level>=PRINT_MOST) printf("Nonlinear AMLI level %d, post-smoother %d.\n", level, smoother);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_amg_amli_coef (const REAL lambda_max, const REAL lambda_min, 
 *                              const INT degree, REAL *coef)
 *
 * \brief Compute the coefficients of the polynomial used by AMLI-cycle
 *
 * \param lambda_max  maximal lambda
 * \param lambda_min  minimal lambda
 * \param degree      degree of polynomial approximation
 * \param coef        coefficient of AMLI (output)
 * 
 * \author Xiaozhe Hu
 * \date 01/23/2011
 */
void fasp_amg_amli_coef (const REAL lambda_max, 
                         const REAL lambda_min, 
                         const INT degree, 
                         REAL *coef)
{
	const REAL mu0 = 1.0/lambda_max, mu1 = 1.0/lambda_min;
	const REAL c = (sqrt(mu0)+sqrt(mu1))*(sqrt(mu0)+sqrt(mu1));
	const REAL a = (4*mu0*mu1)/(c);
	
	const REAL kappa = lambda_max/lambda_min; // condition number
	const REAL delta = (sqrt(kappa) - 1.0)/(sqrt(kappa)+1.0);  
	const REAL b = delta*delta;
	
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
