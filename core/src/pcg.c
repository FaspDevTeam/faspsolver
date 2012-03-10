/*! \file pcg.c
 *  \brief Krylov subspace methods -- Preconditioned Conjugate Gradient.
 *
 *  Abstract algorithm of Krylov method    
 *
 *  Krylov method to solve A*x=b is to generate {x_k} to approximate x, 
 *  where x_k is the optimal solution in Krylov space 
 *
 *     V_k=span{r_0,A*r_0,A^2*r_0,...,A^{k-1}*r_0}, 
 *
 *  under some inner product. 
 *
 *  For the implementation, we generate a series of {p_k} such that V_k=span{p_1,...,p_k}. Details: 
 *
 *  Step 0. Given A, b, x_0, M  
 *  
 *  Step 1. Compute residual r_0 = b-A*x_0 and convergence check;  
 *  
 *  Step 2. Initialization z_0 = M^{-1}*r_0, p_0=z_0;  
 *  
 *  Step 3. Main loop ...
 *
 *  FOR k = 0:MaxIt  	
 *      - get step size alpha = f(r_k,z_k,p_k);  	
 *      - update solution: x_{k+1} = x_k + alpha*p_k;  	
 *      - perform stagnation check;  	
 *      - update residual: r_{k+1} = r_k - alpha*(A*p_k);    	
 *      - perform residual check;  	
 *      - obtain p_{k+1} using {p_0, p_1, ... , p_k};  	
 *      - prepare for next iteration;  	
 *      - print the result of k-th iteration; 
 *  END FOR
 * 
 *  Convergence check is: norm(r)/norm(b) < tol  
 *  
 *  Stagnation check is like following:    
 *      - IF norm(alpha*p_k)/norm(x_{k+1}) < tol_stag 
 *          -# compute r=b-A*x_{k+1}; 
 *          -# convergence check; 
 *          -# IF ( not converged & restart_number < Max_Stag_Check ) restart;
 *      - END IF  
 *  
 *  Residual check is like following:     
 *      - IF norm(r_{k+1})/norm(b) < tol             
 *          -# compute the real residual r = b-A*x_{k+1}; 
 *          -# convergence check; 
 *          -# IF ( not converged & restart_number < Max_Res_Check ) restart;
 *      - END IF 
 *
 *  Ref: Iterative methods for sparse linear systems (2nd Edition) 
 *  By Y. Saad, SIAM, 2003.
 *
 */  

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "its_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dcsr_pcg (dCSRmat *A, dvector *b, dvector *u, const INT MaxIt, 
 *                               const REAL tol, precond *pre, const SHORT print_level, 
 *                               const SHORT stop_type)
 *
 * \brief Preconditioned conjugate gradient (CG) method for solving Au=b 
 *
 * \param A	 pointer to the coefficient matrix
 * \param b	 pointer to the dvector of right hand side
 * \param u	 pointer to the dvector of DOFs
 * \param MaxIt integer, maximal number of iterations
 * \param tol REAL float, the tolerance for stopage
 * \param pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 * \param stop_type stopping criteria type
 *
 * \return the number of iterations
 * 
 * \author Chensong Zhang, Xiaozhe Hu, Shiquan Zhang
 * \date 05/06/2010
 */
INT fasp_solver_dcsr_pcg (dCSRmat *A, 
                          dvector *b, 
                          dvector *u, 
                          const INT MaxIt, 
                          const REAL tol,
                          precond *pre, 
                          const SHORT print_level, 
                          const SHORT stop_type)
{
	const SHORT  MaxStag=20, MaxRestartStep=20;
	const REAL   maxdiff = tol*1e-4; // staganation tolerance
	const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
	INT          iter=0, m=A->row, stag, more_step, restart_step;
	REAL         absres0=BIGREAL, absres, relres=BIGREAL, reldiff, factor;
	REAL         alpha, beta, temp1, temp2, tempr, normb=BIGREAL, normu, infnormu;
	INT          status = SUCCESS;
	
	// allocate temp memory (need 4*m REAL)
	REAL *work=(REAL *)fasp_mem_calloc(4*m,sizeof(REAL));	
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
	if (status<0) goto FINISHED;
	
#if CHMEM_MODE		
	total_alloc_mem += 4*m*sizeof(REAL);
#endif
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pcg ...... [Start]\n");
#endif	
	
	// initialization counters
	stag=1; more_step=1; restart_step=1;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pre != NULL) 
				pre->fct(b->val,t,pre->data); /* Preconditioning */
			else 
				fasp_array_cp(m,b->val,t); /* No preconditioner, B=I */
			normb=sqrt(ABS(fasp_blas_array_dotprod(m,b->val,t)));
			break;
		case STOP_MOD_REL_RES:
			break;
		default:
			normb=fasp_blas_array_norm2(m,b->val); // norm(b)
			break;
	}
	normb=MAX(SMALLREAL,normb);
	normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,u->val));
	
	// r = b-A*u
	fasp_array_cp(m,b->val,r);
	fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
	
	if (pre != NULL)
		pre->fct(r,z,pre->data); /* Preconditioning */
	else
		fasp_array_cp(m,r,z); /* No preconditioner, B=I */
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			temp2=sqrt(fasp_blas_array_dotprod(m,r,z));
			relres=temp2/normb; 
			break;
		case STOP_MOD_REL_RES:
			tempr=fasp_blas_array_norm2(m,r);
			relres=tempr/normu; 
			break;
		default:
			tempr=fasp_blas_array_norm2(m,r);
			relres=tempr/normb; 
			break;
	}
	
	if (iter<0 || relres<tol) goto FINISHED;
	
	fasp_array_cp(m,z,p);
	
	temp1=fasp_blas_array_dotprod(m,z,r);
	
	while ( iter++ < MaxIt )
	{		
		// t=A*p
		fasp_blas_dcsr_mxv(A,p,t);
		
		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=fasp_blas_array_dotprod(m,t,p);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		fasp_blas_array_axpy(m,alpha,p,u->val);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		fasp_blas_array_axpy(m,-alpha,t,r);
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absres0;
		
		// compute relative difference
		normu=fasp_blas_dvec_norm2(u);
		reldiff=ABS(alpha)*fasp_blas_array_norm2(m,p)/normu;	
		
		// compute relative residual 
		switch (stop_type) {
			case STOP_REL_PRECRES:
				// z = B*r
				if (pre == NULL)
					fasp_array_cp(m,r,z); /* No preconditioner, B=I */
				else
					pre->fct(r,z,pre->data); /* Preconditioning */
				temp2=fasp_blas_array_dotprod(m,z,r);
				relres=sqrt(ABS(temp2))/normb;
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu;
				break;
			default:
				relres=absres/normb;	
				break;				
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, u->val); 
		if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
			break;
		}
		
		// stagnation check
		if ((stag<=MaxStag) & (reldiff<maxdiff)) {
			
			if (print_level>=PRINT_MORE) { 
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}
			
			fasp_array_cp(m,b->val,r);
			fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						fasp_array_cp(m,r,z); /* No preconditioner, B=I */
					else
						pre->fct(r,z,pre->data); /* Preconditioning */
					temp2=fasp_blas_array_dotprod(m,z,r);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu;
					break;
				default:
					relres=absres/normb;	
					break;					
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol) 
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					break;
				}							
				fasp_array_set(m,p,0.0);
				++stag;
				++restart_step;
			}
		} // end of staggnation check!
		
		// safe-guard check
		if (relres<tol) 
		{
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,b->val,r);
			fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						pre->fct(r,z,pre->data); /* Preconditioning */
					else
						fasp_array_cp(m,r,z); /* No preconditioner, B=I */
					temp2=fasp_blas_array_dotprod(m,z,r);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					absres=fasp_blas_array_norm2(m,r);
					relres=absres/normu;
					break;
				default:
					absres=fasp_blas_array_norm2(m,r);
					relres=absres/normb;	
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
                        
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
				break;
			}
			
			// prepare for restarting the method
			fasp_array_set(m,p,0.0); 			
			++more_step;
			++restart_step;
			
		} // end of safe-guard check!
		
		// update relative residual here
		absres0 = absres;
		
		// compute z_k = B*r_k
		if (stop_type!=STOP_REL_PRECRES) {
			if (pre != NULL)
				pre->fct(r,z,pre->data); /* preconditioning */
			else
				fasp_array_cp(m,r,z);	 /* No preconditioner, B=I */
		}
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=fasp_blas_array_dotprod(m,z,r);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		fasp_blas_array_axpby(m,1.0,z,beta,p);
		
	} // end of main PCG loop.
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pcg ...... [Finish]\n");
#endif	
	
	if (iter>MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/**
 * \fn INT fasp_solver_bdcsr_pcg (block_dCSRmat *A, dvector *b, dvector *u, const INT MaxIt,
 *                                const REAL tol, precond *pre, const SHORT print_level, 
 *                                const SHORT stop_type)
 *
 * \brief Preconditioned conjugate gradient (CG) method for solving Au=b 
 *
 * \param A pointer to the coefficient matrix
 * \param b pointer to the dvector of right hand side
 * \param u pointer to the dvector of DOFs
 * \param MaxIt integer, maximal number of iterations
 * \param tol REAL float, the tolerance for stopage
 * \param pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 * \param stop_type stopping criteria type
 *
 * \return the number of iterations
 *
 * \author Xiaozhe Hu
 * \date 05/24/2010
 */
INT fasp_solver_bdcsr_pcg (block_dCSRmat *A, 
                           dvector *b, 
                           dvector *u, 
                           const INT MaxIt, 
                           const REAL tol,
                           precond *pre, 
                           const SHORT print_level, 
                           const SHORT stop_type)
{
	const SHORT  MaxStag=20, MaxRestartStep=20;
	const REAL   maxdiff = tol*1e-4; // staganation tolerance
	const REAL   sol_inf_tol = 1e-16; // infinity norm tolrance
    
    // local variables
	INT          iter=0, m=b->row, stag, more_step, restart_step;
	REAL         absres0=BIGREAL, absres, relres=BIGREAL, reldiff, factor;
	REAL         alpha, beta, temp1, temp2, tempr, normb=BIGREAL, normu, infnormu;
	
	// allocate temp memory (need 4*m REAL)
	REAL *work=(REAL *)fasp_mem_calloc(4*m,sizeof(REAL));		
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m;

#if CHMEM_MODE		
	total_alloc_mem += 4*m*sizeof(REAL);
#endif
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_bdcsr_pcg ...... [Start]\n");
#endif	
    
	// initialization counters
	stag=1; more_step=1; restart_step=1;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pre != NULL) 
				pre->fct(b->val,t,pre->data); /* Preconditioning */
			else 
				fasp_array_cp(m,b->val,t); /* No preconditioner, B=I */
			normb=sqrt(ABS(fasp_blas_array_dotprod(m,b->val,t)));
			break;
		case STOP_MOD_REL_RES:
			break;
		default: // STOP_REL_RES
			normb=fasp_blas_array_norm2(m,b->val); // norm(b)
			break;
	}
	normb=MAX(SMALLREAL,normb);
	normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,u->val));
	
	// r = b-A*u
	fasp_array_cp(m,b->val,r);
	fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
	
	if (pre != NULL)
		pre->fct(r,z,pre->data); /* Preconditioning */
	else
		fasp_array_cp(m,r,z); /* No preconditioner, B=I */
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			temp2=sqrt(fasp_blas_array_dotprod(m,r,z));
			relres=temp2/normb; 
			break;
		case STOP_MOD_REL_RES:
			tempr=fasp_blas_array_norm2(m,r);
			relres=tempr/normu; 
			break;
		default: // STOP_REL_RES
			tempr=fasp_blas_array_norm2(m,r);
			relres=tempr/normb; 
			break;
	}
	
	if ( iter<0 || relres<tol ) goto FINISHED;
		
	fasp_array_cp(m,z,p);
	
	temp1=fasp_blas_array_dotprod(m,z,r);
	
	while ( iter++ < MaxIt )
	{		
		// t=A*p
		fasp_array_set(m,t,0.0);
		fasp_blas_bdcsr_aAxpy(1.0,A,p,t);
		
		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=fasp_blas_array_dotprod(m,t,p);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		fasp_blas_array_axpy(m,alpha,p,u->val);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		fasp_blas_array_axpy(m,-alpha,t,r);
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absres0;
		
		// compute relative difference
		normu=fasp_blas_dvec_norm2(u);
		reldiff=ABS(alpha)*fasp_blas_array_norm2(m,p)/normu;	
		
		// compute relative residual 
		switch (stop_type) {
			case STOP_REL_PRECRES:
				// z = B*r
				if (pre == NULL)
					fasp_array_cp(m,r,z); /* No preconditioner, B=I */
				else
					pre->fct(r,z,pre->data); /* Preconditioning */
				temp2=fasp_blas_array_dotprod(m,z,r);
				relres=sqrt(ABS(temp2))/normb;
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu;
				break;
			default: // STOP_REL_RES
				relres=absres/normb;
                break;
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, u->val); 
		if (infnormu <= sol_inf_tol)
		{
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			break;
		}
		
		// stagnation check
		if ((stag<=MaxStag) & (reldiff<maxdiff)) {
			
			if (print_level>=PRINT_MORE) { 
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}
			
			fasp_array_cp(m,b->val,r);
			fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_RES:
					relres=absres/normb;	break;
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						fasp_array_cp(m,r,z); /* No preconditioner, B=I */
					else
						pre->fct(r,z,pre->data); /* Preconditioning */
					temp2=fasp_blas_array_dotprod(m,z,r);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu;
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol) 
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					break;
				}							
				fasp_array_set(m,p,0.0);
				++stag;
				++restart_step;
			}
		} // end of staggnation check!
		
		// safe-guard check
		if (relres<tol) 
		{
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,b->val,r);
			fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_RES:
					absres=fasp_blas_array_norm2(m,r);
					relres=absres/normb;	break;
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						pre->fct(r,z,pre->data); /* Preconditioning */
					else
						fasp_array_cp(m,r,z); /* No preconditioner, B=I */
					temp2=fasp_blas_array_dotprod(m,z,r);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					absres=fasp_blas_array_norm2(m,r);
					relres=absres/normu;
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				break;
			}
			
			// prepare for restarting the method
			fasp_array_set(m,p,0.0); 			
			++more_step;
			++restart_step;
			
		} // end of safe-guard check!
		
		// update relative residual here
		absres0 = absres;
		
		// compute z_k = B*r_k
		if (stop_type!=STOP_REL_PRECRES) {
			if (pre != NULL)
				pre->fct(r,z,pre->data); /* preconditioning */
			else
				fasp_array_cp(m,r,z);	 /* No preconditioner, B=I */
		}
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=fasp_blas_array_dotprod(m,z,r);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		fasp_blas_array_axpby(m,1.0,z,beta,p);
		
	} // end of main PCG loop.
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(p);
	fasp_mem_free(r);
	fasp_mem_free(z);
	fasp_mem_free(t);

#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_bdcsr_pcg ...... [Finish]\n");
#endif	
    
	return iter;
}

/**
 * \fn INT fasp_solver_dstr_pcg (dSTRmat *A, dvector *b, dvector *u, const INT MaxIt, 
 *                               const REAL tol, precond *pre, const SHORT print_level, 
 *                               const SHORT stop_type)
 *
 * \brief Preconditioned conjugate gradient (CG) method for solving Au=b 
 *
 * \param A pointer to the coefficient matrix
 * \param b pointer to the dvector of right hand side
 * \param u pointer to the dvector of DOFs
 * \param MaxIt integer, maximal number of iterations
 * \param tol REAL float, the tolerance for stopage
 * \param pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 * \param stop_type stopping criteria type
 *
 * \return the number of iterations
 *
 * \author Zhiyang Zhou
 * \date 04/25/2010
 */
INT fasp_solver_dstr_pcg (dSTRmat *A, 
                          dvector *b, 
                          dvector *u, 
                          const INT MaxIt, 
                          const REAL tol, 
                          precond *pre, 
                          const SHORT print_level, 
                          const SHORT stop_type)
{
	const SHORT MaxStag=20, MaxRestartStep=20;
	const REAL maxdiff = tol*1e-4; // staganation tolerance
	const REAL sol_inf_tol = SMALLREAL; // infinity norm tolerance
	const INT ngrid = A->ngrid; // number of grids
	const INT nc = A->nc; // size of each block (number of components)	
	
    // local variables
	INT iter=0, m=nc*ngrid, stag, morestep, restart_step;
	INT status = SUCCESS;
	REAL absres0=BIGREAL, absres, relres=BIGREAL, reldiff, factor;
	REAL alpha, beta, temp1, temp2, tempr, normb=BIGREAL, normu, infnormu;
	
	// allocate temp memory (need 4*m REAL)
	REAL *work=(REAL *)fasp_mem_calloc(4*m,sizeof(REAL));	
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
	
#if CHMEM_MODE		
	total_alloc_mem += 4*m*sizeof(REAL);
#endif
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dstr_pcg ...... [Start]\n");
#endif	
	
	// initialization counters
	stag=1; morestep=1; restart_step=1;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pre != NULL) 
				pre->fct(b->val,t,pre->data); /* Preconditioning */
			else 
				fasp_array_cp(m,b->val,t); /* No preconditioner, B=I */
			normb=sqrt(ABS(fasp_blas_array_dotprod(m,b->val,t)));
			break;
		case STOP_MOD_REL_RES:
			break;
		default:
			normb=fasp_blas_array_norm2(m,b->val); // norm(b)
			break;
	}
	normb=MAX(SMALLREAL,normb);
	normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,u->val));
	
	// r = b-A*u
	fasp_array_cp(m,b->val,r);
	fasp_blas_dstr_aAxpy(-1.0,A,u->val,r);
	
	if (pre != NULL)
		pre->fct(r,z,pre->data); /* Preconditioning */
	else
		fasp_array_cp(m,r,z); /* No preconditioner, B=I */
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			temp2=sqrt(fasp_blas_array_dotprod(m,r,z));
			relres=temp2/normb; 
			break;
		case STOP_MOD_REL_RES:
			tempr=fasp_blas_array_norm2(m,r);
			relres=tempr/normu; 
			break;
		default:
			tempr=fasp_blas_array_norm2(m,r);
			relres=tempr/normb; 
			break;
	}
	
	if (iter<0 || relres<tol) goto FINISHED;
	
	fasp_array_cp(m,z,p);
	
	temp1=fasp_blas_array_dotprod(m,z,r);
	
	while ( iter++ < MaxIt )
	{		
		// t=A*p
		fasp_array_set(m,t,0.0);
        fasp_blas_dstr_aAxpy(1.0,A,p,t);
		// spmxv_str(A,p,t); This shall be added. --Chensong
		
		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=fasp_blas_array_dotprod(m,t,p);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		fasp_blas_array_axpy(m,alpha,p,u->val);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		fasp_blas_array_axpy(m,-alpha,t,r);
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absres0;
		
		// compute relative difference
		normu=fasp_blas_dvec_norm2(u);
		reldiff=ABS(alpha)*fasp_blas_array_norm2(m,p)/normu;	
		
		// compute relative residual 
		switch (stop_type) {
			case STOP_REL_PRECRES:
				// z = B*r
				if (pre == NULL)
					fasp_array_cp(m,r,z); /* No preconditioner, B=I */
				else
					pre->fct(r,z,pre->data); /* Preconditioning */
				temp2=fasp_blas_array_dotprod(m,z,r);
				relres=sqrt(ABS(temp2))/normb;
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu;
				break;
			default:
				relres=absres/normb;	
				break;				
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, u->val); 
		if (infnormu <= sol_inf_tol)
		{
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			break;
		}
		
		// stagnation check
		if ((stag<=MaxStag) & (reldiff<maxdiff)) {
			
			if (print_level>=PRINT_MORE) { 
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}
			
			fasp_array_cp(m,b->val,r);
			fasp_blas_dstr_aAxpy(-1.0,A,u->val,r);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						fasp_array_cp(m,r,z); /* No preconditioner, B=I */
					else
						pre->fct(r,z,pre->data); /* Preconditioning */
					temp2=fasp_blas_array_dotprod(m,z,r);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu;
					break;
				default:
					relres=absres/normb;	
					break;					
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol) 
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					break;
				}							
				fasp_array_set(m,p,0.0);
				stag++;
				restart_step++;
			}
		} // end of staggnation check!
		
		// safe-guard check
		if (relres<tol) 
		{
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,b->val,r);
			fasp_blas_dstr_aAxpy(-1.0,A,u->val,r);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						pre->fct(r,z,pre->data); /* Preconditioning */
					else
						fasp_array_cp(m,r,z); /* No preconditioner, B=I */
					temp2=fasp_blas_array_dotprod(m,z,r);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					absres=fasp_blas_array_norm2(m,r);
					relres=absres/normu;
					break;
				default:
					absres=fasp_blas_array_norm2(m,r);
					relres=absres/normb;	
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);

			// check convergence
			if (relres<tol) break;
			
			if (morestep>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				break;
			}
			
			// prepare for restarting the method
			fasp_array_set(m,p,0.0); 			
			morestep++;
			restart_step++;
			
		} // end of safe-guard check!
		
		// update relative residual here
		absres0 = absres;
		
		// compute z_k = B*r_k
		if (stop_type!=STOP_REL_PRECRES) {
			if (pre != NULL)
				pre->fct(r,z,pre->data); /* preconditioning */
			else
				fasp_array_cp(m,r,z);	 /* No preconditioner, B=I */
		}
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=fasp_blas_array_dotprod(m,z,r);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		fasp_blas_array_axpby(m,1.0,z,beta,p);
		
	} // end of main PCG loop.
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dstr_pcg ...... [Finish]\n");
#endif	
	
	if (iter>MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
