/*! \file pcg_omp.c
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

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*-----------------------------------omp--------------------------------------*/

/**
 * \fn int fasp_solver_dcsr_pcg_omp (dCSRmat *A, dvector *b, dvector *u, const int MaxIt, const double tol, \
 *             precond *pre, const int print_level, const int stop_type, int nthreads, int openmp_holds)
 *	 \brief Preconditioned conjugate gradient (CG) method for solving Au=b 
 *	 \param A	 pointer to the coefficient matrix
 *	 \param b	 pointer to the dvector of right hand side
 *	 \param u	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *	 \return the number of iterations
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/06/2011
 */
int fasp_solver_dcsr_pcg_omp (dCSRmat *A,
                              dvector *b,
                              dvector *u,
                              const int MaxIt,
                              const double tol,
                              precond *pre,
                              const int print_level,
                              const int stop_type,
                              int nthreads,
                              int openmp_holds)
{
	int status = SUCCESS;
#if FASP_USE_OPENMP
	const int MaxStag=20, MaxRestartStep=20;
	const double maxdiff = tol*1e-4; // staganation tolerance
	const double sol_inf_tol = SMALLREAL; // infinity norm tolerance
	int iter=0, m=A->row, stag, more_step, restart_step;
	double absres0=BIGREAL, absres, relres=BIGREAL, reldiff, factor;
	double alpha, beta, temp1, temp2, tempr, normb=BIGREAL, normu, infnormu;
	
	// allocate temp memory (need 4*m double)
	double *work=(double *)fasp_mem_calloc(4*m,sizeof(double));	
	double *p=work, *z=work+m, *r=z+m, *t=r+m;
	if (status<0) goto FINISHED;
	
#if CHMEM_MODE		
	total_alloc_mem += 4*m*sizeof(double);
#endif
	
#if DEBUG_MODE
	printf("fasp_solver_dcsr_pcg ...... [Start]\n");
#endif	
	
	// initialization counters
	stag=1; more_step=1; restart_step=1;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pre != NULL) 
				pre->fct_omp(b->val,t,pre->data,nthreads,openmp_holds); /* Preconditioning */
			else 
				fasp_array_cp_omp(m,b->val,t,nthreads,openmp_holds); /* No preconditioner, B=I */
			normb=sqrt(ABS(fasp_blas_array_dotprod_omp(m,b->val,t,nthreads,openmp_holds)));
			break;
		case STOP_MOD_REL_RES:
			break;
		default:
			normb=fasp_blas_array_norm2_omp(m,b->val,nthreads,openmp_holds); // norm(b)
			break;
	}
	normb=MAX(SMALLREAL,normb);
	normu=MAX(SMALLREAL,fasp_blas_array_norm2_omp(m,u->val,nthreads,openmp_holds));
	
	// r = b-A*u
	fasp_array_cp_omp(m,b->val,r,nthreads,openmp_holds);
	fasp_blas_dcsr_aAxpy_omp(-1.0,A,u->val,r,nthreads,openmp_holds);
	
	if (pre != NULL)
		pre->fct_omp(r,z,pre->data,nthreads,openmp_holds); /* Preconditioning */
	else
		fasp_array_cp_omp(m,r,z,nthreads,openmp_holds); /* No preconditioner, B=I */
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			temp2=sqrt(fasp_blas_array_dotprod_omp(m,r,z,nthreads,openmp_holds));
			relres=temp2/normb; 
			break;
		case STOP_MOD_REL_RES:
			tempr=fasp_blas_array_norm2_omp(m,r,nthreads,openmp_holds);
			relres=tempr/normu; 
			break;
		default:
			tempr=fasp_blas_array_norm2_omp(m,r,nthreads,openmp_holds);
			relres=tempr/normb; 
			break;
	}
	
	if (iter<0 || relres<tol) goto FINISHED;
	
	fasp_array_cp_omp(m,z,p,nthreads,openmp_holds);
	
	temp1=fasp_blas_array_dotprod_omp(m,z,r,nthreads,openmp_holds);
	
	while ( iter++ < MaxIt )
	{		
		// t=A*p
		fasp_blas_dcsr_mxv_omp(A,p,t,nthreads,openmp_holds);
		
		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=fasp_blas_array_dotprod_omp(m,t,p,nthreads,openmp_holds);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		fasp_blas_array_axpy_omp(m,alpha,p,u->val,nthreads,openmp_holds);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		fasp_blas_array_axpy_omp(m,-alpha,t,r,nthreads,openmp_holds);
		absres=fasp_blas_array_norm2_omp(m,r,nthreads,openmp_holds);
		factor=absres/absres0;
		
		// compute relative difference
		normu=fasp_blas_dvec_norm2_omp(u,nthreads,openmp_holds);
		reldiff=ABS(alpha)*fasp_blas_array_norm2_omp(m,p,nthreads,openmp_holds)/normu;
		
		// compute relative residual 
		switch (stop_type) {
			case STOP_REL_PRECRES:
				// z = B*r
				if (pre == NULL)
					fasp_array_cp_omp(m,r,z,nthreads,openmp_holds); /* No preconditioner, B=I */
				else
					pre->fct_omp(r,z,pre->data,nthreads,openmp_holds); /* Preconditioning */
				temp2=fasp_blas_array_dotprod_omp(m,z,r,nthreads,openmp_holds);
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
			print_message(print_level, "PCG stops: infinity norm of the solution is too small!\n");
			iter = ERROR_SOLVER_SOLSTAG;
			break;
		}
		
		// stagnation check
		if ((stag<=MaxStag) & (reldiff<maxdiff)) {
			
			if (print_level>4) { 
				printf("PCG: restart %d caused by stagnation.\n", restart_step);
				printf("||u-u'||/||u|| = %e and computed rel. res. = %e.\n",reldiff,relres);
			}
			
			fasp_array_cp_omp(m,b->val,r,nthreads,openmp_holds);
			fasp_blas_dcsr_aAxpy_omp(-1.0,A,u->val,r,nthreads,openmp_holds);
			absres=fasp_blas_array_norm2_omp(m,r,nthreads,openmp_holds);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						fasp_array_cp_omp(m,r,z,nthreads,openmp_holds); /* No preconditioner, B=I */
					else
						pre->fct_omp(r,z,pre->data,nthreads,openmp_holds); /* Preconditioning */
					temp2=fasp_blas_array_dotprod_omp(m,z,r,nthreads,openmp_holds);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu;
					break;
				default:
					relres=absres/normb;	
					break;					
			}
			
			if (print_level>4) printf("The actual relative residual = %e\n",relres);
			
			if (relres<tol) 
				break;
			else {
				if (stag>=MaxStag) {
					print_message(print_level,"CG does not converge: staggnation error!\n");
					iter = ERROR_SOLVER_STAG;
					break;
				}							
				fasp_array_set_omp(m,p,0.0,nthreads,openmp_holds);
				++stag;
				++restart_step;
			}
		} // end of staggnation check!
		
		// safe-guard check
		if (relres<tol) 
		{
			if (print_level>4) printf("The computed relative residual = %e\n",relres);
			
			fasp_array_cp_omp(m,b->val,r,nthreads,openmp_holds);
			fasp_blas_dcsr_aAxpy_omp(-1.0,A,u->val,r,nthreads,openmp_holds);
			
			// relative residual 
			switch (stop_type) {
				case STOP_REL_PRECRES:
					// z = B*r
					if (pre == NULL)
						pre->fct_omp(r,z,pre->data,nthreads,openmp_holds); /* Preconditioning */
					else
						fasp_array_cp_omp(m,r,z,nthreads,openmp_holds); /* No preconditioner, B=I */
					temp2=fasp_blas_array_dotprod_omp(m,z,r,nthreads,openmp_holds);
					relres=sqrt(ABS(temp2))/normb;
					break;
				case STOP_MOD_REL_RES:
					absres=fasp_blas_array_norm2_omp(m,r,nthreads,openmp_holds);
					relres=absres/normu;
					break;
				default:
					absres=fasp_blas_array_norm2_omp(m,r,nthreads,openmp_holds);
					relres=absres/normb;	
					break;
			}
			
			if (print_level>4) printf("The actual relative residual = %e\n",relres);
			
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
				print_message(print_level,"Warning: the tolerence may be too small!\n");
				iter = ERROR_SOLVER_TOLSMALL;
				break;
			}
			
			// prepare for restarting the method
			fasp_array_set_omp(m,p,0.0,nthreads,openmp_holds);
			++more_step;
			++restart_step;
			
		} // end of safe-guard check!
		
		// update relative residual here
		absres0 = absres;
		
		// compute z_k = B*r_k
		if (stop_type!=STOP_REL_PRECRES) {
			if (pre != NULL)
				pre->fct_omp(r,z,pre->data,nthreads,openmp_holds); /* preconditioning */
			else
				fasp_array_cp_omp(m,r,z,nthreads,openmp_holds);	 /* No preconditioner, B=I */
		}
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=fasp_blas_array_dotprod_omp(m,z,r,nthreads,openmp_holds);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		fasp_blas_array_axpby_omp(m,1.0,z,beta,p,nthreads,openmp_holds);
		
	} // end of main PCG loop.
	
FINISHED:  // finish the iterative method
	if (print_level>0) {
		if (iter>MaxIt){
			printf("Maximal iteration %d reached with relative residual %e.\n", MaxIt, relres);
		}
		else if (iter >= 0)
			printf("Number of iterations = %d with relative residual %e.\n", iter, relres);
	}
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("fasp_solver_dcsr_pcg ...... [Finish]\n");
#endif	
	
	if (iter>MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
#else
	return status;
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
