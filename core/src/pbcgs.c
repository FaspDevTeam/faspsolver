/*! \file pbcgs.c
 *  \brief Krylov subspace methods -- Preconditioned BiCGstab.
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
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
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
 * \fn INT fasp_solver_dcsr_pbcgs (dCSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, 
 *                                 const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b 
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param u	           Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Shuo Zhang, Chensong Zhang
 * \date   09/24/2009
 *
 * Modified by Chensong Zhang on 09/09/2011
 * Modified by Chensong Zhang on 04/30/2012
 */
INT fasp_solver_dcsr_pbcgs (dCSRmat *A, 
                            dvector *b, 
                            dvector *u, 
                            precond *pc, 
                            const REAL tol,
                            const INT MaxIt, 
                            const SHORT stop_type,
                            const SHORT print_level)
{
	const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m=b->row;
	const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
	const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
	INT          iter = 0, stag, more_step, restart_step;	
    REAL         alpha, beta, relres, normr0, fenzi, fenmu, omega, reldiff;
	REAL         factor=0.0, absres,  absresp, normd, normu, temp2, infnormu;
	REAL         *uval=u->val, *bval=b->val;
    
	// allocate temp memory (need 9*m double)
	REAL *work=(REAL *)fasp_mem_calloc(9*m,sizeof(double));	
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m, *rho=t+m;
	REAL *vec=rho+m, *pp=vec+m, *sp=pp+m, *rp=sp+m;
	    
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pbcgs ...... [Start]\n");
#endif	
	
	// initialize counters
	stag=more_step=restart_step=1;
	
	// r = b-A*u
	fasp_array_cp(m,bval,r);
	fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
	
	// save the previous step
	fasp_array_cp(m,r,rp);
	
	// fenzi=(r,rho)
	fasp_array_cp(m,r,rho);	
	fenzi=fasp_blas_array_dotprod(m,r,rho);
	
	// pp=precond(p)
	fasp_array_cp(m,r,p);
	if (pc != NULL) 
        pc->fct(p,pp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,p,pp); /* No preconditioner */
	
	// t=A*pp;
	fasp_blas_dcsr_mxv(A,pp,t);
	
	// fenmu=(A*pp,rho)
	fenmu=fasp_blas_array_dotprod(m,t,rho);
    
	// (r,rho)/(A*pp,rho)
    if (ABS(fenmu)>SMALLREAL) {
        alpha=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        return ERROR_SOLVER_MISC;	        
    }
	
	// s=r-alpha t
	fasp_blas_array_axpy(m,-alpha,t,r);
	
	// sp=precond(s)
	if (pc != NULL) 
        pc->fct(r,sp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,r,sp); /* No preconditioner */
	
	// t=A*sp;
	fasp_blas_dcsr_mxv(A,sp,t);
    
	fenzi=fasp_blas_array_dotprod(m,r,t);
	fenmu=fasp_blas_array_dotprod(m,t,t);
    
    if (ABS(fenmu)>SMALLREAL) {
        omega=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        omega=1.0; 
        // We use an approximate stepsize. Take the risky. --Chensong
        // For safety, we could just quit BiCGstab:
        //     return ERROR_SOLVER_MISC;	        
    }
    
	//sol=sol+alpha pp+omega sp
	fasp_blas_array_axpy(m,alpha,pp,uval);
	fasp_blas_array_axpy(m,omega,sp,uval);	
    iter++;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(bval,t,pc->data); /* Preconditioning */
			else fasp_array_cp(m,bval,t); /* No preconditioner, B=I */
			normr0=sqrt(ABS(fasp_blas_array_dotprod(m,bval,t)));
            normr0=MAX(SMALLREAL,normr0);
			break;
		case STOP_MOD_REL_RES:
            normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
			break;
		default:
			normr0=fasp_blas_array_norm2(m,bval);
            normr0=MAX(SMALLREAL,normr0);
			break;
	}
	
	// r=s-omega t
	fasp_blas_array_axpy(m,-omega,t,r);
	
	absres=fasp_blas_array_norm2(m,r);
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(r,z,pc->data);
			else fasp_array_cp(m,r,z);
			temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
			relres=temp2/normr0; 
			break;
		case STOP_MOD_REL_RES:
			relres=absres/normu; 
			break;
		default:
			relres=absres/normr0; 
			break;
	}
	
	absresp=absres;
	
    // output iteration information if needed	
    print_itinfo(print_level,stop_type,iter,relres,absres,factor);
    
	if (relres<tol) goto FINISHED;
	
	while ( iter++ < MaxIt )
	{
		// beta = (r,rho)/(rp,rho)
		fenzi=fasp_blas_array_dotprod(m,r,rho);
		fenmu=fasp_blas_array_dotprod(m,rp,rho);
		
        if (ABS(fenmu)>SMALLREAL) {
            beta=(fenzi*alpha)/(fenmu*omega);
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// p - omega A*pp
		fasp_blas_dcsr_aAxpy(-omega,A,pp,p);
		
		// r + beta p
		fasp_blas_array_axpby(m,1.0,r,beta,p);
		
		// pp=precond(p)
        if (pc != NULL) 
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else 
            fasp_array_cp(m,p,pp); /* No preconditioner */
		
		// t = A*pp
		fasp_blas_dcsr_mxv(A,pp,t);
		
		fenmu=fasp_blas_array_dotprod(m,t,rho);	
        
        if (ABS(fenmu)>SMALLREAL) {
            alpha=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// save previous r
		fasp_array_cp(m,r,rp);
		
		// s = r - alpha t (still use r-space)
		fasp_blas_array_axpy(m,-alpha,t,r);
		
		// sp = precond(s)
		if (pc != NULL) 
            pc->fct(r,sp,pc->data); /* Apply preconditioner */
		else 
            fasp_array_cp(m,r,sp); /* No preconditioner */
	
		// t = A sp;
		fasp_blas_dcsr_mxv(A,sp,t);
		
		fenzi=fasp_blas_array_dotprod(m,r,t);
		fenmu=fasp_blas_array_dotprod(m,t,t);
        
        if (ABS(fenmu)>SMALLREAL) {
            omega=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// sol = sol + alpha pp + omega sp
		fasp_blas_array_axpby(m,alpha,pp,omega,sp);
		fasp_blas_array_axpy(m,1.0,sp,uval);
		
		// r = s - omega t
		fasp_blas_array_axpy(m,-omega,t,r);
		
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absresp;
		
		// relative difference
		normd = fasp_blas_array_norm2(m,sp);
		normu = fasp_blas_array_norm2(m,uval);
		reldiff = normd/normu;	
		
		// relative residual
		switch (stop_type) {
			case STOP_REL_PRECRES:
				if (pc == NULL) fasp_array_cp(m,r,z);
				else pc->fct(r,z,pc->data);
				temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
				relres=temp2/normr0; 
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu; 
				break;
			default:
				relres=absres/normr0; 
				break;
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, uval); 
		if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			goto FINISHED;
		}
		
		// stagnation check, if restarted too many times, return -2.
		if ( (stag<=MaxStag) && (reldiff<maxdiff) )
		{				
			if (print_level>=PRINT_MORE) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}	
			
			fasp_array_cp(m,bval,r);
			fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL)
						pc->fct(r,z,pc->data);
					else
						fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol)
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					goto FINISHED;
				}
				++stag;
				++restart_step;
			}
			
		} // end of stagnation check!
		
		// safe-guard check
		if (relres<tol) {
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,bval,r);
			fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL) 
                        pc->fct(r,z,pc->data);
					else 
                        fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				goto FINISHED;
			}
			else {
				if (print_level>PRINT_NONE) ITS_RESTART;
			}
			
			++more_step;
			++restart_step;
		} // end if safe guard
		
		absresp=absres;		
	}
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pbcgs ...... [Finish]\n");
#endif
	
	if (iter>MaxIt) 
        return ERROR_SOLVER_MAXIT;
	else 
        return iter;
}

/**
 * \fn INT fasp_solver_dbsr_pbcgs (dBSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, 
 *                                 const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b 
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param u	           Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou
 * \date   2010/10/26
 */
INT fasp_solver_dbsr_pbcgs(dBSRmat *A, 
                           dvector *b, 
                           dvector *u, 
                           const INT MaxIt, 
                           const REAL tol, 
                           precond *pc, 
                           const INT print_level, 
                           const INT stop_type)
{
	const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m=b->row;
	const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
	const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
	INT          iter = 0, stag, more_step, restart_step;	
    REAL         alpha, beta, relres, normr0, fenzi, fenmu, omega, reldiff;
	REAL         factor=0.0, absres,  absresp, normd, normu, temp2, infnormu;
	REAL         *uval=u->val, *bval=b->val;
    
	// allocate temp memory (need 9*m double)
	REAL *work=(REAL *)fasp_mem_calloc(9*m,sizeof(double));	
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m, *rho=t+m;
	REAL *vec=rho+m, *pp=vec+m, *sp=pp+m, *rp=sp+m;
    
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dbsr_pbcgs ...... [Start]\n");
#endif	
	
	// initialize counters
	stag=more_step=restart_step=1;
	
	// r = b-A*u
	fasp_array_cp(m,bval,r);
	fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
	
	// save the previous step
	fasp_array_cp(m,r,rp);
	
	// fenzi=(r,rho)
	fasp_array_cp(m,r,rho);	
	fenzi=fasp_blas_array_dotprod(m,r,rho);
	
	// pp=precond(p)
	fasp_array_cp(m,r,p);
	if (pc != NULL) 
        pc->fct(p,pp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,p,pp); /* No preconditioner */
	
	// t=A*pp;
	fasp_blas_dbsr_mxv(A,pp,t);
	
	// fenmu=(A*pp,rho)
	fenmu=fasp_blas_array_dotprod(m,t,rho);
    
	// (r,rho)/(A*pp,rho)
    if (ABS(fenmu)>SMALLREAL) {
        alpha=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        return ERROR_SOLVER_MISC;	        
    }
	
	// s=r-alpha t
	fasp_blas_array_axpy(m,-alpha,t,r);
	
	// sp=precond(s)
	if (pc != NULL) 
        pc->fct(r,sp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,r,sp); /* No preconditioner */
	
	// t=A*sp;
	fasp_blas_dbsr_mxv(A,sp,t);
    
	fenzi=fasp_blas_array_dotprod(m,r,t);
	fenmu=fasp_blas_array_dotprod(m,t,t);
    
    if (ABS(fenmu)>SMALLREAL) {
        omega=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        omega=1.0; 
        // We use an approximate stepsize. Take the risky. --Chensong
        // For safety, we could just quit BiCGstab:
        //     return ERROR_SOLVER_MISC;	        
    }
    
	//sol=sol+alpha pp+omega sp
	fasp_blas_array_axpy(m,alpha,pp,uval);
	fasp_blas_array_axpy(m,omega,sp,uval);	
    iter++;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(bval,t,pc->data); /* Preconditioning */
			else fasp_array_cp(m,bval,t); /* No preconditioner, B=I */
			normr0=sqrt(ABS(fasp_blas_array_dotprod(m,bval,t)));
            normr0=MAX(SMALLREAL,normr0);
			break;
		case STOP_MOD_REL_RES:
            normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
			break;
		default:
			normr0=fasp_blas_array_norm2(m,bval);
            normr0=MAX(SMALLREAL,normr0);
			break;
	}
	
	// r=s-omega t
	fasp_blas_array_axpy(m,-omega,t,r);
	
	absres=fasp_blas_array_norm2(m,r);
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(r,z,pc->data);
			else fasp_array_cp(m,r,z);
			temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
			relres=temp2/normr0; 
			break;
		case STOP_MOD_REL_RES:
			relres=absres/normu; 
			break;
		default:
			relres=absres/normr0; 
			break;
	}
	
	absresp=absres;
	
    // output iteration information if needed	
    print_itinfo(print_level,stop_type,iter,relres,absres,factor);
    
	if (relres<tol) goto FINISHED;
	
	while ( iter++ < MaxIt )
	{
		// beta = (r,rho)/(rp,rho)
		fenzi=fasp_blas_array_dotprod(m,r,rho);
		fenmu=fasp_blas_array_dotprod(m,rp,rho);
		
        if (ABS(fenmu)>SMALLREAL) {
            beta=(fenzi*alpha)/(fenmu*omega);
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// p - omega A*pp
		fasp_blas_dbsr_aAxpy(-omega,A,pp,p);
		
		// r + beta p
		fasp_blas_array_axpby(m,1.0,r,beta,p);
		
		// pp=precond(p)
        if (pc != NULL) 
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else 
            fasp_array_cp(m,p,pp); /* No preconditioner */
		
		// t = A*pp
		fasp_blas_dbsr_mxv(A,pp,t);
		
		fenmu=fasp_blas_array_dotprod(m,t,rho);	
        
        if (ABS(fenmu)>SMALLREAL) {
            alpha=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// save previous r
		fasp_array_cp(m,r,rp);
		
		// s = r - alpha t (still use r-space)
		fasp_blas_array_axpy(m,-alpha,t,r);
		
		// sp = precond(s)
		if (pc != NULL) 
            pc->fct(r,sp,pc->data); /* Apply preconditioner */
		else 
            fasp_array_cp(m,r,sp); /* No preconditioner */
        
		// t = A sp;
		fasp_blas_dbsr_mxv(A,sp,t);
		
		fenzi=fasp_blas_array_dotprod(m,r,t);
		fenmu=fasp_blas_array_dotprod(m,t,t);
        
        if (ABS(fenmu)>SMALLREAL) {
            omega=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// sol = sol + alpha pp + omega sp
		fasp_blas_array_axpby(m,alpha,pp,omega,sp);
		fasp_blas_array_axpy(m,1.0,sp,uval);
		
		// r = s - omega t
		fasp_blas_array_axpy(m,-omega,t,r);
		
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absresp;
		
		// relative difference
		normd = fasp_blas_array_norm2(m,sp);
		normu = fasp_blas_array_norm2(m,uval);
		reldiff = normd/normu;	
		
		// relative residual
		switch (stop_type) {
			case STOP_REL_PRECRES:
				if (pc == NULL) fasp_array_cp(m,r,z);
				else pc->fct(r,z,pc->data);
				temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
				relres=temp2/normr0; 
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu; 
				break;
			default:
				relres=absres/normr0; 
				break;
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, uval); 
		if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			goto FINISHED;
		}
		
		// stagnation check, if restarted too many times, return -2.
		if ( (stag<=MaxStag) && (reldiff<maxdiff) )
		{				
			if (print_level>=PRINT_MORE) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}	
			
			fasp_array_cp(m,bval,r);
			fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL)
						pc->fct(r,z,pc->data);
					else
						fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol)
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					goto FINISHED;
				}
				++stag;
				++restart_step;
			}
			
		} // end of stagnation check!
		
		// safe-guard check
		if (relres<tol) {
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,bval,r);
			fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL) 
                        pc->fct(r,z,pc->data);
					else 
                        fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				goto FINISHED;
			}
			else {
				if (print_level>PRINT_NONE) ITS_RESTART;
			}
			
			++more_step;
			++restart_step;
		} // end if safe guard
		
		absresp=absres;		
	}
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dbsr_pbcgs ...... [Finish]\n");
#endif
	
	if (iter>MaxIt) 
        return ERROR_SOLVER_MAXIT;
	else 
        return iter;
}

/**
 * \fn INT fasp_solver_bdcsr_pbcgs (block_dCSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                  const REAL tol, const INT MaxIt, 
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \brief A preconditioned BiCGstab method for solving Au=b 
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param u	           Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/24/2010
 *
 * Modified by Chensong Zhang on 04/30/2012
 */
INT fasp_solver_bdcsr_pbcgs (block_dCSRmat *A, 
                             dvector *b, 
                             dvector *u, 
                             precond *pc, 
                             const REAL tol,
                             const INT MaxIt, 
                             const INT stop_type,
                             const INT print_level) 
{
	const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m=b->row;
	const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
	const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
	INT          iter = 0, stag, more_step, restart_step;	
    REAL         alpha, beta, relres, normr0, fenzi, fenmu, omega, reldiff;
	REAL         factor=0.0, absres,  absresp, normd, normu, temp2, infnormu;
	REAL         *uval=u->val, *bval=b->val;
    
	// allocate temp memory (need 9*m double)
	REAL *work=(REAL *)fasp_mem_calloc(9*m,sizeof(double));	
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m, *rho=t+m;
	REAL *vec=rho+m, *pp=vec+m, *sp=pp+m, *rp=sp+m;
    
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_bdcsr_pbcgs ...... [Start]\n");
#endif	
	
	// initialize counters
	stag=more_step=restart_step=1;
	
	// r = b-A*u
	fasp_array_cp(m,bval,r);
	fasp_blas_bdcsr_aAxpy(-1.0,A,uval,r);
	
	// save the previous step
	fasp_array_cp(m,r,rp);
	
	// fenzi=(r,rho)
	fasp_array_cp(m,r,rho);	
	fenzi=fasp_blas_array_dotprod(m,r,rho);
	
	// pp=precond(p)
	fasp_array_cp(m,r,p);
	if (pc != NULL) 
        pc->fct(p,pp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,p,pp); /* No preconditioner */
	
	// t=A*pp;
	fasp_array_set(m,t,0.0);
	fasp_blas_bdcsr_aAxpy(1.0,A,pp,t);
	
	// fenmu=(A*pp,rho)
	fenmu=fasp_blas_array_dotprod(m,t,rho);
    
	// (r,rho)/(A*pp,rho)
    if (ABS(fenmu)>SMALLREAL) {
        alpha=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        return ERROR_SOLVER_MISC;	        
    }
	
	// s=r-alpha t
	fasp_blas_array_axpy(m,-alpha,t,r);
	
	// sp=precond(s)
	if (pc != NULL) 
        pc->fct(r,sp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,r,sp); /* No preconditioner */
	
	// t=A*sp;
	fasp_array_set(m,t,0.0);
	fasp_blas_bdcsr_aAxpy(1.0,A,sp,t);
    
	fenzi=fasp_blas_array_dotprod(m,r,t);
	fenmu=fasp_blas_array_dotprod(m,t,t);
    
    if (ABS(fenmu)>SMALLREAL) {
        omega=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        omega=1.0; 
        // We use an approximate stepsize. Take the risky. --Chensong
        // For safety, we could just quit BiCGstab:
        //     return ERROR_SOLVER_MISC;	        
    }
    
	//sol=sol+alpha pp+omega sp
	fasp_blas_array_axpy(m,alpha,pp,uval);
	fasp_blas_array_axpy(m,omega,sp,uval);	
    iter++;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(bval,t,pc->data); /* Preconditioning */
			else fasp_array_cp(m,bval,t); /* No preconditioner, B=I */
			normr0=sqrt(ABS(fasp_blas_array_dotprod(m,bval,t)));
            normr0=MAX(SMALLREAL,normr0);
			break;
		case STOP_MOD_REL_RES:
            normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
			break;
		default:
			normr0=fasp_blas_array_norm2(m,bval);
            normr0=MAX(SMALLREAL,normr0);
			break;
	}
	
	// r=s-omega t
	fasp_blas_array_axpy(m,-omega,t,r);
	
	absres=fasp_blas_array_norm2(m,r);
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(r,z,pc->data);
			else fasp_array_cp(m,r,z);
			temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
			relres=temp2/normr0; 
			break;
		case STOP_MOD_REL_RES:
			relres=absres/normu; 
			break;
		default:
			relres=absres/normr0; 
			break;
	}
	
	absresp=absres;
	
    // output iteration information if needed	
    print_itinfo(print_level,stop_type,iter,relres,absres,factor);
    
	if (relres<tol) goto FINISHED;
	
	while ( iter++ < MaxIt )
	{
		// beta = (r,rho)/(rp,rho)
		fenzi=fasp_blas_array_dotprod(m,r,rho);
		fenmu=fasp_blas_array_dotprod(m,rp,rho);
		
        if (ABS(fenmu)>SMALLREAL) {
            beta=(fenzi*alpha)/(fenmu*omega);
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// p - omega A*pp
		fasp_blas_bdcsr_aAxpy(-omega,A,pp,p);
		
		// r + beta p
		fasp_blas_array_axpby(m,1.0,r,beta,p);
		
		// pp=precond(p)
        if (pc != NULL) 
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else 
            fasp_array_cp(m,p,pp); /* No preconditioner */
		
		// t = A*pp
        fasp_array_set(m,t,0.0);
		fasp_blas_bdcsr_aAxpy(1.0,A,pp,t);
		
		fenmu=fasp_blas_array_dotprod(m,t,rho);	
        
        if (ABS(fenmu)>SMALLREAL) {
            alpha=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// save previous r
		fasp_array_cp(m,r,rp);
		
		// s = r - alpha t (still use r-space)
		fasp_blas_array_axpy(m,-alpha,t,r);
		
		// sp = precond(s)
		if (pc != NULL) 
            pc->fct(r,sp,pc->data); /* Apply preconditioner */
		else 
            fasp_array_cp(m,r,sp); /* No preconditioner */
        
		// t = A sp;
		fasp_array_set(m,t,0.0);
		fasp_blas_bdcsr_aAxpy(1.0,A,sp,t);
		
		fenzi=fasp_blas_array_dotprod(m,r,t);
		fenmu=fasp_blas_array_dotprod(m,t,t);
        
        if (ABS(fenmu)>SMALLREAL) {
            omega=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// sol = sol + alpha pp + omega sp
		fasp_blas_array_axpby(m,alpha,pp,omega,sp);
		fasp_blas_array_axpy(m,1.0,sp,uval);
		
		// r = s - omega t
		fasp_blas_array_axpy(m,-omega,t,r);
		
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absresp;
		
		// relative difference
		normd = fasp_blas_array_norm2(m,sp);
		normu = fasp_blas_array_norm2(m,uval);
		reldiff = normd/normu;	
		
		// relative residual
		switch (stop_type) {
			case STOP_REL_PRECRES:
				if (pc == NULL) fasp_array_cp(m,r,z);
				else pc->fct(r,z,pc->data);
				temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
				relres=temp2/normr0; 
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu; 
				break;
			default:
				relres=absres/normr0; 
				break;
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, uval); 
		if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			goto FINISHED;
		}
		
		// stagnation check, if restarted too many times, return -2.
		if ( (stag<=MaxStag) && (reldiff<maxdiff) )
		{				
			if (print_level>=PRINT_MORE) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}	
			
			fasp_array_cp(m,bval,r);
			fasp_blas_bdcsr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL)
						pc->fct(r,z,pc->data);
					else
						fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol)
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					goto FINISHED;
				}
				++stag;
				++restart_step;
			}
			
		} // end of stagnation check!
		
		// safe-guard check
		if (relres<tol) {
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,bval,r);
			fasp_blas_bdcsr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL) 
                        pc->fct(r,z,pc->data);
					else 
                        fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				goto FINISHED;
			}
			else {
				if (print_level>PRINT_NONE) ITS_RESTART;
			}
			
			++more_step;
			++restart_step;
		} // end if safe guard
		
		absresp=absres;		
	}
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_bdcsr_pbcgs ...... [Finish]\n");
#endif
	
	if (iter>MaxIt) 
        return ERROR_SOLVER_MAXIT;
	else 
        return iter;
}

/**
 * \fn INT fasp_solver_dstr_pbcgs (dSTRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, 
 *                                 const SHORT stop_type, const SHORT print_level)
 * 
 * \brief Preconditioned BiCGstab method for solving Au=b 
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param u	           Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou
 * \date   04/25/2010
 */
INT fasp_solver_dstr_pbcgs (dSTRmat *A, 
                            dvector *b, 
                            dvector *u, 
                            const INT MaxIt, 
                            const REAL tol, 
                            precond *pc, 
                            const INT print_level, 
                            const INT stop_type)
{
	const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m=b->row;
	const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
	const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
	INT          iter = 0, stag, more_step, restart_step;	
    REAL         alpha, beta, relres, normr0, fenzi, fenmu, omega, reldiff;
	REAL         factor=0.0, absres,  absresp, normd, normu, temp2, infnormu;
	REAL         *uval=u->val, *bval=b->val;
    
	// allocate temp memory (need 9*m double)
	REAL *work=(REAL *)fasp_mem_calloc(9*m,sizeof(double));	
	REAL *p=work, *z=work+m, *r=z+m, *t=r+m, *rho=t+m;
	REAL *vec=rho+m, *pp=vec+m, *sp=pp+m, *rp=sp+m;
    
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dstr_pbcgs ...... [Start]\n");
#endif	
	
	// initialize counters
	stag=more_step=restart_step=1;
	
	// r = b-A*u
	fasp_array_cp(m,bval,r);
	fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
	
	// save the previous step
	fasp_array_cp(m,r,rp);
	
	// fenzi=(r,rho)
	fasp_array_cp(m,r,rho);	
	fenzi=fasp_blas_array_dotprod(m,r,rho);
	
	// pp=precond(p)
	fasp_array_cp(m,r,p);
	if (pc != NULL) 
        pc->fct(p,pp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,p,pp); /* No preconditioner */
	
	// t=A*pp;
	fasp_array_set(m,t,0.0);
	fasp_blas_dstr_aAxpy(1.0,A,pp,t);
	
	// fenmu=(A*pp,rho)
	fenmu=fasp_blas_array_dotprod(m,t,rho);
    
	// (r,rho)/(A*pp,rho)
    if (ABS(fenmu)>SMALLREAL) {
        alpha=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        return ERROR_SOLVER_MISC;	        
    }
	
	// s=r-alpha t
	fasp_blas_array_axpy(m,-alpha,t,r);
	
	// sp=precond(s)
	if (pc != NULL) 
        pc->fct(r,sp,pc->data); /* Apply preconditioner */
	else 
        fasp_array_cp(m,r,sp); /* No preconditioner */
	
	// t=A*sp;
	fasp_array_set(m,t,0.0);
	fasp_blas_dstr_aAxpy(1.0,A,sp,t);
    
	fenzi=fasp_blas_array_dotprod(m,r,t);
	fenmu=fasp_blas_array_dotprod(m,t,t);
    
    if (ABS(fenmu)>SMALLREAL) {
        omega=fenzi/fenmu;
    }
    else {
        if (print_level>=PRINT_SOME) ITS_DIVZERO;
        omega=1.0; 
        // We use an approximate stepsize. Take the risky. --Chensong
        // For safety, we could just quit BiCGstab:
        //     return ERROR_SOLVER_MISC;	        
    }
    
	//sol=sol+alpha pp+omega sp
	fasp_blas_array_axpy(m,alpha,pp,uval);
	fasp_blas_array_axpy(m,omega,sp,uval);	
    iter++;
	
	// compute initial relative residual 
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(bval,t,pc->data); /* Preconditioning */
			else fasp_array_cp(m,bval,t); /* No preconditioner, B=I */
			normr0=sqrt(ABS(fasp_blas_array_dotprod(m,bval,t)));
            normr0=MAX(SMALLREAL,normr0);
			break;
		case STOP_MOD_REL_RES:
            normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
			break;
		default:
			normr0=fasp_blas_array_norm2(m,bval);
            normr0=MAX(SMALLREAL,normr0);
			break;
	}
	
	// r=s-omega t
	fasp_blas_array_axpy(m,-omega,t,r);
	
	absres=fasp_blas_array_norm2(m,r);
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) pc->fct(r,z,pc->data);
			else fasp_array_cp(m,r,z);
			temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
			relres=temp2/normr0; 
			break;
		case STOP_MOD_REL_RES:
			relres=absres/normu; 
			break;
		default:
			relres=absres/normr0; 
			break;
	}
	
	absresp=absres;
	
    // output iteration information if needed	
    print_itinfo(print_level,stop_type,iter,relres,absres,factor);
    
	if (relres<tol) goto FINISHED;
	
	while ( iter++ < MaxIt )
	{
		// beta = (r,rho)/(rp,rho)
		fenzi=fasp_blas_array_dotprod(m,r,rho);
		fenmu=fasp_blas_array_dotprod(m,rp,rho);
		
        if (ABS(fenmu)>SMALLREAL) {
            beta=(fenzi*alpha)/(fenmu*omega);
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// p - omega A*pp
		fasp_blas_dstr_aAxpy(-omega,A,pp,p);
		
		// r + beta p
		fasp_blas_array_axpby(m,1.0,r,beta,p);
		
		// pp=precond(p)
        if (pc != NULL) 
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else 
            fasp_array_cp(m,p,pp); /* No preconditioner */
		
		// t = A*pp
        fasp_array_set(m,t,0.0);
		fasp_blas_dstr_aAxpy(1.0,A,pp,t);
		
		fenmu=fasp_blas_array_dotprod(m,t,rho);	
        
        if (ABS(fenmu)>SMALLREAL) {
            alpha=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// save previous r
		fasp_array_cp(m,r,rp);
		
		// s = r - alpha t (still use r-space)
		fasp_blas_array_axpy(m,-alpha,t,r);
		
		// sp = precond(s)
		if (pc != NULL) 
            pc->fct(r,sp,pc->data); /* Apply preconditioner */
		else 
            fasp_array_cp(m,r,sp); /* No preconditioner */
        
		// t = A sp;
		fasp_array_set(m,t,0.0);
		fasp_blas_dstr_aAxpy(1.0,A,sp,t);
		
		fenzi=fasp_blas_array_dotprod(m,r,t);
		fenmu=fasp_blas_array_dotprod(m,t,t);
        
        if (ABS(fenmu)>SMALLREAL) {
            omega=fenzi/fenmu;
        }
        else {
            if (print_level>=PRINT_SOME) ITS_DIVZERO;
            return ERROR_SOLVER_MISC;	        
        }
		
		// sol = sol + alpha pp + omega sp
		fasp_blas_array_axpby(m,alpha,pp,omega,sp);
		fasp_blas_array_axpy(m,1.0,sp,uval);
		
		// r = s - omega t
		fasp_blas_array_axpy(m,-omega,t,r);
		
		absres=fasp_blas_array_norm2(m,r);
		factor=absres/absresp;
		
		// relative difference
		normd = fasp_blas_array_norm2(m,sp);
		normu = fasp_blas_array_norm2(m,uval);
		reldiff = normd/normu;	
		
		// relative residual
		switch (stop_type) {
			case STOP_REL_PRECRES:
				if (pc == NULL) fasp_array_cp(m,r,z);
				else pc->fct(r,z,pc->data);
				temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
				relres=temp2/normr0; 
				break;
			case STOP_MOD_REL_RES:
				relres=absres/normu; 
				break;
			default:
				relres=absres/normr0; 
				break;
		}
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, uval); 
		if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			goto FINISHED;
		}
		
		// stagnation check, if restarted too many times, return -2.
		if ( (stag<=MaxStag) && (reldiff<maxdiff) )
		{				
			if (print_level>=PRINT_MORE) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
			}	
			
			fasp_array_cp(m,bval,r);
			fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL)
						pc->fct(r,z,pc->data);
					else
						fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			if (relres<tol)
				break;
			else {
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					goto FINISHED;
				}
				++stag;
				++restart_step;
			}
			
		} // end of stagnation check!
		
		// safe-guard check
		if (relres<tol) {
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
			
			fasp_array_cp(m,bval,r);
			fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
			fasp_array_set(m,p,0.0);
			fasp_array_cp(m,p,pp);
			absres=fasp_blas_array_norm2(m,r);
			
			// relative residual
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc != NULL) 
                        pc->fct(r,z,pc->data);
					else 
                        fasp_array_cp(m,r,z);
					temp2=sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
					relres=temp2/normr0; 
					break;
				case STOP_MOD_REL_RES:
					relres=absres/normu; 
					break;
				default:
					relres=absres/normr0; 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres);
			
			// check convergence
			if (relres<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				goto FINISHED;
			}
			else {
				if (print_level>PRINT_NONE) ITS_RESTART;
			}
			
			++more_step;
			++restart_step;
		} // end if safe guard
		
		absresp=absres;		
	}
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dstr_pbcgs ...... [Finish]\n");
#endif
	
	if (iter>MaxIt) 
        return ERROR_SOLVER_MAXIT;
	else 
        return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
