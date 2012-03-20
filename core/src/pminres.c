/*! \file pminres.c
 *  \brief Krylov subspace methods -- Preconditioned Minimal Residual.
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
 * \fn INT fasp_solver_dcsr_pminres (dCSRmat *A, dvector *b, dvector *u, const INT MaxIt, 
 *                                   const REAL tol, precond *pc, const SHORT print_level, 
 *                                   const SHORT stop_type)
 *
 * \brief A preconditioned minimal residual (Minres) method for solving Au=b 
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param u	           Pointer to the dvector of DOFs
 * \param MaxIt        Maximal number of iterations
 * \param tol          Tolerance for stopping
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param print_level  How much information to print out
 * \param stop_type    Stopping criteria type
 *
 * \return             Number of iterations if converged, error message otherwise
 * 
 * \author Shiquan Zhang
 * \date   10/24/2010
 */
INT fasp_solver_dcsr_pminres (dCSRmat *A, 
                              dvector *b, 
                              dvector *u, 
                              const INT MaxIt, 
                              const REAL tol,
                              precond *pc, 
                              const SHORT print_level, 
                              const SHORT stop_type)
{
	const SHORT MaxStag=MIN(100, MaxIt), MaxRestartStep=MIN(5, MaxIt);
	const REAL maxdiff = tol*1e-4; // staganation tolerance
	const REAL sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
	INT iter=0, m=A->row, stag, more_step, restart_step;
	REAL absres0, absres, relres1, factor;
	REAL alpha, alpha0, alpha1, temp2, normb2, normu2, normuu, normp, infnormu;
	INT status = SUCCESS;
	
	// allocate temp memory (need 11*m REAL)
	REAL *work=(REAL *)fasp_mem_calloc(11*m,sizeof(REAL));		
	REAL *p0=work, *p1=work+m, *p2=p1+m, *z0=p2+m, *z1=z0+m;
	REAL *t0=z1+m, *t1=t0+m, *t=t1+m, *tp=t+m, *tz=tp+m, *r=tz+m;
		
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pminres ...... [Start]\n");
#endif	
	
	// initialization counters
	stag=1; more_step=1; restart_step=1;
	
	// normb2 = (b,b)
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) 
				pc->fct(b->val,t,pc->data); /* Preconditioning */
			else 
				fasp_array_cp(m,b->val,t); /* No preconditioner, B=I */
			normb2=ABS(fasp_blas_array_dotprod(m,b->val,t));
			break;
		case STOP_MOD_REL_RES:
			break;
        default: // STOP_REL_RES
			normb2=fasp_blas_array_dotprod(m,b->val, b->val); // norm(b)
			break;
	}
	normb2=MAX(SMALLREAL,normb2);
	normu2=MAX(SMALLREAL,fasp_blas_array_dotprod(m,u->val, u->val));
	
	// r = b-A*u
	fasp_array_cp(m,b->val,r);
	fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
	
	temp2=fasp_blas_array_dotprod(m,r,r);
	absres0=sqrt(temp2);
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc == NULL)
				fasp_array_cp(m,r,t);
			else
				pc->fct(r,t,pc->data);
			temp2=ABS(fasp_blas_array_dotprod(m,r,t));
			relres1=sqrt(temp2/normb2); 
			break;
		case STOP_MOD_REL_RES:
			relres1=sqrt(temp2/normu2); 
			break;
        default: // STOP_REL_RES
			relres1=sqrt(temp2/normb2); 
            break;
	}
	
	if (relres1<tol) goto FINISHED;
	
	// p0=0
	fasp_array_set(m,p0,0.0);
	
	// p1 = B*r
	if (pc == NULL) fasp_array_cp(m,r,p1); /* No preconditioner, B=I */
	else pc->fct(r,p1,pc->data); /* Preconditioning */
	
	// tp=A*p1
	fasp_blas_dcsr_mxv(A,p1,tp);
	
	// tz = B*tp
	if (pc == NULL) fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
	else pc->fct(tp,tz,pc->data); /* Preconditioning */
	
	// p1=p1/normp
	//normp=fasp_blas_array_dotprod(m,tz,tp);
	normp=ABS(fasp_blas_array_dotprod(m,tz,tp));
	normp=sqrt(normp);	
	fasp_array_cp(m,p1,t);
	fasp_array_set(m,p1,0.0);
	fasp_blas_array_axpy(m,1/normp,t,p1);
	
	// t0=A*p0=0
	fasp_array_set(m,t0,0.0);
	fasp_array_cp(m,t0,z0);
	fasp_array_cp(m,t0,t1);
	fasp_array_cp(m,t0,z1);
	
	// t1=tp/normp,z1=tz/normp
	fasp_blas_array_axpy(m,1.0/normp,tp,t1);
	fasp_blas_array_axpy(m,1.0/normp,tz,z1);
	
	if (iter < 0) goto FINISHED;
	
	while( iter++ < MaxIt) 
	{		
		// alpha=<r,z1>
		alpha=fasp_blas_array_dotprod(m,r,z1);
		
		// u=u+alpha*p1
		fasp_blas_array_axpy(m,alpha,p1,u->val);
		
		// r=r-alpha*Ap1
		fasp_blas_array_axpy(m,-alpha,t1,r);
		
		// compute t=A*z1 alpha1=<z1,t> 
		fasp_blas_dcsr_mxv(A,z1,t);
		alpha1=fasp_blas_array_dotprod(m,z1,t);
		
		// compute t=A*z0 alpha0=<z1,t> 
		fasp_blas_dcsr_mxv(A,z0,t);
		alpha0=fasp_blas_array_dotprod(m,z1,t);
		
		// p2=z1-alpha1*p1-alpha0*p0
		fasp_array_cp(m,z1,p2);
		fasp_blas_array_axpy(m,-alpha1,p1,p2);
		fasp_blas_array_axpy(m,-alpha0,p0,p2);
		
		// tp=A*p2
		fasp_blas_dcsr_mxv(A,p2,tp);
		
		// tz = B*tp
		if (pc == NULL) fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
		else pc->fct(tp,tz,pc->data); /* Preconditioning */
		
		// p2=p2/normp
		//normp=fasp_blas_array_dotprod(m,tz,tp);
		normp=ABS(fasp_blas_array_dotprod(m,tz,tp));
		normp=sqrt(normp);
		fasp_array_cp(m,p2,t);
		fasp_array_set(m,p2,0.0);
		fasp_blas_array_axpy(m,1/normp,t,p2);
		
		// prepare for the next iteration
		fasp_array_cp(m,p1,p0);
		fasp_array_cp(m,p2,p1);
		fasp_array_cp(m,t1,t0);
		fasp_array_cp(m,z1,z0);
		
		// t1=tp/normp,z1=tz/normp
		fasp_array_set(m,t1,0.0);
		fasp_array_cp(m,t1,z1);
		fasp_blas_array_axpy(m,1/normp,tp,t1);
		fasp_blas_array_axpy(m,1/normp,tz,z1);
		
		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		temp2=fasp_blas_array_dotprod(m,r,r);
		absres=sqrt(temp2);		
		
		//normb2=fasp_blas_array_dotprod(m,b->val,b->val);
		normu2=fasp_blas_dvec_dotprod(u,u);
		
		switch (stop_type) {
			case STOP_REL_PRECRES:
				if (pc == NULL)
					fasp_array_cp(m,r,t);
				else
					pc->fct(r,t,pc->data);
				temp2=ABS(fasp_blas_array_dotprod(m,r,t));
				relres1=sqrt(temp2/normb2); 
				break;
			case STOP_MOD_REL_RES:
				relres1=sqrt(temp2/normu2); 
				break;
			default: // STOP_REL_RES
				relres1=sqrt(temp2/normb2); 
				break;
		}
		
		factor=absres/absres0;
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres1,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, u->val); 
		if (infnormu <= sol_inf_tol)
		{
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			break;
		}
		
		normuu=fasp_blas_array_dotprod(m,p1,p1);
		normuu=ABS(alpha)*sqrt(normuu/normu2);
		
		// check convergence
		if (normuu<maxdiff)
		{
			if (stag<MaxStag) {
                if (print_level>=PRINT_MORE) { 
                    ITS_DIFFRES(normuu,relres1);
                    ITS_RESTART;
                }
			}
			
			fasp_array_cp(m,b->val,r); fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
			
			temp2=fasp_blas_array_dotprod(m,r,r);
			absres=sqrt(temp2);
			switch (stop_type) {
				case STOP_REL_RES:
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_REL_PRECRES:
					if (pc == NULL)
						fasp_array_cp(m,r,t);
					else
						pc->fct(r,t,pc->data);
					temp2=ABS(fasp_blas_array_dotprod(m,r,t));
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_MOD_REL_RES:
					relres1=sqrt(temp2/normu2); 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres1);
			
			if (relres1<tol)
				break;
			else
			{
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					break;
				}
				++stag;
				++restart_step;
				
				fasp_array_set(m,p0,0.0);
				
				// p1 = B*r
				if (pc == NULL)
					fasp_array_cp(m,r,p1); /* No preconditioner, B=I */
				else
					pc->fct(r,p1,pc->data); /* Preconditioning */
				
				// tp=A*p1
				fasp_blas_dcsr_mxv(A,p1,tp);
				
				// tz = B*tp
				if (pc == NULL) fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
				else pc->fct(tp,tz,pc->data); /* Preconditioning */
				
				// p1=p1/normp
				normp=fasp_blas_array_dotprod(m,tz,tp);
				normp=sqrt(normp);
				fasp_array_cp(m,p1,t);
				
				// t0=A*p0=0
				fasp_array_set(m,t0,0.0);
				fasp_array_cp(m,t0,z0);
				fasp_array_cp(m,t0,t1);
				fasp_array_cp(m,t0,z1);
				fasp_array_cp(m,t0,p1);
				
				fasp_blas_array_axpy(m,1/normp,t,p1);
				
				// t1=tp/normp,z1=tz/normp
				fasp_blas_array_axpy(m,1/normp,tp,t1);
				fasp_blas_array_axpy(m,1/normp,tz,z1);
			}
		}
		
		// safe guard
		if (relres1<tol) 
		{
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres1);
			
			fasp_array_cp(m,b->val,r); fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
			
			temp2=fasp_blas_array_dotprod(m,r,r);
			absres=sqrt(temp2);
			switch (stop_type) {
				case STOP_REL_RES:
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_REL_PRECRES:
					if (pc == NULL)
						fasp_array_cp(m,r,t);
					else
						pc->fct(r,t,pc->data);
					temp2=ABS(fasp_blas_array_dotprod(m,r,t));
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_MOD_REL_RES:
					relres1=sqrt(temp2/normu2); 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres1);
			
			// check convergence
			if (relres1<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				break;
			}
			
			if (more_step<MaxRestartStep) {
				if (print_level>PRINT_NONE) ITS_RESTART;
			}
			
			++more_step;
			++restart_step;
			
			fasp_array_set(m,p0,0.0);
			
			// p1 = B*r
			if (pc == NULL)
				fasp_array_cp(m,r,p1); /* No preconditioner, B=I */
			else
				pc->fct(r,p1,pc->data); /* Preconditioning */
			
			// tp = A*p1
			fasp_blas_dcsr_mxv(A,p1,tp);
			
			// tz = B*tp
			if (pc == NULL) fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
			else pc->fct(tp,tz,pc->data); /* Preconditioning */
			
			// p1 = p1/normp
			normp=fasp_blas_array_dotprod(m,tz,tp);
			normp=sqrt(normp);
			fasp_array_cp(m,p1,t);
			
			// t0=A*p0=0
			fasp_array_set(m,t0,0.0);
			fasp_array_cp(m,t0,z0);
			fasp_array_cp(m,t0,t1);
			fasp_array_cp(m,t0,z1);
			fasp_array_cp(m,t0,p1);
			
			fasp_blas_array_axpy(m,1/normp,t,p1);
			
			// t1=tp/normp,z1=tz/normp
			fasp_blas_array_axpy(m,1/normp,tp,t1);
			fasp_blas_array_axpy(m,1/normp,tz,z1);
			
		}
		// update relative residual here
		absres0 = absres;
	}
	
FINISHED:  // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres1);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pminres ...... [Finish]\n");
#endif
	
	if (iter>MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/**
 * \fn INT fasp_solver_bdcsr_pminres (block_dCSRmat *A, dvector *b, dvector *u, const INT MaxIt,
 *                                    const REAL tol, precond *pc, const SHORT print_level, 
 *                                    const SHORT stop_type)
 *
 * \brief A preconditioned minimal residual (Minres) method for solving Au=b 
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param u	           Pointer to the dvector of DOFs
 * \param MaxIt        Maximal number of iterations
 * \param tol          Tolerance for stopping
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param print_level  How much information to print out
 * \param stop_type    Stopping criteria type
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/24/2010
 */
INT fasp_solver_bdcsr_pminres (block_dCSRmat *A, 
                               dvector *b, 
                               dvector *u, 
                               const INT MaxIt, 
                               const REAL tol,
                               precond *pc, 
                               const SHORT print_level, 
                               const SHORT stop_type)
{
	const SHORT MaxStag=100, MaxRestartStep=20;
	const REAL maxdiff = tol*1e-4; // staganation tolerance
	const REAL sol_inf_tol = 1e-16; // infinity norm tolerance 
    
    // local variables
	INT iter=0, m=b->row, stag, more_step, restart_step;
	REAL absres0, absres, relres1, factor;
	REAL alpha, alpha0, alpha1, temp2, normb2, normu2, normuu, normp, infnormu;
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_bdcsr_pminres ...... [Start]\n");
#endif	

	REAL *p0=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *p1=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *p2=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *z0=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *z1=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *t0=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *t1=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *t =(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *tp=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *tz=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	REAL *r=(REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
	// initialization counters
	stag=1; more_step=1; restart_step=1;
	
	// normb2 = (b,b)
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc != NULL) 
				pc->fct(b->val,t,pc->data); /* Preconditioning */
			else 
				fasp_array_cp(m,b->val,t); /* No preconditioner, B=I */
			normb2=ABS(fasp_blas_array_dotprod(m,b->val,t));
			break;
		case STOP_MOD_REL_RES:
			break;
		default: // STOP_REL_RES
			normb2=fasp_blas_array_dotprod(m,b->val, b->val); // norm(b)
			break;
	}
	normb2=MAX(SMALLREAL,normb2);
	normu2=MAX(SMALLREAL,fasp_blas_array_dotprod(m,u->val, u->val));
	
	// r = b-A*u
	fasp_array_cp(m,b->val,r);
	fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
	
	temp2=fasp_blas_array_dotprod(m,r,r);
	absres0=sqrt(temp2);
	
	switch (stop_type) {
		case STOP_REL_PRECRES:
			if (pc == NULL)
				fasp_array_cp(m,r,t);
			else
				pc->fct(r,t,pc->data);
			temp2=ABS(fasp_blas_array_dotprod(m,r,t));
			relres1=sqrt(temp2/normb2); 
			break;
		case STOP_MOD_REL_RES:
			relres1=sqrt(temp2/normu2); 
			break;
		default: // STOP_REL_RES
			relres1=sqrt(temp2/normb2); 
			break;
	}
	
	if (relres1<tol) goto FINISHED;
	
	// p0=0
	fasp_array_set(m,p0,0.0);
	
	// p1 = B*r
	if (pc == NULL)
		fasp_array_cp(m,r,p1); /* No preconditioner, B=I */
	else
		pc->fct(r,p1,pc->data); /* Preconditioning */
	
	// tp=A*p1
	fasp_array_set(m,tp,0.0);
	fasp_blas_bdcsr_aAxpy(1.0,A,p1,tp);
	
	// tz = B*tp
	if (pc == NULL)
		fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
	else
		pc->fct(tp,tz,pc->data); /* Preconditioning */
	
	// p1=p1/normp
	normp=ABS(fasp_blas_array_dotprod(m,tz,tp));
	normp=sqrt(normp);	
	fasp_array_cp(m,p1,t);
	fasp_array_set(m,p1,0.0);
	fasp_blas_array_axpy(m,1/normp,t,p1);
	
	// t0=A*p0=0
	fasp_array_set(m,t0,0.0);
	fasp_array_cp(m,t0,z0);
	fasp_array_cp(m,t0,t1);
	fasp_array_cp(m,t0,z1);
	
	// t1=tp/normp,z1=tz/normp
	fasp_blas_array_axpy(m,1.0/normp,tp,t1);
	fasp_blas_array_axpy(m,1.0/normp,tz,z1);
	
	if (iter < 0) goto FINISHED;
	
	while( iter++ < MaxIt) 
	{		
		// alpha=<r,z1>
		alpha=fasp_blas_array_dotprod(m,r,z1);
		
		// u=u+alpha*p1
		fasp_blas_array_axpy(m,alpha,p1,u->val);
		
		// r=r-alpha*Ap1
		fasp_blas_array_axpy(m,-alpha,t1,r);
		
		// compute t=A*z1 alpha1=<z1,t> 
		fasp_array_set(m,t,0.0);
		fasp_blas_bdcsr_aAxpy(1.0,A,z1,t);
		alpha1=fasp_blas_array_dotprod(m,z1,t);
		
		// compute t=A*z0 alpha0=<z1,t> 
		fasp_array_set(m,t,0.0);
		fasp_blas_bdcsr_aAxpy(1.0,A,z0,t);
		alpha0=fasp_blas_array_dotprod(m,z1,t);
		
		// p2=z1-alpha1*p1-alpha0*p0
		fasp_array_cp(m,z1,p2);
		fasp_blas_array_axpy(m,-alpha1,p1,p2);
		fasp_blas_array_axpy(m,-alpha0,p0,p2);
		
		// tp=A*p2
		fasp_array_set(m,tp,0.0); fasp_blas_bdcsr_aAxpy(1.0,A,p2,tp);
		
		// tz = B*tp
		if (pc == NULL)
			fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
		else
			pc->fct(tp,tz,pc->data); /* Preconditioning */
		
		// p2=p2/normp
		//normp=fasp_blas_array_dotprod(m,tz,tp);
		normp=ABS(fasp_blas_array_dotprod(m,tz,tp));
		normp=sqrt(normp);
		fasp_array_cp(m,p2,t);
		fasp_array_set(m,p2,0.0);
		fasp_blas_array_axpy(m,1/normp,t,p2);
		
		// prepare for the next iteration
		fasp_array_cp(m,p1,p0);
		fasp_array_cp(m,p2,p1);
		fasp_array_cp(m,t1,t0);
		fasp_array_cp(m,z1,z0);
		
		// t1=tp/normp,z1=tz/normp
		fasp_array_set(m,t1,0.0);
		fasp_array_cp(m,t1,z1);
		fasp_blas_array_axpy(m,1/normp,tp,t1);
		fasp_blas_array_axpy(m,1/normp,tz,z1);
		
		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		temp2=fasp_blas_array_dotprod(m,r,r);
		absres=sqrt(temp2);		
		//normb2=fasp_blas_array_dotprod(m,b->val,b->val);
		normu2=fasp_blas_dvec_dotprod(u,u);
		
		switch (stop_type) {
			case STOP_REL_RES:
				relres1=sqrt(temp2/normb2); 
				break;
			case STOP_REL_PRECRES:
				if (pc == NULL)
					fasp_array_cp(m,r,t);
				else
					pc->fct(r,t,pc->data);
				temp2=ABS(fasp_blas_array_dotprod(m,r,t));
				relres1=sqrt(temp2/normb2); 
				break;
			case STOP_MOD_REL_RES:
				relres1=sqrt(temp2/normu2); 
				break;
		}
		
		factor=absres/absres0;
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter,relres1,absres,factor);
		
		// solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
		infnormu = fasp_blas_array_norminf(m, u->val); 
		if (infnormu <= sol_inf_tol)
		{
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
			iter = ERROR_SOLVER_SOLSTAG;
			break;
		}
		
		normuu=fasp_blas_array_dotprod(m,p1,p1);
		normuu=ABS(alpha)*sqrt(normuu/normu2);
		
		// check convergence
		if (normuu<maxdiff)
		{
			if (stag<MaxStag) {
                if (print_level>=PRINT_MORE) { 
                    ITS_DIFFRES(normuu,relres1);
                    ITS_RESTART;
                }
			}
			
			fasp_array_cp(m,b->val,r); fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
			
			temp2=fasp_blas_array_dotprod(m,r,r);
			absres=sqrt(temp2);
			switch (stop_type) {
				case STOP_REL_RES:
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_REL_PRECRES:
					if (pc == NULL)
						fasp_array_cp(m,r,t);
					else
						pc->fct(r,t,pc->data);
					temp2=ABS(fasp_blas_array_dotprod(m,r,t));
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_MOD_REL_RES:
					relres1=sqrt(temp2/normu2); 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres1);
			
			if (relres1<tol)
				break;
			else
			{
				if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
					iter = ERROR_SOLVER_STAG;
					break;
				}
				++stag;
				++restart_step;
				
				fasp_array_set(m,p0,0.0);
				
				// p1 = B*r
				if (pc == NULL)
					fasp_array_cp(m,r,p1); /* No preconditioner, B=I */
				else
					pc->fct(r,p1,pc->data); /* Preconditioning */
				
				// tp=A*p1
				fasp_array_set(m,tp,0.0); fasp_blas_bdcsr_aAxpy(1.0,A,p1,tp);
				
				// tz = B*tp
				if (pc == NULL)
					fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
				else
					pc->fct(tp,tz,pc->data); /* Preconditioning */
				
				// p1=p1/normp
				normp=fasp_blas_array_dotprod(m,tz,tp);
				normp=sqrt(normp);
				fasp_array_cp(m,p1,t);
				fasp_array_set(m,p1,0.0); fasp_blas_array_axpy(m,1/normp,t,p1);
				
				// t0=A*p0=0
				fasp_array_set(m,t0,0.0);
				fasp_array_cp(m,t0,z0);
				fasp_array_cp(m,t0,t1);
				fasp_array_cp(m,t0,z1);
				
				// t1=tp/normp,z1=tz/normp
				fasp_blas_array_axpy(m,1/normp,tp,t1);
				fasp_blas_array_axpy(m,1/normp,tz,z1);
			}
		}
		
		// safe guard
		if (relres1<tol) 
		{
			if (print_level>=PRINT_MORE) ITS_COMPRES(relres1);
			
			fasp_array_cp(m,b->val,r); fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
			
			temp2=fasp_blas_array_dotprod(m,r,r);
			absres=sqrt(temp2);
			switch (stop_type) {
				case STOP_REL_PRECRES:
					if (pc == NULL)
						fasp_array_cp(m,r,t);
					else
						pc->fct(r,t,pc->data);
					temp2=ABS(fasp_blas_array_dotprod(m,r,t));
					relres1=sqrt(temp2/normb2); 
					break;
				case STOP_MOD_REL_RES:
					relres1=sqrt(temp2/normu2); 
					break;
				default: // STOP_REL_RES
					relres1=sqrt(temp2/normb2); 
					break;
			}
			
			if (print_level>=PRINT_MORE) ITS_REALRES(relres1);
			
			// check convergence
			if (relres1<tol) break;
			
			if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
				iter = ERROR_SOLVER_TOLSMALL;
				break;
			}
			
			if (more_step<MaxRestartStep) {
				if (print_level>PRINT_NONE) ITS_RESTART;
			}
			
			++more_step;
			++restart_step;
			
			fasp_array_set(m,p0,0.0);
			
			// p1 = B*r
			if (pc == NULL)
				fasp_array_cp(m,r,p1); /* No preconditioner, B=I */
			else
				pc->fct(r,p1,pc->data); /* Preconditioning */
			
			// tp = A*p1
			fasp_array_set(m,tp,0.0);
			fasp_blas_bdcsr_aAxpy(1.0,A,p1,tp);
			
			// tz = B*tp
			if (pc == NULL)
				fasp_array_cp(m,tp,tz); /* No preconditioner, B=I */
			else
				pc->fct(tp,tz,pc->data); /* Preconditioning */
			
			// p1 = p1/normp
			normp=fasp_blas_array_dotprod(m,tz,tp);
			normp=sqrt(normp);
			fasp_array_cp(m,p1,t);
			fasp_array_set(m,p1,0.0); fasp_blas_array_axpy(m,1/normp,t,p1);
			
			// t0=A*p0=0
			fasp_array_set(m,t0,0.0);
			fasp_array_cp(m,t0,z0);
			fasp_array_cp(m,t0,t1);
			fasp_array_cp(m,t0,z1);
			
			// t1=tp/normp,z1=tz/normp
			fasp_blas_array_axpy(m,1/normp,tp,t1);
			fasp_blas_array_axpy(m,1/normp,tz,z1);
			
		}
		// update relative residual here
		absres0 = absres;
	}
	
FINISHED: // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres1);
	
	fasp_mem_free(p0);
	fasp_mem_free(p1);
	fasp_mem_free(p2);
	fasp_mem_free(z0);
	fasp_mem_free(z1);
	fasp_mem_free(t0);
	fasp_mem_free(t1);
	fasp_mem_free(tp);
	fasp_mem_free(tz);
	fasp_mem_free(r);
	fasp_mem_free(t);

#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_bdcsr_pminres ...... [Finish]\n");
#endif	
    
	return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
