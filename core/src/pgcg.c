/*! \file pgcg.c
 *  \brief Krylov subspace methods -- Preconditioned Generalized Conjugate Gradient.
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
 * \fn INT fasp_solver_dcsr_pgcg (dCSRmat *A, dvector *b, dvector *u, const INT MaxIt, 
 *                                const double tol, precond *pc, const SHORT print_level,
 *                                const SHORT stop_type)
 *
 * \brief Preconditioned generilzed conjugate gradient (GCG) method for solving Au=b 
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
 * \author Xiaozhe Hu
 * \date   01/01/2012
 *
 * \note Not completely implemented yet! --Chensong
 */
INT fasp_solver_dcsr_pgcg (dCSRmat *A, 
                           dvector *b, 
                           dvector *u, 
                           const INT MaxIt, 
                           const double tol,
                           precond *pc, 
                           const SHORT print_level,
                           const SHORT stop_type)
{
	INT iter=0, m=A->row, i;
	double absres0=BIGREAL, absres, relres=BIGREAL, factor;
	double alpha, normb=BIGREAL;
	SHORT status = SUCCESS;
    
	// allocate temp memory 
	double *work = (double *)fasp_mem_calloc(2*m+MaxIt+MaxIt*m,sizeof(double));	
	
    double *r, *Br, *beta, *p;
	r = work; Br = r + m; beta = Br + m; p = beta + MaxIt;
	
#if CHMEM_MODE		
	total_alloc_mem += (2*m+MaxIt+MaxIt*m)*sizeof(double);
#endif
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pgcg ...... [Start]\n");
#endif	
    
	normb=fasp_blas_array_norm2(m,b->val);
	
	// -------------------------------------
	// 1st iteration (Steepest descent)
	// -------------------------------------
	// r = b-A*u
	fasp_array_cp(m,b->val,r);
	fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);
	
	//Br 
	if (pc != NULL)
		pc->fct(r,p,pc->data); /* Preconditioning */
	else
		fasp_array_cp(m,r,p); /* No preconditioner, B=I */
	
	// alpha = (p'r)/(p'Ap)
	alpha = fasp_blas_array_dotprod (m,r,p) / fasp_blas_dcsr_vmv (A, p, p);
    
	// u = u + alpha *p
	fasp_blas_array_axpy(m, alpha , p, u->val);
	
    // r = r - alpha *Ap
	fasp_blas_dcsr_aAxpy((-1.0*alpha),A,p,r);
    
    // norm(r), factor
    absres = fasp_blas_array_norm2(m,r); factor = absres/absres0;
    
    // compute relative residual 
    relres = absres/normb;	
    
    // output iteration information if needed	
    print_itinfo(print_level,stop_type,iter+1,relres,absres,factor);
    
    // update relative residual here
    absres0 = absres;
	
	for ( iter = 1; iter < MaxIt ; iter++)
	{		
		// Br
		if (pc != NULL)
			pc->fct(r, Br ,pc->data); // Preconditioning 
		else
			fasp_array_cp(m,r, Br); // No preconditioner, B=I 
        
		// form p
		fasp_array_cp(m, Br, p+iter*m);
		
		for (i=0; i<iter; i++)
		{
			beta[i] = (-1.0) * ( fasp_blas_dcsr_vmv (A, Br, p+i*m)
                                /fasp_blas_dcsr_vmv (A, p+i*m, p+i*m) );
			
			fasp_blas_array_axpy(m, beta[i], p+i*m, p+iter*m);
		}
		
        // -------------------------------------
        // next iteration
        // -------------------------------------

        // alpha = (p'r)/(p'Ap)
		alpha = fasp_blas_array_dotprod(m,r,p+iter*m)
              / fasp_blas_dcsr_vmv (A, p+iter*m, p+iter*m);
		
        // u = u + alpha *p
		fasp_blas_array_axpy(m, alpha , p+iter*m, u->val);
		
        // r = r - alpha *Ap
		fasp_blas_dcsr_aAxpy((-1.0*alpha),A,p+iter*m,r);
		
		// norm(r), factor
		absres = fasp_blas_array_norm2(m,r); factor = absres/absres0;
		
		// compute relative residual 
		relres = absres/normb;	
		
		// output iteration information if needed	
		print_itinfo(print_level,stop_type,iter+1,relres,absres,factor);
		
		if (relres < tol) break;
		
		// update relative residual here
		absres0 = absres;
        
	} // end of main GCG loop.
	
    // finish the iterative method
	if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres);
	
	// clean up temp memory
	fasp_mem_free(work);
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_dcsr_pgcg ...... [Finish]\n");
#endif	
	
	if (iter>MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
