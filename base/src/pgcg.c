/*! \file pgcg.c
 *
 *  \brief Krylov subspace methods -- Preconditioned Generalized CG
 *
 *  \note Refer to Concus, P. and Golub, G.H. and O'Leary, D.P.
 *        A Generalized Conjugate Gradient Method for the Numerical: 
 *        Solution of Elliptic Partial Differential Equations,
 *        Computer Science Department, Stanford University, 1976
 */  

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dcsr_pgcg (dCSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                const REAL tol, const INT MaxIt, 
 *                                const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned generilzed conjugate gradient (GCG) method for solving Au=b 
 *
 * \param A            Pointer to dCSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param u            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu
 * \date   01/01/2012
 *
 * \note Not completely implemented yet! --Chensong
 *
 * Modified by Chensong Zhang on 05/01/2012
 */
INT fasp_solver_dcsr_pgcg (dCSRmat *A, 
                           dvector *b, 
                           dvector *u, 
                           precond *pc, 
                           const REAL tol,
                           const INT MaxIt, 
                           const SHORT stop_type,
                           const SHORT print_level)
{
    INT    iter=0, m=A->row, i;
    REAL   absres0 = BIGREAL, absres = BIGREAL;
    REAL   relres  = BIGREAL, normb  = BIGREAL;
    REAL   alpha, factor;

    
    // allocate temp memory 
    REAL *work = (REAL *)fasp_mem_calloc(2*m+MaxIt+MaxIt*m,sizeof(REAL));    
    
    REAL *r, *Br, *beta, *p;
    r = work; Br = r + m; beta = Br + m; p = beta + MaxIt;
        
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif

    normb=fasp_blas_array_norm2(m,b->val);
    
    // -------------------------------------
    // 1st iteration (Steepest descent)
    // -------------------------------------
    // r = b-A*u
    fasp_array_cp(m,b->val,r);
    fasp_blas_dcsr_aAxpy(-1.0,A,u->val,r);

    // Br 
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
    
    for ( iter = 1; iter < MaxIt ; iter++) {
    
        // Br
        if (pc != NULL)
            pc->fct(r, Br ,pc->data); // Preconditioning 
        else
            fasp_array_cp(m,r, Br); // No preconditioner, B=I 
        
        // form p
        fasp_array_cp(m, Br, p+iter*m);
    
        for (i=0; i<iter; i++) {
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
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if (iter>MaxIt) 
        return ERROR_SOLVER_MAXIT;
    else 
        return iter;
}

/**
 *	\fn int fasp_krylov_cycle_dcsr_pgcg (dCSRmat *A,  dvector *b, dvector *u, precond *pc)
 *
 *	\brief A preconditioned GCR method for solving Au=b 
 *
 *	\param *A	 Pointer to the coefficient matrix
 *	\param *b	 Pointer to the dvector of right hand side
 *	\param *u	 Pointer to the dvector of dofs
 *	\param *pre  Pointer to the structure of precondition (precond) 
 *
 * \author zheng Li, Chensong Zhang 
 * \date   11/09/2014
 *
 * \Note: Specified for unsmoothed aggreagtion cycle.
 */
INT fasp_krylov_cycle_dcsr_pgcg (dCSRmat *A, 
                                 dvector *b, 
                                 dvector *u, 
                                 precond *pc) 
{
    REAL   absres, relres, normb;
    REAL   alpha1, alpha2, gamma1, gamma2, rho1, rho2, beta1, beta2, beta3, beta4; 
    REAL   *work, *r, *x1, *v1, *v2;

    INT    m=A->row;
	REAL   *x = u->val;
    
    // allocate temp memory 
    work = (REAL *)fasp_mem_calloc(4*m,sizeof(REAL));    
    r = work; x1 = r + m; v1 = r + 2*m; v2 = r + 3*m;

    normb = fasp_blas_array_norm2(m, b->val);

	fasp_array_cp(m, b->val, r);

    // Preconditioning 
    if (pc != NULL)
        pc->fct(r, x, pc->data); 
    else
        fasp_array_cp(m, r, x);
    
	//v1 = A*p
	fasp_blas_dcsr_mxv(A, x, v1);

	// rho1 = (p,v1)
	rho1 = fasp_blas_array_dotprod (m, x, v1);

    // alpha1 = (p, r)
	alpha1 = fasp_blas_array_dotprod (m, x, r);

	beta1 = alpha1/rho1;

    // r = r - beta1 *v1
    fasp_blas_array_axpy(m, -beta1, v1, r);

    // norm(r)
    absres = fasp_blas_array_norm2(m, r);

    // compute relative residual 
    relres = absres/normb;    

    // if relres reachs tol(0.2), pgcr will stop,
    // otherwise, another one pgcr iteration will do.
    if(relres < 0.2) {
        fasp_blas_array_ax(m, beta1, x);
		return 0;
	}

    // Preconditioning 
    if (pc != NULL)
        pc->fct(r, x1, pc->data);
    else
        fasp_array_cp(m, r, x1);

	//v2 = A*p
	fasp_blas_dcsr_mxv(A, x1, v2);

	//gamma0 = (x1,v1)
	gamma1 = fasp_blas_array_dotprod (m, x1, v1);

	//alpha2 = (x1,r)
	alpha2  = fasp_blas_array_dotprod(m, x1, r); 

	//rho2 = (x1,v2)
	rho2 = fasp_blas_array_dotprod(m, x1, v2);

	gamma2 = gamma1;

	beta2 = rho2 - gamma1*gamma2/rho1;
	beta3 = (alpha1 - gamma2*alpha2/beta2)/rho1;
    beta4 = alpha2/beta2;

    fasp_blas_array_ax(m, beta3, x);
    
	fasp_blas_array_axpy(m, beta4, x1, x);

    // free 
    fasp_mem_free(work);
}
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
