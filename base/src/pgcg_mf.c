/*! \file pgcg.c
 *  \brief Krylov subspace methods -- Preconditioned Generalized Conjugate Gradient.
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
 * \fn INT fasp_solver_pgcg (mxv_matfree *mf, dvector *b, dvector *u, precond *pc, 
 *                           const REAL tol, const INT MaxIt, 
 *                           const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned generilzed conjugate gradient (GCG) method for solving Au=b 
 *
 * \param mf           Pointer to the mxv_matfree
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type -- Not implemented
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
 * Modified by Feiteng Huang on 09/26/2012, (matrix free)
 */
INT fasp_solver_pgcg (mxv_matfree *mf, 
                      dvector *b, 
                      dvector *u, 
                      precond *pc, 
                      const REAL tol,
                      const INT MaxIt, 
                      const SHORT stop_type,
                      const SHORT print_level)
{
    INT    iter=0, m=b->row, i;
    REAL   absres0=BIGREAL, absres, relres=BIGREAL, factor;
    REAL   alpha, normb=BIGREAL, gama_1, gama_2;
    
    // allocate temp memory 
    REAL *work = (REAL *)fasp_mem_calloc(3*m+MaxIt+MaxIt*m,sizeof(REAL));    
    
    REAL *r, *Br, *beta, *p, *q;
    q = work; r = q + m; Br = r + m; beta = Br + m; p = beta + MaxIt;
        
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_pgcg ...... [Start]\n");
    printf("### DEBUG: maxit = %d, tol = %.4le, stop type = %d\n", MaxIt, tol, stop_type);
#endif    

    normb=fasp_blas_array_norm2(m,b->val);
    
    // -------------------------------------
    // 1st iteration (Steepest descent)
    // -------------------------------------
    // r = b-A*u
    mf->fct(mf->data, u->val, r);
    fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
    
    // Br 
    if (pc != NULL)
        pc->fct(r,p,pc->data); /* Preconditioning */
    else
        fasp_array_cp(m,r,p); /* No preconditioner, B=I */
    
    // alpha = (p'r)/(p'Ap)
    mf->fct(mf->data, p, q);
    alpha = fasp_blas_array_dotprod (m,r,p) / fasp_blas_array_dotprod (m, p, q);
    
    // u = u + alpha *p
    fasp_blas_array_axpy(m, alpha , p, u->val);
    
    // r = r - alpha *Ap
    mf->fct(mf->data, p, q);
    fasp_blas_array_axpby(m, (-1.0*alpha), q, 1.0, r);
    
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
            mf->fct(mf->data, Br, q);
            gama_1 = fasp_blas_array_dotprod(m, p+i*m, q);
            mf->fct(mf->data, p+i*m, q);
            gama_2 = fasp_blas_array_dotprod(m, p+i*m, q);
            beta[i] = (-1.0) * ( gama_1 / gama_2 );
    
            fasp_blas_array_axpy(m, beta[i], p+i*m, p+iter*m);
        }
    
        // -------------------------------------
        // next iteration
        // -------------------------------------

        // alpha = (p'r)/(p'Ap)
        mf->fct(mf->data, p+iter*m, q);
        alpha = fasp_blas_array_dotprod(m,r,p+iter*m)
            / fasp_blas_array_dotprod (m, q, p+iter*m);
    
        // u = u + alpha *p
        fasp_blas_array_axpy(m, alpha , p+iter*m, u->val);
    
        // r = r - alpha *Ap
        mf->fct(mf->data, p+iter*m, q);
        fasp_blas_array_axpby(m, (-1.0*alpha), q, 1.0, r);
    
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
    
    if (iter>MaxIt) 
        return ERROR_SOLVER_MAXIT;
    else 
        return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
