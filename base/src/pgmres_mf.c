/*! \file pgmres.c
 *  \brief Krylov subspace methods -- Preconditioned GMRes.
 *
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 *  \note Refer to A.H. Baker, E.R. Jessup, and Tz.V. Kolev
 *        A Simple Strategy for Varying the Restart Parameter in GMRES(m)
 *        Journal of Computational and Applied Mathematics, 230 (2009)
 *        pp. 751-761. UCRL-JRNL-235266.
 *
 */  

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/*!
 * \fn INT fasp_solver_pgmres (mxv_matfree *mf, dvector *b, dvector *x, precond *pc, 
 *                             const REAL tol, const INT MaxIt, const SHORT restart,
 *                             const SHORT stop_type, const SHORT print_level)
 *
 * \brief Solve "Ax=b" using PGMRES (right preconditioned) iterative method
 *
 * \param mf           Pointer to the mxv_matfree
 * \param b            Pointer to the dvector of right hand side
 * \param x            Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type -- DOES not support this parameter
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou 
 * \date   2010/11/28
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Feiteng Huang on 09/26/2012, (matrix free)
 */ 
INT fasp_solver_pgmres (mxv_matfree *mf, 
                        dvector *b, 
                        dvector *x, 
                        precond *pc, 
                        const REAL tol,
                        const INT MaxIt, 
                        const SHORT restart,
                        const SHORT stop_type, 
                        const SHORT print_level)
{
    const INT n         = b->row;  
    const INT min_iter  = 0;
    
    // local variables
    INT      iter = 0;
    INT      restartplus1 = restart + 1;
    INT      i, j, k;
    
    REAL     epsmac = SMALLREAL; 
    REAL     r_norm, b_norm, den_norm;
    REAL     epsilon, gamma, t;   
    
    REAL    *c = NULL, *s = NULL, *rs = NULL; 
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL   **p = NULL, **hh = NULL;
    REAL    *work = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_pgmres ...... [Start]\n");
    printf("### DEBUG: maxit = %d, tol = %.4le, stop type = %d\n", MaxIt, tol, stop_type);
#endif    

    /* allocate memory */
    work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
    
    p    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    
    hh   = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
    
    if (print_level>PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;    
    
    for (i = 0; i < restartplus1; i ++) p[i] = s + restart + i*n;
    
    for (i = 0; i < restartplus1; i ++) hh[i] = p[restart] + n + i*restart;
    
    /* initialization */
    mf->fct(mf->data, x->val, p[0]);
    fasp_blas_array_axpby(n, 1.0, b->val, -1.0, p[0]);
    
    b_norm = fasp_blas_array_norm2(n, b->val);
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
    if ( print_level > PRINT_NONE) {
        norms[0] = r_norm;
        if ( print_level >= PRINT_SOME ) {
            ITS_PUTNORM("right-hand side", b_norm);
            ITS_PUTNORM("residual", r_norm);
        }
    }
    
    if (b_norm > 0.0)  den_norm = b_norm;
    else               den_norm = r_norm;
    
    epsilon = tol*den_norm;
    
    /* outer iteration cycle */
    while (iter < MaxIt) {
  
        rs[0] = r_norm;
        if (r_norm == 0.0) {
            fasp_mem_free(work); 
            fasp_mem_free(p); 
            fasp_mem_free(hh);
            fasp_mem_free(norms);
            return iter; 
        }
    
        if (r_norm <= epsilon && iter >= min_iter) {
            mf->fct(mf->data, x->val, r);
            fasp_blas_array_axpby(n, 1.0, b->val, -1.0, r);
            r_norm = fasp_blas_array_norm2(n, r);
    
            if (r_norm <= epsilon) {
                break;
            }
            else {
                if (print_level >= PRINT_SOME) ITS_FACONV;
            }
        }
    
        t = 1.0 / r_norm;
        for (j = 0; j < n; j ++) p[0][j] *= t;
    
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while (i < restart && iter < MaxIt) {

            i ++;  iter ++;
    
            fasp_array_set(n, r, 0.0);
    
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);          
    
            mf->fct(mf->data, r, p[i]);
    
            /* modified Gram_Schmidt */
            for (j = 0; j < i; j ++) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;    
            if (t != 0.0) {
                t = 1.0/t;
                for (j = 0; j < n; j ++) p[i][j] *= t;
            }
    
            for (j = 1; j < i; ++j) {
                t = hh[j-1][i-1];
                hh[j-1][i-1] = s[j-1]*hh[j][i-1] + c[j-1]*t;
                hh[j][i-1] = -s[j-1]*t + c[j-1]*hh[j][i-1];
            }
            t= hh[i][i-1]*hh[i][i-1];
            t+= hh[i-1][i-1]*hh[i-1][i-1];
            gamma = sqrt(t);
            if (gamma == 0.0) gamma = epsmac;
            c[i-1]  = hh[i-1][i-1] / gamma;
            s[i-1]  = hh[i][i-1] / gamma;
            rs[i]   = -s[i-1]*rs[i-1];
            rs[i-1] = c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            r_norm = fabs(rs[i]);
    
            if (print_level > PRINT_NONE) norms[iter] = r_norm;
    
            if (b_norm > 0 ) {
                if (print_level > PRINT_NONE) print_itinfo(print_level,stop_type,iter,norms[iter]/b_norm,norms[iter],norms[iter]/norms[iter-1]);
            }
            else {
                if (print_level > PRINT_NONE) print_itinfo(print_level,stop_type,iter,norms[iter],norms[iter],norms[iter]/norms[iter-1]);
            }
    
            /* should we exit the restart cycle? */
            if (r_norm <= epsilon && iter >= min_iter) {
                break;
            }         
        } /* end of restart cycle */
    
        /* now compute solution, first solve upper triangular system */    
        rs[i-1] = rs[i-1] / hh[i-1][i-1];
        for (k = i-2; k >= 0; k --) {
            t = 0.0;
            for (j = k+1; j < i; j ++)  t -= hh[k][j]*rs[j];
    
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        fasp_array_cp(n, p[i-1], w);
        for (j = 0; j < n; j ++) w[j] *= rs[i-1];
        for (j = i-2; j >= 0; j --)  fasp_blas_array_axpy(n, rs[j], p[j], w);
        fasp_array_set(n, r, 0.0);
    
        /* apply the preconditioner */
        if (pc == NULL)
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
    
        fasp_blas_array_axpy(n, 1.0, r, x->val);
    
        if (r_norm  <= epsilon && iter >= min_iter) {
            mf->fct(mf->data, x->val, r);
            fasp_blas_array_axpby(n, 1.0, b->val, -1.0, r);
            r_norm = fasp_blas_array_norm2(n, r);
    
            if (r_norm  <= epsilon) {
                break;
            }
            else {
                if (print_level >= PRINT_SOME) ITS_FACONV;
                fasp_array_cp(n, r, p[0]); i = 0;
            }
        } /* end of convergence check */
    
        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j] = c[j-1]*rs[j];
        }
    
        if (i) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
    
        for (j = i-1 ; j > 0; j --) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
    
        if (i) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }        
    } /* end of iteration while loop */
    
    if (print_level > PRINT_NONE) ITS_FINAL(iter,MaxIt,r_norm);

    /*-------------------------------------------
     * Clean up workspace
     *------------------------------------------*/
    fasp_mem_free(work); 
    fasp_mem_free(p); 
    fasp_mem_free(hh);
    fasp_mem_free(norms);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_pgmres ...... [Finish]\n");
#endif    

    if (iter>=MaxIt) 
        return ERROR_SOLVER_MAXIT;
    else 
        return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
