/*! \file pgmres_mf.c
 *
 *  \brief Krylov subspace methods -- Preconditioned GMRes (matrix free)
 *
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 *  \note Refer to A.H. Baker, E.R. Jessup, and Tz.V. Kolev
 *        A Simple Strategy for Varying the Restart Parameter in GMRES(m)
 *        Journal of Computational and Applied Mathematics, 230 (2009)
 *        pp. 751-761. UCRL-JRNL-235266.
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
 *                             const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Solve "Ax=b" using PGMRES (right preconditioned) iterative method
 *
 * \param mf           Pointer to mxv_matfree: the spmv operation
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type -- DOES not support this parameter
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Zhiyang Zhou
 * \date   2010/11/28
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Feiteng Huang on 09/26/2012: matrix free
 * Modified by Chunsheng Feng on 07/22/2013: Add adapt memory allocate
 */
INT fasp_solver_pgmres (mxv_matfree *mf,
                        dvector *b,
                        dvector *x,
                        precond *pc,
                        const REAL tol,
                        const INT MaxIt,
                        const SHORT restart,
                        const SHORT stop_type,
                        const SHORT prtlvl)
{
    const INT n         = b->row;
    const INT min_iter  = 0;
    
    // local variables
    INT      iter = 0;
    INT      i, j, k;
    
    REAL     epsmac = SMALLREAL;
    REAL     r_norm, b_norm, den_norm;
    REAL     epsilon, gamma, t;
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL    *c = NULL, *s = NULL, *rs = NULL;
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL    *work = NULL;
    REAL    **p = NULL, **hh = NULL;
    
    unsigned INT  Restart  = restart;
    unsigned INT  Restart1 = Restart + 1;
    unsigned LONG worksize = (Restart+4)*(Restart+n)+1-n;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc(worksize, sizeof(REAL));
    
    /* check whether memory is enough for GMRES */
    while ( (work == NULL) && (Restart > 5) ) {
        Restart = Restart - 5;
        worksize = (Restart+4)*(Restart+n)+1-n;
        work = (REAL *) fasp_mem_calloc(worksize, sizeof(REAL));
        Restart1 = Restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }
    
    if ( prtlvl > PRINT_MIN && Restart < restart ) {
        printf("### WARNING: GMRES restart number set to %d!\n", Restart);
    }
    
    p     = (REAL **)fasp_mem_calloc(Restart1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(Restart1, sizeof(REAL *));
    norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c  = rs + Restart1; s  = c + Restart;
    
    for (i = 0; i < Restart1; i ++) p[i] = s + Restart + i*n;
    
    for (i = 0; i < Restart1; i ++) hh[i] = p[Restart] + n + i*Restart;
    
    /* initialization */
    mf->fct(mf->data, x->val, p[0]);
    fasp_blas_array_axpby(n, 1.0, b->val, -1.0, p[0]);
    
    b_norm = fasp_blas_array_norm2(n, b->val);
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
    if ( prtlvl > PRINT_NONE) {
        norms[0] = r_norm;
        if ( prtlvl >= PRINT_SOME ) {
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
                if (prtlvl >= PRINT_SOME) ITS_FACONV;
            }
        }
        
        t = 1.0 / r_norm;
        //for (j = 0; j < n; j ++) p[0][j] *= t;
        fasp_blas_array_ax(n, t, p[0]);
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while (i < Restart && iter < MaxIt) {
            
            i ++;  iter ++;
            
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
                //for (j = 0; j < n; j ++) p[i][j] *= t;
                fasp_blas_array_ax(n, t, p[i]);
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
            
            norms[iter] = r_norm;
            
            if (b_norm > 0 ) {
                print_itinfo(prtlvl,stop_type,iter,norms[iter]/b_norm,
                             norms[iter],norms[iter]/norms[iter-1]);
            }
            else {
                print_itinfo(prtlvl,stop_type,iter,norms[iter],norms[iter],
                             norms[iter]/norms[iter-1]);
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
        //for (j = 0; j < n; j ++) w[j] *= rs[i-1];
        fasp_blas_array_ax(n, rs[i-1], w);
        for (j = i-2; j >= 0; j --)  fasp_blas_array_axpy(n, rs[j], p[j], w);
        
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
                if (prtlvl >= PRINT_SOME) ITS_FACONV;
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
    
    if (prtlvl > PRINT_NONE) ITS_FINAL(iter,MaxIt,r_norm);
    
    /*-------------------------------------------
     * Clean up workspace
     *------------------------------------------*/
    fasp_mem_free(work);
    fasp_mem_free(p);
    fasp_mem_free(hh);
    fasp_mem_free(norms);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if (iter>=MaxIt)
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
