/*! \file pvfgmres_mf.c
 *
 *  \brief Krylov subspace methods -- Preconditioned variable-restarting
 *         flexible GMRes (matrix free)
 *
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 *  \note Refer to A.H. Baker, E.R. Jessup, and Tz.V. Kolev
 *        A Simple Strategy for Varying the Restart Parameter in GMRES(m)
 *        Journal of Computational and Applied Mathematics, 230 (2009)
 *        pp. 751-761. UCRL-JRNL-235266.
 *
 *  \note This file is modifed from pvgmres.c
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
 * \fn INT fasp_solver_pvfgmres (mxv_matfree *mf, dvector *b, dvector *x, precond *pc, 
 *                               const REAL tol, const INT MaxIt, const SHORT restart,
 *                               const SHORT stop_type, const SHORT print_level)
 *
 * \brief Solve "Ax=b" using PFGMRES(right preconditioned) iterative method in which 
 *        the restart parameter can be adaptively modified during the iteration and 
 *        flexible preconditioner can be used.
 *
 * \param mf           Pointer to mxv_matfree: the spmv operation
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type -- DOES not support this parameter
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu 
 * \date   01/04/2012
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Feiteng Huang on 09/26/2012: matrix free
 * Modified by Chunsheng Feng on 07/22/2013: Add adapt memory allocate
 */ 
INT fasp_solver_pvfgmres (mxv_matfree *mf, 
                          dvector *b, 
                          dvector *x, 
                          precond *pc, 
                          const REAL tol,
                          const INT MaxIt, 
                          const SHORT restart,
                          const SHORT stop_type, 
                          const SHORT print_level)
{
    const INT n                 = b->row;  
    const INT min_iter          = 0;
    
    //--------------------------------------------//
    //   Newly added parameters to monitor when   //
    //   to change the restart parameter          //
    //--------------------------------------------//    
    const REAL cr_max           = 0.99;    // = cos(8^o)  (experimental) 
    const REAL cr_min           = 0.174;   // = cos(80^o) (experimental)
    
    // local variables
    INT    iter                 = 0;
    INT    i,j,k;
    
    REAL   epsmac               = SMALLREAL;
    REAL   r_norm, b_norm, den_norm;
    REAL   epsilon, gamma, t;
    
    REAL  *c = NULL, *s = NULL, *rs = NULL;
    REAL  *norms = NULL, *r = NULL;
    REAL  **p = NULL, **hh = NULL, **z=NULL;
    REAL  *work = NULL;
    
    REAL   cr          = 1.0;     // convergence rate
    REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
    INT    d           = 3;       // reduction for the restart parameter
    INT    restart_max = restart; // upper bound for restart in each restart cycle
    INT    restart_min = 3;       // lower bound for restart in each restart cycle (should be small)
    INT    Restart = restart;     // the real restart in some fixed restarted cycle
    INT    Restartplus1 = Restart + 1;
    LONG   worksize = (Restart+4)*(Restart+n)+1-n+Restart*n;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc(worksize, sizeof(REAL));
    
    while ( (work == NULL) && (Restart > 5 ) ) {
        Restart = Restart - 5;
        worksize = (Restart+4)*(Restart+n)+1-n+Restart*n;
        work = (REAL *) fasp_mem_calloc(worksize, sizeof(REAL));
        printf("### WARNING: vFGMRES restart number changed to %d!\n", Restart);
        Restartplus1 = Restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for vFGMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }
    
    p  = (REAL **)fasp_mem_calloc(Restartplus1, sizeof(REAL *));
    hh = (REAL **)fasp_mem_calloc(Restartplus1, sizeof(REAL *));
    z  = (REAL **)fasp_mem_calloc(Restartplus1, sizeof(REAL *));
    
    if (print_level > PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
    
    r = work; rs = r + n; c = rs + Restartplus1; s = c + Restart;    
    for (i = 0; i < Restartplus1; i ++) p[i] = s + Restart + i*n;
    for (i = 0; i < Restartplus1; i ++) hh[i] = p[Restart] + n + i*Restart;
    for (i = 0; i < Restartplus1; i ++) z[i] = hh[Restart] + Restart + i*n;
    
    /* initialization */
    mf->fct(mf->data, x->val, p[0]);
    fasp_blas_array_axpby(n, 1.0, b->val, -1.0, p[0]);
    
    b_norm = fasp_blas_array_norm2(n, b->val);
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
    if ( print_level > PRINT_NONE) {
        norms[0] = r_norm;
        if ( print_level >= PRINT_SOME) {
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
        r_norm_old = r_norm;
        if (r_norm == 0.0) {
            fasp_mem_free(work); 
            fasp_mem_free(p); 
            fasp_mem_free(hh);
            fasp_mem_free(norms);
            fasp_mem_free(z);
            return iter; 
        }
    
        //-----------------------------------//
        //   adjust the restart parameter    //
        //-----------------------------------//
    
        if (cr > cr_max || iter == 0) {
            Restart = restart_max;
        }
        else if (cr < cr_min) {
            Restart = Restart;
        }
        else {
            if (Restart - d > restart_min) {
                Restart -= d;
            }
            else {
                Restart = restart_max;
            }
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
        //for (j = 0; j < n; j ++) p[0][j] *= t;
        fasp_blas_array_ax(n, t, p[0]);
    
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while (i < Restart && iter < MaxIt) {
            i ++;  iter ++;
    
            fasp_array_set(n, z[i-1], 0.0);
            
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], z[i-1]);
            else
                pc->fct(p[i-1], z[i-1], pc->data);
    
            mf->fct(mf->data, z[i-1], p[i]);
    
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
    
            if (print_level > PRINT_NONE) norms[iter] = r_norm;
    
            if (b_norm > 0 ) {
                if (print_level > PRINT_NONE) 
                    print_itinfo(print_level,stop_type,iter,norms[iter]/b_norm,
                                 norms[iter],norms[iter]/norms[iter-1]);
            }
            else {
                if (print_level > PRINT_NONE) 
                    print_itinfo(print_level,stop_type,iter,norms[iter],norms[iter],
                                 norms[iter]/norms[iter-1]);
            }
    
            /* should we exit the restart cycle? */
            if (r_norm <= epsilon && iter >= min_iter) break;
            
        } /* end of restart cycle */
    
        /* now compute solution, first solve upper triangular system */
    
        rs[i-1] = rs[i-1] / hh[i-1][i-1];
        for (k = i-2; k >= 0; k --) {
            t = 0.0;
            for (j = k+1; j < i; j ++)  t -= hh[k][j]*rs[j];
    
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        
        fasp_array_cp(n, z[i-1], r);
        //for (j = 0; j < n; j ++) r[j] *= rs[i-1];
        fasp_blas_array_ax(n, rs[i-1], r);
        for (j = i-2; j >= 0; j --)  fasp_blas_array_axpy(n, rs[j], z[j], r);
        
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
    
        //-----------------------------------//
        //   compute the convergence rate    //
        //-----------------------------------//    
        cr = r_norm / r_norm_old;
    
    } /* end of iteration while loop */

    if (print_level > PRINT_NONE) ITS_FINAL(iter,MaxIt,r_norm);
    
    /*-------------------------------------------
     * Free some stuff
     *------------------------------------------*/
    fasp_mem_free(work); 
    fasp_mem_free(p); 
    fasp_mem_free(hh);
    fasp_mem_free(norms);
    fasp_mem_free(z);
    
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
