/*! \file pgmres.c
 *
 *  \brief Krylov subspace methods -- Preconditioned GMRes
 *
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 *  \note See also pvgmres.c for a variable restarting version.
 *
 *  \note See spgmres.c for a safer version
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/*!
 * \fn INT fasp_solver_dcsr_pgmres (dCSRmat *A, dvector *b, dvector *x, precond *pc,
 *                                  const REAL tol, const INT MaxIt, SHORT restart,
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned GMRES method for solving Au=b
 *
 * \param A            Pointer to dCSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou
 * \date   2010/11/28
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Chensong Zhang on 04/05/2013: add stop_type and safe check
 * Modified by Chunsheng Feng on 07/22/2013: Add adapt memory allocate
 */
INT fasp_solver_dcsr_pgmres (dCSRmat *A,
                             dvector *b,
                             dvector *x,
                             precond *pc,
                             const REAL tol,
                             const INT MaxIt,
                             SHORT restart,
                             const SHORT stop_type,
                             const SHORT print_level)
{
    const INT   n         = b->row;
    const INT   MIN_ITER  = 0;
    const REAL  epsmac    = SMALLREAL;
    
    // local variables
    INT      iter = 0;
    INT      restartplus1 = restart + 1;
    INT      i, j, k;
    
    REAL     r_norm, r_normb, gamma, t;
    REAL     absres0, absres, relres, normu;
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL    *c = NULL, *s = NULL, *rs = NULL;
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL    *work = NULL;
    REAL    **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
    
    while ( (work == NULL) && (restart > 5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }
    
    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;
    
    for ( i = 0; i < restartplus1; i++ ) p[i]  = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dcsr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm  = fasp_blas_array_norm2(n,p[0]);
    
    // compute initial residuals
    switch (stop_type) {
        case STOP_REL_RES:
            absres0 = MAX(SMALLREAL,r_norm);
            relres  = r_norm/absres0;
            break;
        case STOP_REL_PRECRES:
            if ( pc == NULL )
                fasp_array_cp(n, p[0], r);
            else
                pc->fct(p[0], r, pc->data);
            r_normb = sqrt(fasp_blas_array_dotprod(n,p[0],r));
            absres0 = MAX(SMALLREAL,r_normb);
            relres  = r_normb/absres0;
            break;
        case STOP_MOD_REL_RES:
            normu   = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
            absres0 = r_norm;
            relres  = absres0/normu;
            break;
        default:
            printf("### WARNING: Unrecognized stopping type!\n");
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,0,relres,absres0,0.0);
    
    // store initial residual
    norms[0] = relres;
    
    /* outer iteration cycle */
    while ( iter < MaxIt ) {
        
        rs[0] = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < restart && iter < MaxIt ) {
            
            i++; iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if ( pc == NULL )
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_dcsr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for ( j = 0; j < i; j++ ) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;
            if (t != 0.0) {
                t = 1.0/t;
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
            rs[i-1] =  c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
        } /* end of restart cycle */
        
        /* compute solution, first solve upper triangular system */
        rs[i-1] = rs[i-1] / hh[i-1][i-1];
        for ( k = i-2; k >= 0; k-- ) {
            t = 0.0;
            for (j = k+1; j < i; j++) t -= hh[k][j]*rs[j];
            
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        
        fasp_array_cp(n, p[i-1], w);
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // Check: prevent false convergence
        if ( relres <= tol && iter >= MIN_ITER ) {
            
            fasp_array_cp(n, b->val, r);
            fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
            
            r_norm = fasp_blas_array_norm2(n, r);
            
            switch ( stop_type ) {
                case STOP_REL_RES:
                    absres = r_norm;
                    relres = absres/absres0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc == NULL )
                        fasp_array_cp(n, r, w);
                    else
                        pc->fct(r, w, pc->data);
                    absres = sqrt(fasp_blas_array_dotprod(n,w,r));
                    relres = absres/absres0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = r_norm;
                    normu  = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            if ( relres <= tol ) {
                break;
            }
            else {
                // Need to restart
                fasp_array_cp(n, r, p[0]); i = 0;
            }
            
        } /* end of convergence check */
        
        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
    } /* end of main while loop */
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
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
    
    if ( iter >= MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/**
 * \fn INT fasp_solver_bdcsr_pgmres (block_dCSRmat *A, dvector *b, dvector *x, precond *pc,
 *                                   const REAL tol, const INT MaxIt, SHORT restart,
 *                                   const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned GMRES method for solving Au=b
 *
 * \param A            Pointer to block_dCSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/24/2010
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Chensong Zhang on 04/05/2013: add stop_type and safe check
 */
INT fasp_solver_bdcsr_pgmres (block_dCSRmat *A,
                              dvector *b,
                              dvector *x,
                              precond *pc,
                              const REAL tol,
                              const INT MaxIt,
                              SHORT restart,
                              const SHORT stop_type,
                              const SHORT print_level)
{
    const INT   n         = b->row;
    const INT   MIN_ITER  = 0;
    const REAL  epsmac    = SMALLREAL;
    
    // local variables
    INT      iter = 0;
    INT      restartplus1 = restart + 1;
    INT      i, j, k;
    
    REAL     r_norm, r_normb, gamma, t;
    REAL     absres0, absres, relres, normu;
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL    *c = NULL, *s = NULL, *rs = NULL;
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL    *work = NULL;
    REAL    **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));

    while ( (work == NULL) && (restart > 5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;
    
    for ( i = 0; i < restartplus1; i++ ) p[i]  = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_bdcsr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm  = fasp_blas_array_norm2(n,p[0]);
    
    // compute initial residuals
    switch (stop_type) {
        case STOP_REL_RES:
            absres0 = MAX(SMALLREAL,r_norm);
            relres  = r_norm/absres0;
            break;
        case STOP_REL_PRECRES:
            if ( pc == NULL )
                fasp_array_cp(n, p[0], r);
            else
                pc->fct(p[0], r, pc->data);
            r_normb = sqrt(fasp_blas_array_dotprod(n,p[0],r));
            absres0 = MAX(SMALLREAL,r_normb);
            relres  = r_normb/absres0;
            break;
        case STOP_MOD_REL_RES:
            normu   = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
            absres0 = r_norm;
            relres  = absres0/normu;
            break;
        default:
            printf("### WARNING: Unrecognized stopping type!\n");
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,0,relres,absres0,0.0);
    
    // store initial residual
    norms[0] = relres;
    
    /* outer iteration cycle */
    while ( iter < MaxIt ) {
        
        rs[0] = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < restart && iter < MaxIt ) {
            
            i++; iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if ( pc == NULL )
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_bdcsr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for ( j = 0; j < i; j++ ) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;
            if (t != 0.0) {
                t = 1.0/t;
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
            rs[i-1] =  c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
        } /* end of restart cycle */
        
        /* compute solution, first solve upper triangular system */
        rs[i-1] = rs[i-1] / hh[i-1][i-1];
        for ( k = i-2; k >= 0; k-- ) {
            t = 0.0;
            for (j = k+1; j < i; j++) t -= hh[k][j]*rs[j];
            
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        
        fasp_array_cp(n, p[i-1], w);
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // Check: prevent false convergence
        if ( relres <= tol && iter >= MIN_ITER ) {
            
            fasp_array_cp(n, b->val, r);
            fasp_blas_bdcsr_aAxpy(-1.0, A, x->val, r);
            
            r_norm = fasp_blas_array_norm2(n, r);
            
            switch ( stop_type ) {
                case STOP_REL_RES:
                    absres = r_norm;
                    relres = absres/absres0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc == NULL )
                        fasp_array_cp(n, r, w);
                    else
                        pc->fct(r, w, pc->data);
                    absres = sqrt(fasp_blas_array_dotprod(n,w,r));
                    relres = absres/absres0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = r_norm;
                    normu  = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            if ( relres <= tol ) {
                break;
            }
            else {
                // Need to restart
                fasp_array_cp(n, r, p[0]); i = 0;
            }
            
        } /* end of convergence check */
        
        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
    } /* end of main while loop */
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
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
    
    if ( iter >= MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/*!
 * \fn INT fasp_solver_dbsr_pgmres (dBSRmat *A, dvector *b, dvector *x, precond *pc,
 *                                  const REAL tol, const INT MaxIt, SHORT restart,
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned GMRES method for solving Au=b
 *
 * \param A            Pointer to dBSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou
 * \date   2010/12/21
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Chensong Zhang on 04/05/2013: add stop_type and safe check
 */
INT fasp_solver_dbsr_pgmres (dBSRmat *A,
                             dvector *b,
                             dvector *x,
                             precond *pc,
                             const REAL tol,
                             const INT MaxIt,
                             SHORT restart,
                             const SHORT stop_type,
                             const SHORT print_level)
{
    const INT   n         = b->row;
    const INT   MIN_ITER  = 0;
    const REAL  epsmac    = SMALLREAL;
    
    // local variables
    INT      iter = 0;
    INT      restartplus1 = restart + 1;
    INT      i, j, k;
    
    REAL     r_norm, r_normb, gamma, t;
    REAL     absres0, absres, relres, normu;
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL    *c = NULL, *s = NULL, *rs = NULL;
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL    *work = NULL;
    REAL    **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));

    while ( (work == NULL) && (restart > 5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;
    
    for ( i = 0; i < restartplus1; i++ ) p[i]  = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dbsr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm  = fasp_blas_array_norm2(n,p[0]);
    
    // compute initial residuals
    switch (stop_type) {
        case STOP_REL_RES:
            absres0 = MAX(SMALLREAL,r_norm);
            relres  = r_norm/absres0;
            break;
        case STOP_REL_PRECRES:
            if ( pc == NULL )
                fasp_array_cp(n, p[0], r);
            else
                pc->fct(p[0], r, pc->data);
            r_normb = sqrt(fasp_blas_array_dotprod(n,p[0],r));
            absres0 = MAX(SMALLREAL,r_normb);
            relres  = r_normb/absres0;
            break;
        case STOP_MOD_REL_RES:
            normu   = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
            absres0 = r_norm;
            relres  = absres0/normu;
            break;
        default:
            printf("### WARNING: Unrecognized stopping type!\n");
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,0,relres,absres0,0.0);
    
    // store initial residual
    norms[0] = relres;
    
    /* outer iteration cycle */
    while ( iter < MaxIt ) {
        
        rs[0] = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < restart && iter < MaxIt ) {
            
            i++; iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if ( pc == NULL )
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_dbsr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for ( j = 0; j < i; j++ ) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;
            if (t != 0.0) {
                t = 1.0/t;
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
            rs[i-1] =  c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
        } /* end of restart cycle */
        
        /* compute solution, first solve upper triangular system */
        rs[i-1] = rs[i-1] / hh[i-1][i-1];
        for ( k = i-2; k >= 0; k-- ) {
            t = 0.0;
            for (j = k+1; j < i; j++) t -= hh[k][j]*rs[j];
            
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        
        fasp_array_cp(n, p[i-1], w);
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // Check: prevent false convergence
        if ( relres <= tol && iter >= MIN_ITER ) {
            
            fasp_array_cp(n, b->val, r);
            fasp_blas_dbsr_aAxpy(-1.0, A, x->val, r);
            
            r_norm = fasp_blas_array_norm2(n, r);
            
            switch ( stop_type ) {
                case STOP_REL_RES:
                    absres = r_norm;
                    relres = absres/absres0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc == NULL )
                        fasp_array_cp(n, r, w);
                    else
                        pc->fct(r, w, pc->data);
                    absres = sqrt(fasp_blas_array_dotprod(n,w,r));
                    relres = absres/absres0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = r_norm;
                    normu  = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            if ( relres <= tol ) {
                break;
            }
            else {
                // Need to restart
                fasp_array_cp(n, r, p[0]); i = 0;
            }
            
        } /* end of convergence check */
        
        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
    } /* end of main while loop */
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
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
    
    if ( iter >= MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/*!
 * \fn INT fasp_solver_dstr_pgmres (dSTRmat *A, dvector *b, dvector *x, precond *pc,
 *                                  const REAL tol, const INT MaxIt, SHORT restart,
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned GMRES method for solving Au=b
 *
 * \param A            Pointer to dSTRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou
 * \date   2010/11/28
 *
 * Modified by Chensong Zhang on 05/01/2012
 * Modified by Chensong Zhang on 04/05/2013: add stop_type and safe check
 */
INT fasp_solver_dstr_pgmres (dSTRmat *A,
                             dvector *b,
                             dvector *x,
                             precond *pc,
                             const REAL tol,
                             const INT MaxIt,
                             SHORT restart,
                             const SHORT stop_type,
                             const SHORT print_level)
{
    const INT   n         = b->row;
    const INT   MIN_ITER  = 0;
    const REAL  epsmac    = SMALLREAL;
    
    // local variables
    INT      iter = 0;
    INT      restartplus1 = restart + 1;
    INT      i, j, k;
    
    REAL     r_norm, r_normb, gamma, t;
    REAL     absres0, absres, relres, normu;
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL    *c = NULL, *s = NULL, *rs = NULL;
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL    *work = NULL;
    REAL    **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
    
    while ( (work == NULL) && (restart > 5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;
    
    for ( i = 0; i < restartplus1; i++ ) p[i]  = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dstr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm  = fasp_blas_array_norm2(n,p[0]);
    
    // compute initial residuals
    switch (stop_type) {
        case STOP_REL_RES:
            absres0 = MAX(SMALLREAL,r_norm);
            relres  = r_norm/absres0;
            break;
        case STOP_REL_PRECRES:
            if ( pc == NULL )
                fasp_array_cp(n, p[0], r);
            else
                pc->fct(p[0], r, pc->data);
            r_normb = sqrt(fasp_blas_array_dotprod(n,p[0],r));
            absres0 = MAX(SMALLREAL,r_normb);
            relres  = r_normb/absres0;
            break;
        case STOP_MOD_REL_RES:
            normu   = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
            absres0 = r_norm;
            relres  = absres0/normu;
            break;
        default:
            printf("### WARNING: Unrecognized stopping type!\n");
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,0,relres,absres0,0.0);
    
    // store initial residual
    norms[0] = relres;
    
    /* outer iteration cycle */
    while ( iter < MaxIt ) {
        
        rs[0] = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < restart && iter < MaxIt ) {
            
            i++; iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if ( pc == NULL )
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_dstr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for ( j = 0; j < i; j++ ) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;
            if (t != 0.0) {
                t = 1.0/t;
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
            rs[i-1] =  c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
        } /* end of restart cycle */
        
        /* compute solution, first solve upper triangular system */
        rs[i-1] = rs[i-1] / hh[i-1][i-1];
        for ( k = i-2; k >= 0; k-- ) {
            t = 0.0;
            for (j = k+1; j < i; j++) t -= hh[k][j]*rs[j];
            
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        
        fasp_array_cp(n, p[i-1], w);
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // Check: prevent false convergence
        if ( relres <= tol && iter >= MIN_ITER ) {
            
            fasp_array_cp(n, b->val, r);
            fasp_blas_dstr_aAxpy(-1.0, A, x->val, r);
            
            r_norm = fasp_blas_array_norm2(n, r);
            
            switch ( stop_type ) {
                case STOP_REL_RES:
                    absres = r_norm;
                    relres = absres/absres0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc == NULL )
                        fasp_array_cp(n, r, w);
                    else
                        pc->fct(r, w, pc->data);
                    absres = sqrt(fasp_blas_array_dotprod(n,w,r));
                    relres = absres/absres0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = r_norm;
                    normu  = MAX(SMALLREAL,fasp_blas_array_norm2(n,x->val));
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level > PRINT_NONE ) norms[iter] = relres;
            
            if ( relres <= tol ) {
                break;
            }
            else {
                // Need to restart
                fasp_array_cp(n, r, p[0]); i = 0;
            }
            
        } /* end of convergence check */
        
        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
    } /* end of main while loop */
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
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
    
    if ( iter >= MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

#if 0
static double estimate_spectral_radius(const double **A, int n, size_t k = 20)
{
    double *x = (double *)malloc(n* sizeof(double));
    double *y = (double *)malloc(n* sizeof(double));
	double *z = (double *)malloc(n* sizeof(double));
    double t;
	int i1,j1;
    
    // initialize x to random values in [0,1)
    //    cusp::copy(cusp::detail::random_reals<ValueType>(N), x);
    dvector px;
    px.row = n;
    px.val = x;
    
    fasp_dvec_rand(n, &px);
	
    for(size_t i = 0; i < k; i++)
    {
        //cusp::blas::scal(x, ValueType(1.0) / cusp::blas::nrmmax(x));
		t= 1.0/ fasp_blas_array_norminf(n, px);
		for(i1= 0; i1 <n; i1++) x[i1] *= t;
		
        //cusp::multiply(A, x, y);
		
		for(i1= 0; i1 <n; i1++) {
            t= 0.0
            for(j1= 0; j1 <n; j1++)  t +=  A[i1][j1] * x[j1];
            y[i1] = t;   }
        //		   x.swap(y);
		for(i1= 0; i1 <n; i1++) z[i1] = x[i1];
		for(i1= 0; i1 <n; i1++) x[i1] = y[i1];
        for(i1= 0; i1 <n; i1++) y[i1] = z[i1];
    }
    
    free(x);
	free(y);
	free(z);
	
    if (k == 0)
        return 0;
    else
        //return cusp::blas::nrm2(x) / cusp::blas::nrm2(y);
		return fasp_blas_array_norm2(n,x) / fasp_blas_array_norm2(n,y) ;
}

static double fasp_spectral_radius(dCSRmat *A,
                                   const SHORT restart)
{
    const INT n         = A->row;
    const INT MIN_ITER  = 0;
    
    // local variables
    INT      iter = 0;
    INT      restartplus1 = restart + 1;
    INT      i, j, k;
    
    REAL     epsmac = SMALLREAL;
    REAL     r_norm, den_norm;
    REAL     epsilon, gamma, t;
    
    REAL    *c = NULL, *s = NULL, *rs = NULL;
    REAL    *norms = NULL, *r = NULL, *w = NULL;
    REAL   **p = NULL, **hh = NULL;
    REAL    *work = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dstr_pgmres ...... [Start]\n");
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
    dvector p0;
    p0.row = n;
    p0.val = p[0];
    fasp_dvec_rand(n, &p0);
    // for (i=0;i<n ;i++) p[0][i] = random()
    
	r_norm = fasp_blas_array_norm2(n, p[0]);
    t = 1.0 / r_norm;
    for (j = 0; j < n; j ++) p[0][j] *= t;
    
	int  maxiter = std::min(n, restart) ;
	for(j = 0; j < maxiter; j++)
	{
        //		cusp::multiply(A, V[j], V[j + 1]);
        fasp_blas_dcsr_mxv(A, p[j], p[j+1]);
		
		for( i = 0; i <= j; i++)
		{
            //	H_(i,j) = cusp::blas::dot(V[i], V[j + 1]);
            hh[i][j] = fasp_blas_array_dotprod(n, p[i], p[j+1]);
			fasp_blas_array_axpy(n, -hh[i][j], p[i], p[ j+1 ]);
            //	cusp::blas::axpy(V[i], V[j + 1], -H_(i,j));
		}
        
		//H_(j+1,j) = cusp::blas::nrm2(V[j + 1]);
        hh[j+1][j] =  fasp_blas_array_norm2 (n, p[j+1]);
		if ( hh[j+1][j] < 1e-10) break;
		//cusp::blas::scal(V[j + 1], ValueType(1) / H_(j+1,j));
		t = 1.0/hh[j+1][j];
        for (int  k = 0; k < n; k ++) p[j+1][k] *= t;
	}
    
    //	H.resize(j,j);
    H   = (REAL **)fasp_mem_calloc(j, sizeof(REAL *));
	H[0]   = (REAL *)fasp_mem_calloc(j*j, sizeof(REAL));
	for (i = 1; i < j; i ++) H[i] = H[i-1] + j;
	
	
	for( size_t row = 0; row < j; row++ )
		for( size_t col = 0; col < j; col++ )
			H[row][col] = hh[row][col];
	
    double spectral_radius = estimate_spectral_radius( H, j, 20);
    
    
    /*-------------------------------------------
     * Clean up workspace
     *------------------------------------------*/
    fasp_mem_free(work);
    fasp_mem_free(p);
    fasp_mem_free(hh);
    
    fasp_mem_free(norms);
	fasp_mem_free(H[0]);
	fasp_mem_free(H);
    
    
    return spectral_radius;
}
#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
