/*! \file spvgmres.c
 *
 *  \brief Krylov subspace methods -- Preconditioned variable-restart GMRes with 
 *         safe net
 *
 *  \note Refer to A.H. Baker, E.R. Jessup, and Tz.V. Kolev
 *        A Simple Strategy for Varying the Restart Parameter in GMRES(m)
 *        Journal of Computational and Applied Mathematics, 230 (2009)
 *        pp. 751-761. UCRL-JRNL-235266.
 *
 *  \note See pvgmres.c a version without safe net
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
 * \fn INT fasp_solver_dcsr_spvgmres (dCSRmat *A, dvector *b, dvector *x, precond *pc,
 *                                    const REAL tol, const INT MaxIt, SHORT restart,
 *                                    const SHORT stop_type, const SHORT print_level)
 *
 * \brief Solve "Ax=b" using PGMRES(right preconditioned) iterative method in which the restart
 *        parameter can be adaptively modified during the iteration.
 *
 * \param A            Pointer to dCSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to the structure of precondition (precond)
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Chensong Zhang
 * \date   04/06/2013
 * Modified by Chunsheng Feng on 07/22/2013: Add adapt memory allocate
 */
INT fasp_solver_dcsr_spvgmres (dCSRmat *A,
                               dvector *b,
                               dvector *x,
                               precond *pc,
                               const REAL tol,
                               const INT MaxIt,
                               SHORT restart,
                               const SHORT stop_type,
                               const SHORT print_level)
{
    const INT   n          = b->row;
    const INT   MIN_ITER   = 0;
    const REAL  maxdiff    = tol*STAG_RATIO; // staganation tolerance
    const REAL  epsmac     = SMALLREAL;
    
    //--------------------------------------------//
    //   Newly added parameters to monitor when   //
    //   to change the restart parameter          //
    //--------------------------------------------//
    const REAL cr_max      = 0.99;    // = cos(8^o)  (experimental)
    const REAL cr_min      = 0.174;   // = cos(80^o) (experimental)
    
    // local variables
    INT    iter            = 0;
    INT    restartplus1    = restart + 1;
    INT    i, j, k;
    
    REAL   r_norm, r_normb, gamma, t;
    REAL   absres0, absres, relres, normu;
    
    REAL   cr          = 1.0;     // convergence rate
    REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
    INT    d           = 3;       // reduction for the restart parameter
    INT    restart_max = restart; // upper bound for restart in each restart cycle
    INT    restart_min = 3;       // lower bound for restart in each restart cycle (should be small)
    INT    Restart;               // the real restart in some fixed restarted cycle
    
    INT    iter_best = 0;         // initial best known iteration
    REAL   absres_best = BIGREAL; // initial best known residual
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL  *c = NULL, *s = NULL, *rs = NULL;
    REAL  *norms = NULL, *r = NULL, *w = NULL;
    REAL  *work = NULL, *x_best = NULL;
    REAL  **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));

    while ( (work == NULL) && (restart > 5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s: %d !\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c = rs + restartplus1;
    x_best = c + restart; s = x_best + n;
    
    for ( i = 0; i < restartplus1; i++ ) p[i] = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dcsr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
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
        
        rs[0] = r_norm_old = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        //-----------------------------------//
        //   adjust the restart parameter    //
        //-----------------------------------//
        if ( cr > cr_max || iter == 0 ) {
            Restart = restart_max;
        }
        else if ( cr < cr_min ) {
            Restart = Restart;
        }
        else {
            if ( Restart - d > restart_min ) {
                Restart -= d;
            }
            else {
                Restart = restart_max;
            }
        }
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < Restart && iter < MaxIt ) {
            
            i++;  iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_dcsr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for (j = 0; j < i; j ++) {
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
            rs[i-1] = c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if (print_level > PRINT_NONE) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
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
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- )  fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // safe net check: save the best-so-far solution
        if ( fasp_dvec_isnan(x) ) {
            // If the solution is NAN, restrore the best solution
            absres = BIGREAL;
            goto RESTORE_BESTSOL;
        }
        
        if ( absres < absres_best - maxdiff) {
            absres_best = absres;
            iter_best   = iter;
            fasp_array_cp(n,x->val,x_best);
        }
        
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
        for ( j = i; j > 0; j-- ) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
        //-----------------------------------//
        //   compute the convergence rate    //
        //-----------------------------------//
        cr = r_norm / r_norm_old;
        
    } /* end of iteration while loop */
    
RESTORE_BESTSOL: // restore the best-so-far solution if necessary
    if ( iter != iter_best ) {
        
        // compute best residual
        fasp_array_cp(n,b->val,r);
        fasp_blas_dcsr_aAxpy(-1.0,A,x_best,r);
        
        switch ( stop_type ) {
            case STOP_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
            case STOP_REL_PRECRES:
                // z = B(r)
                if ( pc != NULL )
                    pc->fct(r,w,pc->data); /* Apply preconditioner */
                else
                    fasp_array_cp(n,r,w); /* No preconditioner */
                absres_best = sqrt(ABS(fasp_blas_array_dotprod(n,w,r)));
                break;
            case STOP_MOD_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
        }
        
        if ( absres > absres_best + maxdiff ) {
            if ( print_level > PRINT_NONE ) ITS_RESTORE(iter);
            fasp_array_cp(n,x_best,x->val);
        }
    }
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    /*-------------------------------------------
     * Free some stuff
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

/**
 * \fn INT fasp_solver_bdcsr_spvgmres (block_dCSRmat *A, dvector *b, dvector *x, precond *pc,
 *                                     const REAL tol, const INT MaxIt, SHORT restart,
 *                                     const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned GMRES method for solving Au=b
 *
 * \param A            Pointer to block_dCSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to the structure of precondition (precond)
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Chensong Zhang
 * \date   04/06/2013
 */
INT fasp_solver_bdcsr_spvgmres (block_dCSRmat *A,
                                dvector *b,
                                dvector *x,
                                precond *pc,
                                const REAL tol,
                                const INT MaxIt,
                                SHORT restart,
                                const SHORT stop_type,
                                const SHORT print_level)
{
    const INT   n          = b->row;
    const INT   MIN_ITER   = 0;
    const REAL  maxdiff    = tol*STAG_RATIO; // staganation tolerance
    const REAL  epsmac     = SMALLREAL;
    
    //--------------------------------------------//
    //   Newly added parameters to monitor when   //
    //   to change the restart parameter          //
    //--------------------------------------------//
    const REAL cr_max      = 0.99;    // = cos(8^o)  (experimental)
    const REAL cr_min      = 0.174;   // = cos(80^o) (experimental)
    
    // local variables
    INT    iter            = 0;
    INT    restartplus1    = restart + 1;
    INT    i, j, k;
    
    REAL   r_norm, r_normb, gamma, t;
    REAL   absres0, absres, relres, normu;
    
    REAL   cr          = 1.0;     // convergence rate
    REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
    INT    d           = 3;       // reduction for the restart parameter
    INT    restart_max = restart; // upper bound for restart in each restart cycle
    INT    restart_min = 3;       // lower bound for restart in each restart cycle (should be small)
    INT    Restart;               // the real restart in some fixed restarted cycle
    
    INT    iter_best = 0;         // initial best known iteration
    REAL   absres_best = BIGREAL; // initial best known residual
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL  *c = NULL, *s = NULL, *rs = NULL;
    REAL  *norms = NULL, *r = NULL, *w = NULL;
    REAL  *work = NULL, *x_best = NULL;
    REAL  **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));

    while ( (work == NULL) && (restart > restart_min+5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
        restart_max = restart;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s: %d !\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c = rs + restartplus1;
    x_best = c + restart; s = x_best + n;
    
    for ( i = 0; i < restartplus1; i++ ) p[i] = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_bdcsr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
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
        
        rs[0] = r_norm_old = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        //-----------------------------------//
        //   adjust the restart parameter    //
        //-----------------------------------//
        if ( cr > cr_max || iter == 0 ) {
            Restart = restart_max;
        }
        else if ( cr < cr_min ) {
            Restart = Restart;
        }
        else {
            if ( Restart - d > restart_min ) {
                Restart -= d;
            }
            else {
                Restart = restart_max;
            }
        }
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < Restart && iter < MaxIt ) {
            
            i++;  iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_bdcsr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for (j = 0; j < i; j ++) {
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
            rs[i-1] = c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if (print_level > PRINT_NONE) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
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
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- )  fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // safe net check: save the best-so-far solution
        if ( fasp_dvec_isnan(x) ) {
            // If the solution is NAN, restrore the best solution
            absres = BIGREAL;
            goto RESTORE_BESTSOL;
        }
        
        if ( absres < absres_best - maxdiff) {
            absres_best = absres;
            iter_best   = iter;
            fasp_array_cp(n,x->val,x_best);
        }
        
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
        for ( j = i; j > 0; j-- ) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
        //-----------------------------------//
        //   compute the convergence rate    //
        //-----------------------------------//
        cr = r_norm / r_norm_old;
        
    } /* end of iteration while loop */
    
RESTORE_BESTSOL: // restore the best-so-far solution if necessary
    if ( iter != iter_best ) {
        
        // compute best residual
        fasp_array_cp(n,b->val,r);
        fasp_blas_bdcsr_aAxpy(-1.0,A,x_best,r);
        
        switch ( stop_type ) {
            case STOP_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
            case STOP_REL_PRECRES:
                // z = B(r)
                if ( pc != NULL )
                    pc->fct(r,w,pc->data); /* Apply preconditioner */
                else
                    fasp_array_cp(n,r,w); /* No preconditioner */
                absres_best = sqrt(ABS(fasp_blas_array_dotprod(n,w,r)));
                break;
            case STOP_MOD_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
        }
        
        if ( absres > absres_best + maxdiff ) {
            if ( print_level > PRINT_NONE ) ITS_RESTORE(iter);
            fasp_array_cp(n,x_best,x->val);
        }
    }
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    /*-------------------------------------------
     * Free some stuff
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

/*!
 * \fn INT fasp_solver_dbsr_spvgmres (dBSRmat *A, dvector *b, dvector *x, precond *pc,
 *                                    const REAL tol, const INT MaxIt, SHORT restart,
 *                                    const SHORT stop_type, const SHORT print_level)
 *
 * \brief Solve "Ax=b" using PGMRES(right preconditioned) iterative method in which the restart
 *        parameter can be adaptively modified during the iteration.
 *
 * \param A            Pointer to dBSRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to the structure of precondition (precond)
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Chensong Zhang
 * \date   04/06/2013
 */
INT fasp_solver_dbsr_spvgmres (dBSRmat *A,
                               dvector *b,
                               dvector *x,
                               precond *pc,
                               const REAL tol,
                               const INT MaxIt,
                               SHORT restart,
                               const SHORT stop_type,
                               const SHORT print_level)
{
    const INT   n          = b->row;
    const INT   MIN_ITER   = 0;
    const REAL  maxdiff    = tol*STAG_RATIO; // staganation tolerance
    const REAL  epsmac     = SMALLREAL;
    
    //--------------------------------------------//
    //   Newly added parameters to monitor when   //
    //   to change the restart parameter          //
    //--------------------------------------------//
    const REAL cr_max      = 0.99;    // = cos(8^o)  (experimental)
    const REAL cr_min      = 0.174;   // = cos(80^o) (experimental)
    
    // local variables
    INT    iter            = 0;
    INT    restartplus1    = restart + 1;
    INT    i, j, k;
    
    REAL   r_norm, r_normb, gamma, t;
    REAL   absres0, absres, relres, normu;
    
    REAL   cr          = 1.0;     // convergence rate
    REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
    INT    d           = 3;       // reduction for the restart parameter
    INT    restart_max = restart; // upper bound for restart in each restart cycle
    INT    restart_min = 3;       // lower bound for restart in each restart cycle (should be small)
    INT    Restart;               // the real restart in some fixed restarted cycle
    
    INT    iter_best = 0;         // initial best known iteration
    REAL   absres_best = BIGREAL; // initial best known residual
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL  *c = NULL, *s = NULL, *rs = NULL;
    REAL  *norms = NULL, *r = NULL, *w = NULL;
    REAL  *work = NULL, *x_best = NULL;
    REAL  **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));

    while ( (work == NULL) && (restart > restart_min+5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
        restart_max = restart;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s: %d !\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c = rs + restartplus1;
    x_best = c + restart; s = x_best + n;
    
    for ( i = 0; i < restartplus1; i++ ) p[i] = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dbsr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
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
        
        rs[0] = r_norm_old = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        //-----------------------------------//
        //   adjust the restart parameter    //
        //-----------------------------------//
        if ( cr > cr_max || iter == 0 ) {
            Restart = restart_max;
        }
        else if ( cr < cr_min ) {
            Restart = Restart;
        }
        else {
            if ( Restart - d > restart_min ) {
                Restart -= d;
            }
            else {
                Restart = restart_max;
            }
        }
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < Restart && iter < MaxIt ) {
            
            i++;  iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_dbsr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for (j = 0; j < i; j ++) {
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
            rs[i-1] = c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if (print_level > PRINT_NONE) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
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
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- )  fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // safe net check: save the best-so-far solution
        if ( fasp_dvec_isnan(x) ) {
            // If the solution is NAN, restrore the best solution
            absres = BIGREAL;
            goto RESTORE_BESTSOL;
        }
        
        if ( absres < absres_best - maxdiff) {
            absres_best = absres;
            iter_best   = iter;
            fasp_array_cp(n,x->val,x_best);
        }
        
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
        for ( j = i; j > 0; j-- ) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
        //-----------------------------------//
        //   compute the convergence rate    //
        //-----------------------------------//
        cr = r_norm / r_norm_old;
        
    } /* end of iteration while loop */
    
RESTORE_BESTSOL: // restore the best-so-far solution if necessary
    if ( iter != iter_best ) {
        
        // compute best residual
        fasp_array_cp(n,b->val,r);
        fasp_blas_dbsr_aAxpy(-1.0,A,x_best,r);
        
        switch ( stop_type ) {
            case STOP_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
            case STOP_REL_PRECRES:
                // z = B(r)
                if ( pc != NULL )
                    pc->fct(r,w,pc->data); /* Apply preconditioner */
                else
                    fasp_array_cp(n,r,w); /* No preconditioner */
                absres_best = sqrt(ABS(fasp_blas_array_dotprod(n,w,r)));
                break;
            case STOP_MOD_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
        }
        
        if ( absres > absres_best + maxdiff ) {
            if ( print_level > PRINT_NONE ) ITS_RESTORE(iter);
            fasp_array_cp(n,x_best,x->val);
        }
    }
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    /*-------------------------------------------
     * Free some stuff
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

/*!
 * \fn INT fasp_solver_dstr_spvgmres (dSTRmat *A, dvector *b, dvector *x, precond *pc,
 *                                    const REAL tol, const INT MaxIt, SHORT restart,
 *                                    const SHORT stop_type, const SHORT print_level)
 *
 * \brief Solve "Ax=b" using PGMRES(right preconditioned) iterative method in which the restart
 *        parameter can be adaptively modified during the iteration.
 *
 * \param A            Pointer to dSTRmat: the coefficient matrix
 * \param b            Pointer to dvector: the right hand side
 * \param x            Pointer to dvector: the unknowns
 * \param pc           Pointer to the structure of precondition (precond)
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Chensong Zhang
 * \date   04/06/2013
 */
INT fasp_solver_dstr_spvgmres (dSTRmat *A,
                               dvector *b,
                               dvector *x,
                               precond *pc,
                               const REAL tol,
                               const INT MaxIt,
                               SHORT restart,
                               const SHORT stop_type,
                               const SHORT print_level)
{
    const INT   n          = b->row;
    const INT   MIN_ITER   = 0;
    const REAL  maxdiff    = tol*STAG_RATIO; // staganation tolerance
    const REAL  epsmac     = SMALLREAL;
    
    //--------------------------------------------//
    //   Newly added parameters to monitor when   //
    //   to change the restart parameter          //
    //--------------------------------------------//
    const REAL cr_max      = 0.99;    // = cos(8^o)  (experimental)
    const REAL cr_min      = 0.174;   // = cos(80^o) (experimental)
    
    // local variables
    INT    iter            = 0;
    INT    restartplus1    = restart + 1;
    INT    i, j, k;
    
    REAL   r_norm, r_normb, gamma, t;
    REAL   absres0, absres, relres, normu;
    
    REAL   cr          = 1.0;     // convergence rate
    REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
    INT    d           = 3;       // reduction for the restart parameter
    INT    restart_max = restart; // upper bound for restart in each restart cycle
    INT    restart_min = 3;       // lower bound for restart in each restart cycle (should be small)
    INT    Restart;               // the real restart in some fixed restarted cycle
    
    INT    iter_best = 0;         // initial best known iteration
    REAL   absres_best = BIGREAL; // initial best known residual
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL  *c = NULL, *s = NULL, *rs = NULL;
    REAL  *norms = NULL, *r = NULL, *w = NULL;
    REAL  *work = NULL, *x_best = NULL;
    REAL  **p = NULL, **hh = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    /* allocate memory and setup temp work space */
    work  = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));

    while ( (work == NULL) && (restart > restart_min+5 ) ) {
        restart = restart - 5 ;
        work = (REAL *) fasp_mem_calloc((restart+4)*(restart+n)+1, sizeof(REAL));
        printf("### WARNING: GMRES restart number becomes %d!\n", restart );
        restartplus1 = restart + 1;
        restart_max = restart;
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s: %d !\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }

    p     = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    hh    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    r = work; w = r + n; rs = w + n; c = rs + restartplus1;
    x_best = c + restart; s = x_best + n;
    
    for ( i = 0; i < restartplus1; i++ ) p[i] = s + restart + i*n;
    
    for ( i = 0; i < restartplus1; i++ ) hh[i] = p[restart] + n + i*restart;
    
    // r = b-A*x
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dstr_aAxpy(-1.0, A, x->val, p[0]);
    
    r_norm = fasp_blas_array_norm2(n, p[0]);
    
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
        
        rs[0] = r_norm_old = r_norm;
        
        t = 1.0 / r_norm;
        
        fasp_blas_array_ax(n, t, p[0]);
        
        //-----------------------------------//
        //   adjust the restart parameter    //
        //-----------------------------------//
        if ( cr > cr_max || iter == 0 ) {
            Restart = restart_max;
        }
        else if ( cr < cr_min ) {
            Restart = Restart;
        }
        else {
            if ( Restart - d > restart_min ) {
                Restart -= d;
            }
            else {
                Restart = restart_max;
            }
        }
        
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while ( i < Restart && iter < MaxIt ) {
            
            i++;  iter++;
            
            fasp_array_set(n, r, 0.0);
            
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);
            
            fasp_blas_dstr_mxv(A, r, p[i]);
            
            /* modified Gram_Schmidt */
            for (j = 0; j < i; j ++) {
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
            rs[i-1] = c[i-1]*rs[i-1];
            hh[i-1][i-1] = s[i-1]*hh[i][i-1] + c[i-1]*hh[i-1][i-1];
            
            absres = r_norm = fabs(rs[i]);
            
            relres = absres/absres0;
            
            if (print_level > PRINT_NONE) norms[iter] = relres;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, iter, relres, absres,
                         norms[iter]/norms[iter-1]);
            
            // should we exit the restart cycle
            if ( relres <= tol && iter >= MIN_ITER ) break;
            
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
        
        fasp_blas_array_ax(n, rs[i-1], w);
        
        for ( j = i-2; j >= 0; j-- )  fasp_blas_array_axpy(n, rs[j], p[j], w);
        
        fasp_array_set(n, r, 0.0);
        
        /* apply the preconditioner */
        if ( pc == NULL )
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
        
        fasp_blas_array_axpy(n, 1.0, r, x->val);
        
        // safe net check: save the best-so-far solution
        if ( fasp_dvec_isnan(x) ) {
            // If the solution is NAN, restrore the best solution
            absres = BIGREAL;
            goto RESTORE_BESTSOL;
        }
        
        if ( absres < absres_best - maxdiff) {
            absres_best = absres;
            iter_best   = iter;
            fasp_array_cp(n,x->val,x_best);
        }
        
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
        for ( j = i; j > 0; j-- ) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j]   = c[j-1]*rs[j];
        }
        
        if ( i ) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
        
        for ( j = i-1 ; j > 0; j-- ) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
        
        if ( i ) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
        
        //-----------------------------------//
        //   compute the convergence rate    //
        //-----------------------------------//
        cr = r_norm / r_norm_old;
        
    } /* end of iteration while loop */
    
RESTORE_BESTSOL: // restore the best-so-far solution if necessary
    if ( iter != iter_best ) {
        
        // compute best residual
        fasp_array_cp(n,b->val,r);
        fasp_blas_dstr_aAxpy(-1.0,A,x_best,r);
        
        switch ( stop_type ) {
            case STOP_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
            case STOP_REL_PRECRES:
                // z = B(r)
                if ( pc != NULL )
                    pc->fct(r,w,pc->data); /* Apply preconditioner */
                else
                    fasp_array_cp(n,r,w); /* No preconditioner */
                absres_best = sqrt(ABS(fasp_blas_array_dotprod(n,w,r)));
                break;
            case STOP_MOD_REL_RES:
                absres_best = fasp_blas_array_norm2(n,r);
                break;
        }
        
        if ( absres > absres_best + maxdiff ) {
            if ( print_level > PRINT_NONE ) ITS_RESTORE(iter);
            fasp_array_cp(n,x_best,x->val);
        }
    }
    
FINISHED:
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    /*-------------------------------------------
     * Free some stuff
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
