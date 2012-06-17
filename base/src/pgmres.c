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
#include "its_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/*!
 * \fn INT fasp_solver_dcsr_pgmres (dCSRmat *A, dvector *b, dvector *x, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, const SHORT restart,
 *                                 const SHORT stop_type, const SHORT print_level)
 *
 * \brief Solve "Ax=b" using PGMRES (right preconditioned) iterative method
 *
 * \param A            Pointer to the coefficient matrix
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
 */ 
INT fasp_solver_dcsr_pgmres (dCSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             precond *pc, 
                             const REAL tol,
                             const INT MaxIt, 
                             const SHORT restart,
                             const SHORT stop_type, 
                             const SHORT print_level)
{
    const INT n         = A->row;  
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

    if ( print_level>PRINT_NONE ) printf("Calling GMRes solver ...\n");    

    /* allocate memory */
    work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
    
    p    = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    
    hh   = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
    
    if (print_level>PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;    
    
    for (i = 0; i < restartplus1; i ++) p[i] = s + restart + i*n;
    
    for (i = 0; i < restartplus1; i ++) hh[i] = p[restart] + n + i*restart;
    
    /* initialization */
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dcsr_aAxpy(-1.0, A, x->val, p[0]);
    
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
            fasp_array_cp(n, b->val, r);
            fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
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
            fasp_array_cp(n, b->val, r);
            fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
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

/**
 * \fn INT fasp_solver_bdcsr_pgmres (block_dCSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                   const REAL tol, const INT MaxIt, const SHORT restart,
 *                                   const SHORT stop_type, const SHORT print_level)
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param restart      Restarting steps
 * \param stop_type    Stopping criteria type -- DOES not support this parameter
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/24/2010
 *
 * Modified by Chensong Zhang on 05/01/2012
 */
INT fasp_solver_bdcsr_pgmres (block_dCSRmat *A, 
                              dvector *b, 
                              dvector *u, 
                              precond *pc, 
                              const REAL tol,
                              const INT MaxIt, 
                              const SHORT restart,
                              const SHORT stop_type, 
                              const SHORT print_level)
{
    const INT   nrow=b->row, nrow_1=nrow-1;
    const REAL  sol_inf_tol = SMALLREAL; // infinity norm tolrance

    // local variables
    INT         i, j, j1, index, iter=0, m=restart;
    REAL        beta, betai, tempe, tempb, tempu, hij, temp2;
    REAL        absres0=BIGREAL, absres, relres1, infnormu, factor;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_bdcsr_pgmres ...... [Start]\n");
#endif    

    if ((m<1)||(m>nrow)||(m>150)) m=20; // default restart level
    
    REAL *tmp=(REAL *)fasp_mem_calloc(m+1,sizeof(REAL));
    
    dvector *v=(dvector *)fasp_mem_calloc(m+1,sizeof(dvector));
    
    for (i=0;i<=m;++i) {
        v[i].row=nrow;
        v[i].val=(REAL*)fasp_mem_calloc(nrow,sizeof(REAL));
    }
    
    dvector y = fasp_dvec_create(m);
    
    REAL *r=(REAL*)fasp_mem_calloc(nrow,sizeof(REAL));
    
    REAL *z=(REAL *)fasp_mem_calloc(nrow,sizeof(REAL));
    
    REAL *w=(REAL *)fasp_mem_calloc(nrow,sizeof(REAL));    
    
    // generate the structure of H, i.e. H->IA, H->JA
    dCSRmat H=fasp_dcsr_create(m+1, m, m*(m+3)/2);
    
    H.IA[1]=m;
    for (i=2;i<=H.row;++i) H.IA[i]=H.IA[i-1]+m+2-i;
    
    for (i=0;i<H.row;++i) {
        if (i==0)
            index=0;
        else
            index=i-1;
    
        unsigned INT begin_row=H.IA[i], end_row1=H.IA[i+1];    
        for(j=begin_row;j<end_row1;++j) {
            H.JA[j]=index;
            index++;
        }
    }
    
    // compute norm for right hand side
    switch (stop_type) {
    case STOP_REL_PRECRES:
        if (pc != NULL) 
            pc->fct(b->val,z,pc->data); /* Preconditioning */
        else 
            fasp_array_cp(nrow,b->val,z); /* No preconditioner, B=I */
        tempb=sqrt(ABS(fasp_blas_array_dotprod(nrow,b->val,z)));
        break;
    case STOP_MOD_REL_RES:
        break;
    default: // STOP_REL_RES
        tempb=fasp_blas_array_norm2(nrow,b->val); // norm(b)
        break;
    }
    tempb=MAX(SMALLREAL,tempb);
    tempu=MAX(SMALLREAL,fasp_blas_array_norm2(nrow,u->val));
    
    // r = b-A*u
    fasp_array_cp(nrow,b->val,r); 
    fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
    tempe=fasp_blas_array_norm2(nrow,r);
    
    switch (stop_type) {
    case STOP_REL_PRECRES:
        if (pc == NULL)
            fasp_array_cp(nrow,r,z);
        else
            pc->fct(r,z,pc->data);
        temp2=sqrt(ABS(fasp_blas_array_dotprod(nrow,r,z)));
        relres1=temp2/tempb; 
        break;
    case STOP_MOD_REL_RES:
        relres1=tempe/tempu; 
        break;
    default: // STOP_REL_RES
        relres1=tempe/tempb; 
        break;
    }
    
    if (relres1<tol) { fasp_mem_free(r); goto FINISHED; }
    
    if  (iter < 0) goto FINISHED;
    
    while (iter++<MaxIt) {

        // z = B*r
        if (pc == NULL) {
            fasp_array_cp(nrow,r,z);  /* No preconditioner, B=I */
        }
        else {
            pc->fct(r,z,pc->data); /* Preconditioning */
        }
    
        beta=fasp_blas_array_norm2(nrow,z); betai=1.0/beta;
    
        // v_0 = z/beta
        for (i=0;i<=nrow_1;++i) v[0].val[i]=betai*z[i];
    
        for (j=0;j<m;++j) {
            // r = Av_j
            fasp_array_set(nrow, r, 0.0); 
            fasp_blas_bdcsr_aAxpy(1.0,A,v[j].val,r);
    
            // w = B*r
            if (pc == NULL) {
                fasp_array_cp(nrow,r,w);  /* No preconditioner, B=I */
            }
            else {
                pc->fct(r,w,pc->data); /* Preconditioning */
            }
    
            for (i=0;i<=j;++i) {
                if (i==0) index=0;
                else index=i-1;
                hij=fasp_blas_array_dotprod(nrow,w,v[i].val);
                H.val[H.IA[i]+j-index]=hij;
                fasp_blas_array_axpy(nrow, -hij, v[i].val, w);
            }
    
            j1=j+1;
            hij=fasp_blas_array_norm2(nrow,w);  // h_{j+1,j}=\|w\|_2
            H.val[H.IA[j1]]=hij;
    
            // v_{j+1}=w/h_{j+1,j}
            hij=1.0/hij;
            for (i=0;i<=nrow_1;++i) v[j1].val[i]=w[i]*hij;
        }
    
        fasp_aux_givens(beta, &H, &y, tmp);
    
        // u_m=u_0 + V_m*y_m
        for (i=0;i<m;++i) fasp_blas_array_axpy(nrow, y.val[i], v[i].val, u->val); // maybe we can check the residual for every iteration ?!
    
        // r = b-A*u
        fasp_array_cp(nrow,b->val,r); 
        fasp_blas_bdcsr_aAxpy(-1.0,A,u->val,r);
    
        // absolute and relative residuals
        absres=sqrt(fasp_blas_array_dotprod(nrow,r,r));
    
        tempu=sqrt(fasp_blas_dvec_dotprod(u,u));    
    
        switch (stop_type) {
        case STOP_REL_PRECRES:
            if (pc == NULL)
                fasp_array_cp(nrow,r,z);
            else
                pc->fct(r,z,pc->data);
            temp2=sqrt(ABS(fasp_blas_array_dotprod(nrow,r,z)));
            relres1=temp2/tempb; 
            break;
        case STOP_MOD_REL_RES:
            relres1=absres/tempu; 
            break;
        default: // STOP_REL_RES
            relres1=absres/tempb; 
            break;
        }
    
        // contraction factor
        factor=absres/absres0;
    
        // output iteration information if needed    
        print_itinfo(print_level,stop_type,iter,relres1,absres,factor);
    
        absres0=absres;
    
        // solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
        infnormu = fasp_blas_array_norminf(nrow, u->val); 
        if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            break;
        }
    
        if (relres1<tol) break;
    }
    
 FINISHED:
    if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,relres1);
    
    fasp_dvec_free(&y);
    fasp_dcsr_free(&H);
    fasp_mem_free(tmp);
    fasp_mem_free(z); 
    fasp_mem_free(w);
    fasp_mem_free(r);    
    for (i=0;i<=m;++i) fasp_mem_free(v[i].val);
    fasp_mem_free(v);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_bdcsr_pgmres ...... [Finish]\n");
#endif    
    
    return iter;
}

/*!
 * \fn INT fasp_solver_dbsr_pgmres (dBSRmat *A, dvector *b, dvector *x, precond *pc, 
 *                                  const REAL tol, const INT MaxIt, const SHORT restart,
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \param A            Pointer to the coefficient matrix
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
 * \date   2010/12/21
 *
 * Modified by Chensong Zhang on 05/01/2012
 */ 
INT fasp_solver_dbsr_pgmres (dBSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             precond *pc, 
                             const REAL tol,
                             const INT MaxIt, 
                             const SHORT restart,
                             const SHORT stop_type, 
                             const SHORT print_level)
{
    const INT n = A->ROW*A->nb;;  
    const INT min_iter = 0;
    
    // local variables
    INT       iter = 0;
    INT       restartplus1 = restart + 1;
    INT       i, j, k;
    
    REAL      epsmac = SMALLREAL; 
    REAL      r_norm, b_norm, den_norm;
    REAL      epsilon, gamma, t;   
    
    REAL     *c = NULL, *s = NULL, *rs = NULL; 
    REAL     *norms = NULL, *r = NULL, *w = NULL;
    REAL    **p = NULL, **hh = NULL;
    REAL     *work = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dbsr_pgmres ...... [Start]\n");
#endif    
    
    /* allocate memory */
    work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
    
    p  = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    
    hh = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
    
    if (print_level > PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;    
    for (i = 0; i < restartplus1; ++i) p[i] = s + restart + i*n;
    for (i = 0; i < restartplus1; ++i) hh[i] = p[restart] + n + i*restart;
    
    /* initialization */
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dbsr_aAxpy(-1.0, A, x->val, p[0]);
    
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
            fasp_array_cp(n, b->val, r);
            fasp_blas_dbsr_aAxpy(-1.0, A, x->val, r);
            r_norm = fasp_blas_array_norm2(n, r);
    
            if (r_norm <= epsilon) {
                break;
            }
            else {
                if (print_level >= PRINT_SOME) ITS_FACONV;
            }
        }
    
        t = 1.0 / r_norm;
        for (j = 0; j < n; ++j) p[0][j] *= t;
    
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while (i < restart && iter < MaxIt) {
            ++i;  ++iter;
    
            fasp_array_set(n, r, 0.0);
    
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);          
    
            fasp_blas_dbsr_mxv(A, r, p[i]);
    
            /* modified Gram_Schmidt */
            for (j = 0; j < i; ++j) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;    
            if (t != 0.0) {
                t = 1.0/t;
                for (j = 0; j < n; ++j) p[i][j] *= t;
            }
    
            for (j = 1; j < i; j++) {
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
        for (k = i-2; k >= 0; --k) {
            t = 0.0;
            for (j = k+1; j < i; ++j)  t -= hh[k][j]*rs[j];
    
            t += rs[k];
            rs[k] = t / hh[k][k];
        }
        fasp_array_cp(n, p[i-1], w);
        for (j = 0; j < n; ++j) w[j] *= rs[i-1];
        for (j = i-2; j >= 0; --j)  fasp_blas_array_axpy(n, rs[j], p[j], w);
        fasp_array_set(n, r, 0.0);
    
        /* apply the preconditioner */
        if (pc == NULL)
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
    
        fasp_blas_array_axpy(n, 1.0, r, x->val);
    
        if (r_norm  <= epsilon && iter >= min_iter) {
            fasp_array_cp(n, b->val, r);
            fasp_blas_dbsr_aAxpy(-1.0, A, x->val, r);
            r_norm = fasp_blas_array_norm2(n, r);
    
            if (r_norm  <= epsilon) {
                break;
            }
            else {
                if ( print_level >= PRINT_SOME ) ITS_FACONV;
                fasp_array_cp(n, r, p[0]); i = 0;
            }
        } /* end of convergence check */
    
        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j-1] = -s[j-1]*rs[j];
            rs[j] = c[j-1]*rs[j];
        }
    
        if (i) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
    
        for (j = i-1 ; j > 0; --j) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
    
        if (i) {
            fasp_blas_array_axpy(n, rs[0]-1.0, p[0], p[0]);
            fasp_blas_array_axpy(n, 1.0, p[i], p[0]);
        }        
    } /* end of iteration while loop */
    
    if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,r_norm);
    
    /*-------------------------------------------
     * Clean up workspace
     *------------------------------------------*/
    fasp_mem_free(work); 
    fasp_mem_free(p); 
    fasp_mem_free(hh);
    fasp_mem_free(norms);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dbsr_pgmres ...... [Finish]\n");
#endif    
    
    if (iter>=MaxIt) 
        return ERROR_SOLVER_MAXIT;
    else 
        return iter;
}

/*!
 * \fn INT fasp_solver_dstr_pgmres (dSTRmat *A, dvector *b, dvector *x, precond *pc, 
 *                                  const REAL tol, const INT MaxIt, const SHORT restart,
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \param A            Pointer to the coefficient matrix
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
 */ 
INT fasp_solver_dstr_pgmres (dSTRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             precond *pc, 
                             const REAL tol,
                             const INT MaxIt, 
                             const SHORT restart,
                             const SHORT stop_type, 
                             const SHORT print_level)
{
    const INT n = A->nc*A->ngrid;  
    const INT min_iter = 0;
    
    // local varialbes
    INT       iter = 0;
    INT       restartplus1 = restart + 1;
    INT       i, j, k;
    
    REAL      epsmac = SMALLREAL; 
    REAL      r_norm, b_norm, den_norm;
    REAL      epsilon, gamma, t;   
    
    REAL     *c = NULL, *s = NULL, *rs = NULL; 
    REAL     *norms = NULL, *r = NULL, *w = NULL;
    REAL    **p = NULL, **hh = NULL;
    REAL     *work = NULL;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dstr_pgmres ...... [Start]\n");
#endif    
    
    /* allocate memory */
    work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));
    
    p  = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));
    
    hh = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
    
    if (print_level > PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
    
    r = work; w = r + n; rs = w + n; c  = rs + restartplus1; s  = c + restart;    
    for (i = 0; i < restartplus1; ++i) p[i] = s + restart + i*n;
    for (i = 0; i < restartplus1; ++i) hh[i] = p[restart] + n + i*restart;
    
    /* initialization */
    fasp_array_cp(n, b->val, p[0]);
    fasp_blas_dstr_aAxpy(-1.0, A, x->val, p[0]);
    
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
            fasp_array_cp(n, b->val, r);
            fasp_blas_dstr_aAxpy(-1.0, A, x->val, r);
            r_norm = fasp_blas_array_norm2(n, r);
    
            if (r_norm <= epsilon) {
                break;
            }
            else {
                if (print_level >= PRINT_SOME) ITS_FACONV;
            }
        }
    
        t = 1.0 / r_norm;
        for (j = 0; j < n; ++j) p[0][j] *= t;
    
        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while (i < restart && iter < MaxIt) {
            ++i;  ++iter;
    
            fasp_array_set(n, r, 0.0);
    
            /* apply the preconditioner */
            if (pc == NULL)
                fasp_array_cp(n, p[i-1], r);
            else
                pc->fct(p[i-1], r, pc->data);          
    
            //fasp_blas_dbsr_mxv(A, r, p[i]);
            fasp_blas_dstr_aAxpy(1.0, A, r, p[i]); 
    
            /* modified Gram_Schmidt */
            for (j = 0; j < i; ++j) {
                hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
                fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
            }
            t = fasp_blas_array_norm2(n, p[i]);
            hh[i][i-1] = t;    
            if (t != 0.0) {
                t = 1.0/t;
                for (j = 0; j < n; ++j) p[i][j] *= t;
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
        for (k = i-2; k >= 0; k --)
            {
                t = 0.0;
                for (j = k+1; j < i; ++j)  t -= hh[k][j]*rs[j];
    
                t += rs[k];
                rs[k] = t / hh[k][k];
            }
        fasp_array_cp(n, p[i-1], w);
        for (j = 0; j < n; ++j) w[j] *= rs[i-1];
        for (j = i-2; j >= 0; j --)  fasp_blas_array_axpy(n, rs[j], p[j], w);
        fasp_array_set(n, r, 0.0);
    
        /* apply the preconditioner */
        if (pc == NULL)
            fasp_array_cp(n, w, r);
        else
            pc->fct(w, r, pc->data);
    
        fasp_blas_array_axpy(n, 1.0, r, x->val);
    
        if (r_norm  <= epsilon && iter >= min_iter) {
            fasp_array_cp(n, b->val, r);
            fasp_blas_dstr_aAxpy(-1.0, A, x->val, r);
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
    
    if (print_level>PRINT_NONE) ITS_FINAL(iter,MaxIt,r_norm);
    
    /*-------------------------------------------
     * Clean up workspace
     *------------------------------------------*/
    fasp_mem_free(work); 
    fasp_mem_free(p); 
    fasp_mem_free(hh);
    fasp_mem_free(norms);

#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dstr_pgmres ...... [Finish]\n");
#endif    
    
    if (iter>=MaxIt) 
        return ERROR_SOLVER_MAXIT;
    else 
        return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
