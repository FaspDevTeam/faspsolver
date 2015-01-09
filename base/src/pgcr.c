/*! \file pgcr.c
 *  \brief Krylov subspace methods -- Preconditioned GCR.
 */
#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

static void dense_aAxpby (INT, INT, REAL *, REAL, REAL *, REAL, REAL *);

/**
 * \fn int fasp_solver_dcsr_pgcr1 (dCSRmat *A,  dvector *b, dvector *u, \
 *                                  precond *pc, const REAL tol, const INT MaxIt, \
 *                                  const INT restart, const INT stop_type, \
 *                                  const INT print_level)
 *
 * \brief A preconditioned GCR method for solving Au=b
 *
 * \param *A     Pointer to the coefficient matrix
 * \param *b     Pointer to the dvector of right hand side
 * \param *u     Pointer to the dvector of dofs
 * \param MaxIt Maximal number of iterations
 * \param tol   Tolerance for stopage
 * \param *pre  Pointer to the structure of precondition (precond)
 * \param print_level How much information to print out
 *
 * \return the number of iterations
 *
 * \author Lu Wang
 * \date   11/02/2014
 */
INT fasp_solver_dcsr_pgcr1 (dCSRmat *A,
                            dvector *b,
                            dvector *x,
                            precond *pc,
                            const REAL tol,
                            const INT MaxIt,
                            const INT restart,
                            const INT stop_type,
                            const INT print_level)
{
    INT i, j;
    INT iter = 0;
    INT m = restart;
    REAL alpha, beta, gamma, tempe, tempb, tempu,temp2;
    REAL absres0 = BIGREAL, absres, relres1, infnormu, factor;
    
    const INT nrow = b->row;
    const REAL sol_inf_tol = 1e-16;
    
    // default restart number
    if ((m < 1)||(m > nrow)||(m >150)) m=10;
    
    dvector *v = (dvector *)fasp_mem_calloc(m,sizeof(dvector));
    
    for (i=0; i<=m-1; ++i) {
        v[i].row = nrow;
        v[i].val = (REAL*)fasp_mem_calloc(nrow,sizeof(REAL));
    }
    
    dvector *s = (dvector *)fasp_mem_calloc(m,sizeof(dvector));
    
    for (i=0; i<=m-1;++i) {
        s[i].row = nrow;
        s[i].val = (REAL*)fasp_mem_calloc(nrow,sizeof(REAL));
    }
    
    REAL *r = (REAL*)fasp_mem_calloc(nrow,sizeof(REAL));
    
    // compute norm for right hand side
    switch (stop_type) {
        case STOP_REL_RES:
            tempb = fasp_blas_array_norm2(nrow,b->val);
            break;
            
        case STOP_REL_PRECRES:
            if (pc != NULL)
                pc->fct(b->val,s[0].val,pc->data);
            else
                fasp_array_cp(nrow,b->val,s[0].val);
            tempb = sqrt(ABS(fasp_blas_array_dotprod(nrow,b->val,s[0].val)));
            break;
            
        case STOP_MOD_REL_RES:
            break;
            
        default:
            printf("Error: Unknown stopping criteria!\n");
            iter = ERROR_INPUT_PAR;
            goto FINISHED;
    }
    
    tempb = MAX(SMALLREAL,tempb);
    tempu = MAX(SMALLREAL,fasp_blas_array_norm2(nrow,x->val));
    
    // r = b-A*u
    fasp_array_cp(nrow, b->val, r);
    fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
    tempe = fasp_blas_array_norm2(nrow, r);
    tempb = MAX(SMALLREAL, tempe);
    
    switch (stop_type) {
        case STOP_REL_RES:
            relres1 = tempe/tempb;
            break;
            
        case STOP_REL_PRECRES:
            if (pc == NULL)
                fasp_array_cp(nrow, r, s[0].val);
            else
                pc->fct(r, s[0].val, pc->data);
            temp2 = sqrt(ABS(fasp_blas_array_dotprod(nrow, r, s[0].val)));
            relres1 = temp2/tempb;
            break;
            
        case STOP_MOD_REL_RES:
            relres1 = tempe/tempu;
            break;
    }
    
    if (relres1<tol) { fasp_mem_free(r); goto FINISHED; }
    
    if (iter < 0) goto FINISHED;
    
    while (iter++ < MaxIt)
    {
        for (j=0; j<m; ++j)
        {
            if (pc == NULL) {
                fasp_array_cp(nrow, r, s[j].val);
            }
            else {
                pc->fct(r, s[j].val, pc->data);
            }
            
            fasp_blas_dcsr_aAxpy(1.0, A, s[j].val, v[j].val);
            
            for (i=0; i<j; ++i)
            {
                alpha = fasp_blas_array_dotprod(nrow, v[j].val, v[i].val);
                fasp_blas_array_axpy(nrow, -alpha, v[i].val, v[j].val);
                fasp_blas_array_axpy(nrow, -alpha, s[i].val, s[j].val);
            }
            
            beta = fasp_blas_array_norm2(nrow, v[j].val);
            fasp_blas_array_ax(nrow, 1.0/beta, v[j].val);
            fasp_blas_array_ax(nrow, 1.0/beta, s[j].val);
            
            gamma = fasp_blas_array_dotprod(nrow, v[j].val, r);
            fasp_blas_array_axpy(nrow, gamma, s[j].val, x->val);
            //fasp_blas_array_axpy(nrow, -gamma, v[j].val, r);
            
            // r = b-A*u
            fasp_array_cp(nrow, b->val, r);
            fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
            
            // absolute and relative residuals
            absres = sqrt(fasp_blas_array_dotprod(nrow, r, r));
            tempu = sqrt(fasp_blas_dvec_dotprod(x, x));
            
            switch (stop_type) {
                case STOP_REL_RES:
                    relres1 = absres/tempb;
                    break;
                    
                case STOP_REL_PRECRES:
                    if (pc == NULL)
                        fasp_array_cp(nrow, r, s[j].val);
                    else
                        pc->fct(r, s[j].val, pc->data);
                    temp2 = sqrt(ABS(fasp_blas_array_dotprod(nrow, r, s[j].val)));
                    relres1 = temp2/tempb;
                    break;
                    
                case STOP_MOD_REL_RES:
                    relres1 = absres/tempu;
                    break;
            }
            
            // contraction factor
            factor = absres/absres0;
            
            // output iteration information if needed
            print_itinfo(print_level, stop_type, (iter-1)*m+j+1, relres1, absres, factor);
            
            absres0 = absres;
            
            // solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
            infnormu = fasp_blas_array_norminf(nrow, x->val);
            
            if (infnormu <= sol_inf_tol)
            {
                print_message(print_level, "GMRes stops: infinity norm of the solution is too small!\n");
                iter = ERROR_SOLVER_SOLSTAG;
                goto FINISHED;
            }
            
            if (relres1<tol) goto FINISHED;
        }
    }
    
FINISHED:
    
    if (print_level > 0) {
        if (iter > MaxIt){
            printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres1);
            iter = ERROR_SOLVER_MAXIT;
        }
        else
            printf("Number of iterations = %d with relative residual %e.\n", (iter-1)*m+j+1, relres1);
    }
    
    fasp_mem_free(r);
    for (i=0; i<=m-1; ++i) fasp_mem_free(v[i].val);
    fasp_mem_free(v);
    for (i=0; i<=m-1; ++i) fasp_mem_free(s[i].val);
    fasp_mem_free(s);
    
    return iter;
}

/**
 * \fn int fasp_solver_dcsr_pgcr (dCSRmat *A,  dvector *b, dvector *x, \
 *                                precond *pc, const REAL tol, const INT MaxIt, \
 *                                const INT restart, const INT stop_type, \
 *                                const INT print_level)
 *
 * \brief A preconditioned GCR method for solving Au=b
 *
 * \param *A     Pointer to the coefficient matrix
 * \param *b     Pointer to the dvector of right hand side
 * \param *u     Pointer to the dvector of dofs
 * \param MaxIt  Maximal number of iterations
 * \param tol    Tolerance for stopage
 * \param *pre   Pointer to the structure of precondition (precond)
 * \param print_level How much information to print out
 *
 * \return the number of iterations
 *
 * \author Zheng Li
 * \date   12/23/2014
 */
INT fasp_solver_dcsr_pgcr (dCSRmat *A,
                           dvector *b,
                           dvector *x,
                           precond *pc,
                           const REAL tol,
                           const INT MaxIt,
                           const SHORT restart,
                           const SHORT stop_type,
                           const SHORT print_level)
{
    const INT   n         = b->row;
    
    // local variables
    INT      iter = 0, rst = -1;
    INT      i, j, k;
    
    REAL     gamma, alpha, beta, checktol;
    REAL     absres0 = BIGREAL, absres = BIGREAL;
    REAL     relres  = BIGREAL;
    
    // allocate temp memory (need about (restart+4)*n REAL numbers)
    REAL    *c = NULL, *z = NULL, *alp = NULL, *tmpx = NULL;
    REAL    *norms = NULL, *r = NULL, *work = NULL;
    REAL    **h = NULL;
    
    INT      Restart = MIN(restart, MaxIt);
    LONG     worksize = n+2*Restart*n+Restart+Restart;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    work = (REAL *) fasp_mem_calloc(worksize, sizeof(REAL));
    
    /* check whether memory is enough for GMRES */
    while ( (work == NULL) && (Restart > 5) ) {
        Restart = Restart - 5;
        worksize = n+2*Restart*n+Restart+Restart;
        work = (REAL *) fasp_mem_calloc(worksize, sizeof(REAL));
    }
    
    if ( work == NULL ) {
        printf("### ERROR: No enough memory for GMRES %s : %s : %d!\n",
               __FILE__, __FUNCTION__, __LINE__ );
        exit(ERROR_ALLOC_MEM);
    }
    
    if ( print_level > PRINT_MIN && Restart < restart ) {
        printf("### WARNING: GMRES restart number set to %d!\n", Restart);
    }
    
    r = work, z = r+n, c = z + Restart*n, alp = c + Restart*n, tmpx = alp + Restart;
    
    h = (REAL **)fasp_mem_calloc(Restart, sizeof(REAL *));
    for (i = 0; i < Restart; i++) h[i] = (REAL*)fasp_mem_calloc(Restart, sizeof(REAL));
    
    norms = (REAL *) fasp_mem_calloc(MaxIt+1, sizeof(REAL));
    
    // r = b-A*x
    fasp_array_cp(n, b->val, r);
    fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
    
    absres = fasp_blas_array_dotprod(n, r, r);
    
    absres0 = MAX(SMALLREAL,absres);
    
    relres  = absres/absres0;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,0,relres,sqrt(absres0),0.0);
    
    // store initial residual
    norms[0] = relres;
    
    checktol = MAX(tol*tol*absres0, absres*1.0e-4);
    
    while ( iter < MaxIt && sqrt(relres) > tol ) {
        
        i = -1; rst ++;
        
        while ( i < Restart-1 && iter < MaxIt ) {
            
            i++; iter++;
            
            // z = B^-1r
            if ( pc == NULL )
                fasp_array_cp(n, r, &z[i*n]);
            else
                pc->fct(r, &z[i*n], pc->data);
            
            // c = Az
            fasp_blas_dcsr_mxv(A, &z[i*n], &c[i*n]);
            
            /* Modified Gram_Schmidt orthogonalization */
            for ( j = 0; j < i; j++ ) {
                gamma = fasp_blas_array_dotprod(n, &c[j*n], &c[i*n]);
                h[i][j] = gamma/h[j][j];
                fasp_blas_array_axpy(n, -h[i][j], &c[j*n], &c[i*n]);
            }
            // gamma = (c,c)
            gamma = fasp_blas_array_dotprod(n, &c[i*n], &c[i*n]);
            
            h[i][i] = gamma;
            
            // alpha = (c, r)
            alpha = fasp_blas_array_dotprod(n, &c[i*n], r);
            
            beta = alpha/gamma;
            
            alp[i] = beta;
            
            // r = r - beta*c
            fasp_blas_array_axpy(n, -beta, &c[i*n], r);
            
            // equivalent to ||r||_2
            absres = absres - alpha*alpha/gamma;
            
            if (absres < checktol) {
                absres = fasp_blas_array_dotprod(n, r, r);
                checktol = MAX(tol*tol*absres0, absres*1.0e-4);
            }
            
            relres = absres / absres0;
            
            norms[iter] = relres;
            
            print_itinfo(print_level, stop_type, iter, sqrt(relres), sqrt(absres),
                         sqrt(norms[iter]/norms[iter-1]));
            
            if (sqrt(relres) < tol)  break;
        }
        
        for ( k = i; k >=0; k-- ) {
            tmpx[k] = alp[k];
            for (j=0; j<k; ++j) {
                alp[j] -= h[k][j]*tmpx[k];
            }
        }
        
        if (rst==0) dense_aAxpby(n, i+1, z, 1.0, tmpx, 0.0, x->val);
        else dense_aAxpby(n, i+1, z, 1.0, tmpx, 1.0, x->val);
        
    }
    
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,sqrt(relres));
    
    // clean up memory
    for (i = 0; i < Restart; i++) fasp_mem_free(h[i]);
    fasp_mem_free(work);
    fasp_mem_free(norms);
    
    if ( iter >= MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/*---------------------------------*/
/*--    Private Functions        --*/
/*---------------------------------*/

/**
 * \fn static void dense_aAxpby (INT n, INT m, REAL *A, REAL alpha,
 *                               REAL *x, REAL beta, REAL *y)
 *
 * \brief  y = alpha*A^Tx + beta*y
 *
 * \param n     Pointer to row
 * \param m     Pointer to col
 * \param A     Pointer to CSR matrix
 * \param alpha Real factor alpha
 * \param x     Pointer to the dvector of right hand side
 * \param beta  Real factor beta
 * \param y     Maximal number of iterations
 *
 * \author zheng Li, Chensong Zhang
 * \date   12/23/2014
 */
static void dense_aAxpby (INT n,
                          INT m,
                          REAL *A,
                          REAL alpha,
                          REAL *x,
                          REAL beta,
                          REAL *y)
{
    INT i, j;
    
    for (i=0; i<m; i++) fasp_blas_array_ax(n, x[i], &A[i*n]);
    
    for (j=1; j<m; j++) {
        for (i=0; i<n; i++) {
            A[i] += A[i+j*n];
        }
    }
    
    fasp_blas_array_axpby(n, alpha, A, beta, y);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
