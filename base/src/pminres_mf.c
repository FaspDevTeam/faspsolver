/*! \file pminres_mf.c
 *
 *  \brief Krylov subspace methods -- Preconditioned minimal residual (matrix free)
 *
 *  Abstract algorithm of Krylov method    
 *
 *  Krylov method to solve A*x=b is to generate {x_k} to approximate x, 
 *  where x_k is the optimal solution in Krylov space 
 *
 *     V_k=span{r_0,A*r_0,A^2*r_0,...,A^{k-1}*r_0}, 
 *
 *  under some inner product. 
 *
 *  For the implementation, we generate a series of {p_k} such that V_k=span{p_1,...,p_k}. Details: 
 *
 *  Step 0. Given A, b, x_0, M  
 *  
 *  Step 1. Compute residual r_0 = b-A*x_0 and convergence check;  
 *  
 *  Step 2. Initialization z_0 = M^{-1}*r_0, p_0=z_0;  
 *  
 *  Step 3. Main loop ...
 *
 *  FOR k = 0:MaxIt      
 *      - get step size alpha = f(r_k,z_k,p_k);      
 *      - update solution: x_{k+1} = x_k + alpha*p_k;      
 *      - perform stagnation check;      
 *      - update residual: r_{k+1} = r_k - alpha*(A*p_k);        
 *      - perform residual check;      
 *      - obtain p_{k+1} using {p_0, p_1, ... , p_k};      
 *      - prepare for next iteration;      
 *      - print the result of k-th iteration; 
 *  END FOR
 * 
 *  Convergence check is: norm(r)/norm(b) < tol  
 *  
 *  Stagnation check is like following:    
 *      - IF norm(alpha*p_k)/norm(x_{k+1}) < tol_stag 
 *          -# compute r=b-A*x_{k+1}; 
 *          -# convergence check; 
 *          -# IF ( not converged & restart_number < Max_Stag_Check ) restart;
 *      - END IF  
 *  
 *  Residual check is like following:     
 *      - IF norm(r_{k+1})/norm(b) < tol             
 *          -# compute the real residual r = b-A*x_{k+1}; 
 *          -# convergence check; 
 *          -# IF ( not converged & restart_number < Max_Res_Check ) restart;
 *      - END IF 
 *
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 */  

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_pminres (mxv_matfree *mf, dvector *b, dvector *u, precond *pc, 
 *                              const REAL tol, const INT MaxIt, 
 *                              const SHORT stop_type, const SHORT print_level) 
 *
 * \brief A preconditioned minimal residual (Minres) method for solving Au=b 
 *
 * \param mf           Pointer to mxv_matfree: the spmv operation
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
 * \author Shiquan Zhang
 * \date   10/24/2010
 *
 * Rewritten by Chensong Zhang on 05/01/2012
 * Modified by Feiteng Huang on 09/26/2012: matrix free
 */
INT fasp_solver_pminres (mxv_matfree *mf, 
                         dvector *b, 
                         dvector *u, 
                         precond *pc, 
                         const REAL tol,
                         const INT MaxIt, 
                         const SHORT stop_type,
                         const SHORT print_level) 
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m=b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance

    // local variables
    INT          iter=0, stag, more_step, restart_step;
    REAL         absres0, absres, relres, factor;
    REAL         alpha, alpha0, alpha1, temp2; 
    REAL         normr0, normu2, normuu, normp, infnormu;
    
    // allocate temp memory (need 11*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(11*m,sizeof(REAL));    
    REAL *p0=work, *p1=work+m, *p2=p1+m, *z0=p2+m, *z1=z0+m;
    REAL *t0=z1+m, *t1=t0+m, *t=t1+m, *tp=t+m, *tz=tp+m, *r=tz+m;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif

    // initialization counters
    stag=1; more_step=1; restart_step=1;
    
    // p0=0
    fasp_array_set(m,p0,0.0);
    
    // r = b-A*u
    mf->fct(mf->data, u->val, r);
    fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
    
    // p1 = B(r)
    if (pc != NULL) 
        pc->fct(r,p1,pc->data); /* Apply preconditioner */
    else 
        fasp_array_cp(m,r,p1); /* No preconditioner */

    // compute initial relative residual 
    switch (stop_type) {
    case STOP_REL_PRECRES:
        absres0=sqrt(ABS(fasp_blas_array_dotprod(m,r,p1)));
        normr0=MAX(SMALLREAL,absres0);
        relres=absres0/normr0; 
        break;
    case STOP_MOD_REL_RES:
        absres0=fasp_blas_array_norm2(m,r);
        normu2=MAX(SMALLREAL,fasp_blas_array_norm2(m,u->val));
        relres=absres0/normu2; 
        break;
    default: // STOP_REL_RES
        absres0=fasp_blas_array_norm2(m,r);
        normr0=MAX(SMALLREAL,absres0);
        relres=absres0/normr0; 
        break;
    }
    
    if (relres<tol) goto FINISHED;
    
    // tp=A*p1
    mf->fct(mf->data, p1, tp);
    
    // tz = B(tp)
    if (pc != NULL) 
        pc->fct(tp,tz,pc->data); /* Apply preconditioner */
    else 
        fasp_array_cp(m,tp,tz); /* No preconditioner */
    
    // p1=p1/normp
    normp=ABS(fasp_blas_array_dotprod(m,tz,tp));
    normp=sqrt(normp);    
    fasp_array_cp(m,p1,t);
    fasp_array_set(m,p1,0.0);
    fasp_blas_array_axpy(m,1/normp,t,p1);
    
    // t0=A*p0=0
    fasp_array_set(m,t0,0.0);
    fasp_array_cp(m,t0,z0);
    fasp_array_cp(m,t0,t1);
    fasp_array_cp(m,t0,z1);
    
    // t1=tp/normp,z1=tz/normp
    fasp_blas_array_axpy(m,1.0/normp,tp,t1);
    fasp_blas_array_axpy(m,1.0/normp,tz,z1);
    
    while( iter++ < MaxIt) {
    
        // alpha=<r,z1>
        alpha=fasp_blas_array_dotprod(m,r,z1);
    
        // u=u+alpha*p1
        fasp_blas_array_axpy(m,alpha,p1,u->val);
    
        // r=r-alpha*Ap1
        fasp_blas_array_axpy(m,-alpha,t1,r);
    
        // compute t=A*z1 alpha1=<z1,t> 
        mf->fct(mf->data, z1, t);
        alpha1=fasp_blas_array_dotprod(m,z1,t);
    
        // compute t=A*z0 alpha0=<z1,t> 
        mf->fct(mf->data, z0, t);
        alpha0=fasp_blas_array_dotprod(m,z1,t);
    
        // p2=z1-alpha1*p1-alpha0*p0
        fasp_array_cp(m,z1,p2);
        fasp_blas_array_axpy(m,-alpha1,p1,p2);
        fasp_blas_array_axpy(m,-alpha0,p0,p2);
    
        // tp=A*p2
        mf->fct(mf->data, p2, tp);
    
        // tz = B(tp)
        if (pc != NULL) 
            pc->fct(tp,tz,pc->data); /* Apply preconditioner */
        else 
            fasp_array_cp(m,tp,tz); /* No preconditioner */
    
        // p2=p2/normp
        normp=ABS(fasp_blas_array_dotprod(m,tz,tp));
        normp=sqrt(normp);
        fasp_array_cp(m,p2,t);
        fasp_array_set(m,p2,0.0);
        fasp_blas_array_axpy(m,1/normp,t,p2);
    
        // prepare for the next iteration
        fasp_array_cp(m,p1,p0);
        fasp_array_cp(m,p2,p1);
        fasp_array_cp(m,t1,t0);
        fasp_array_cp(m,z1,z0);
    
        // t1=tp/normp,z1=tz/normp
        fasp_array_set(m,t1,0.0);
        fasp_array_cp(m,t1,z1);
        fasp_blas_array_axpy(m,1/normp,tp,t1);
        fasp_blas_array_axpy(m,1/normp,tz,z1);
    
        // relative residual = ||r||/||r0||
        temp2=fasp_blas_array_dotprod(m,r,r);
        absres=sqrt(temp2);    
    
        normu2=fasp_blas_array_norm2(m,u->val);
    
        switch (stop_type) {
        case STOP_REL_PRECRES:
            if (pc == NULL)
                fasp_array_cp(m,r,t);
            else
                pc->fct(r,t,pc->data);
            temp2=ABS(fasp_blas_array_dotprod(m,r,t));
            relres=sqrt(temp2)/normr0; 
            break;
        case STOP_MOD_REL_RES:
            relres=sqrt(temp2)/normu2; 
            break;
        default: // STOP_REL_RES
            relres=sqrt(temp2)/normr0; 
            break;
        }
    
        // compute reducation factor of residual ||r||
        factor=absres/absres0;
    
        // output iteration information if needed    
        print_itinfo(print_level,stop_type,iter,relres,absres,factor);
    
        // solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
        infnormu = fasp_blas_array_norminf(m, u->val); 
        if (infnormu <= sol_inf_tol) {
            if (print_level>PRINT_MIN) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            break;
        }
    
        normuu=fasp_blas_array_norm2(m,p1);
        normuu=ABS(alpha)*(normuu/normu2);
    
        // check convergence
        if (normuu<maxdiff) {
            if (stag<MaxStag) {
                if (print_level>=PRINT_MORE) { 
                    ITS_DIFFRES(normuu,relres);
                    ITS_RESTART;
                }
            }
    
            mf->fct(mf->data, u->val, r);
            fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
    
            temp2=fasp_blas_array_dotprod(m,r,r);
            absres=sqrt(temp2);
            switch (stop_type) {
            case STOP_REL_RES:
                relres=sqrt(temp2)/normr0; 
                break;
            case STOP_REL_PRECRES:
                if (pc == NULL)
                    fasp_array_cp(m,r,t);
                else
                    pc->fct(r,t,pc->data);
                temp2=ABS(fasp_blas_array_dotprod(m,r,t));
                relres=sqrt(temp2)/normr0; 
                break;
            case STOP_MOD_REL_RES:
                relres=sqrt(temp2)/normu2; 
                break;
            }
    
            if (print_level>=PRINT_MORE) ITS_REALRES(relres);
    
            if (relres<tol)
                break;
            else {
                if (stag>=MaxStag) {
                    if (print_level>PRINT_MIN) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    break;
                }
                ++stag;
                ++restart_step;
    
                fasp_array_set(m,p0,0.0);
    
                // p1 = B(r)
                if (pc != NULL) 
                    pc->fct(r,p1,pc->data); /* Apply preconditioner */
                else 
                    fasp_array_cp(m,r,p1); /* No preconditioner */
    
                // tp=A*p1
                mf->fct(mf->data, p1, tp);
    
                // tz = B(tp)
                if (pc == NULL) 
                    pc->fct(tp,tz,pc->data); /* Apply rreconditioner */
                else 
                    fasp_array_cp(m,tp,tz); /* No preconditioner */
    
                // p1=p1/normp
                normp=fasp_blas_array_dotprod(m,tz,tp);
                normp=sqrt(normp);
                fasp_array_cp(m,p1,t);
    
                // t0=A*p0=0
                fasp_array_set(m,t0,0.0);
                fasp_array_cp(m,t0,z0);
                fasp_array_cp(m,t0,t1);
                fasp_array_cp(m,t0,z1);
                fasp_array_cp(m,t0,p1);
    
                fasp_blas_array_axpy(m,1/normp,t,p1);
    
                // t1=tp/normp,z1=tz/normp
                fasp_blas_array_axpy(m,1/normp,tp,t1);
                fasp_blas_array_axpy(m,1/normp,tz,z1);
            }
        }
    
        // safe guard
        if (relres<tol) {
            if (print_level>=PRINT_MORE) ITS_COMPRES(relres);
    
            mf->fct(mf->data, u->val, r);
            fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
    
            temp2=fasp_blas_array_dotprod(m,r,r);
            absres=sqrt(temp2);
            switch (stop_type) {
            case STOP_REL_RES:
                relres=sqrt(temp2)/normr0; 
                break;
            case STOP_REL_PRECRES:
                if (pc == NULL)
                    fasp_array_cp(m,r,t);
                else
                    pc->fct(r,t,pc->data);
                temp2=ABS(fasp_blas_array_dotprod(m,r,t));
                relres=sqrt(temp2)/normr0; 
                break;
            case STOP_MOD_REL_RES:
                relres=sqrt(temp2)/normu2; 
                break;
            }
    
            if (print_level>=PRINT_MORE) ITS_REALRES(relres);
    
            // check convergence
            if (relres<tol) break;
    
            if (more_step>=MaxRestartStep) {
                if (print_level>PRINT_MIN) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                break;
            }
    
            if (more_step<MaxRestartStep) {
                if (print_level>PRINT_NONE) ITS_RESTART;
            }
    
            ++more_step;
            ++restart_step;
    
            fasp_array_set(m,p0,0.0);
    
            // p1 = B(r)
            if (pc != NULL) 
                pc->fct(r,p1,pc->data); /* Apply preconditioner */
            else 
                fasp_array_cp(m,r,p1); /* No preconditioner */
    
            // tp = A*p1
            mf->fct(mf->data, p1, tp);
    
            // tz = B(tp)
            if (pc == NULL) 
                pc->fct(tp,tz,pc->data); /* Apply rreconditioner */
            else 
                fasp_array_cp(m,tp,tz); /* No preconditioner */
    
            // p1 = p1/normp
            normp=fasp_blas_array_dotprod(m,tz,tp);
            normp=sqrt(normp);
            fasp_array_cp(m,p1,t);
    
            // t0=A*p0=0
            fasp_array_set(m,t0,0.0);
            fasp_array_cp(m,t0,z0);
            fasp_array_cp(m,t0,t1);
            fasp_array_cp(m,t0,z1);
            fasp_array_cp(m,t0,p1);
    
            fasp_blas_array_axpy(m,1/normp,t,p1);
    
            // t1=tp/normp,z1=tz/normp
            fasp_blas_array_axpy(m,1/normp,tp,t1);
            fasp_blas_array_axpy(m,1/normp,tz,z1);
    
        }
        // update relative residual here
        absres0 = absres;
    }
    
 FINISHED:  // finish the iterative method
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
