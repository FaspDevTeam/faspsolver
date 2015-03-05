/*! \file pcg_mf.c
 *
 *  \brief Krylov subspace methods -- Preconditioned conjugate gradient (matrix free)
 *
 *  Abstract algorithm
 *
 *  PCG method to solve A*x=b is to generate {x_k} to approximate x
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
 */

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_pcg (mxv_matfree *mf, dvector *b, dvector *u, precond *pc,
 *                          const REAL tol, const INT MaxIt,
 *                          const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned conjugate gradient (CG) method for solving Au=b
 *
 * \param mf           Pointer to mxv_matfree: the spmv operation
 * \param b            Pointer to dvector: the right hand side
 * \param u            Pointer to dvector: the unknowns
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chensong Zhang, Xiaozhe Hu, Shiquan Zhang
 * \date   05/06/2010
 *
 * Modified by Chensong Zhang on 04/30/2012
 * Modified by Feiteng Huang on 09/19/2012: matrix free
 */
INT fasp_solver_pcg (mxv_matfree *mf,
                     dvector *b,
                     dvector *u,
                     precond *pc,
                     const REAL tol,
                     const INT MaxIt,
                     const SHORT stop_type,
                     const SHORT prtlvl)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m=b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    
    // local variables
    INT          iter = 0, stag, more_step, restart_step;
    REAL         absres0 = BIGREAL, absres = BIGREAL;
    REAL         relres  = BIGREAL, normu  = BIGREAL, normr0 = BIGREAL;
    REAL         reldiff, factor, infnormu;
    REAL         alpha, beta, temp1, temp2;
    
    // allocate temp memory (need 4*m REAL numbers)
    REAL *work=(REAL *)fasp_mem_calloc(4*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // initialize counters
    stag=1; more_step=1; restart_step=1;
    
    // r = b-A*u
    mf->fct(mf->data, u->val, r);
    fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
    
    if (pc != NULL)
        pc->fct(r,z,pc->data); /* Apply preconditioner */
    else
        fasp_array_cp(m,r,z); /* No preconditioner */
    
    // compute initial relative residual
    switch (stop_type) {
        case STOP_REL_PRECRES:
            absres0=sqrt(fasp_blas_array_dotprod(m,r,z));
            normr0=MAX(SMALLREAL,absres0);
            relres=absres0/normr0;
            break;
        case STOP_MOD_REL_RES:
            absres0=fasp_blas_array_norm2(m,r);
            normu=MAX(SMALLREAL,fasp_blas_array_norm2(m,u->val));
            relres=absres0/normu;
            break;
        default:
            absres0=fasp_blas_array_norm2(m,r);
            normr0=MAX(SMALLREAL,absres0);
            relres=absres0/normr0;
            break;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol || absres0 < 1e-3*tol ) goto FINISHED;
    
    fasp_array_cp(m,z,p);
    temp1=fasp_blas_array_dotprod(m,z,r);
    
    while ( iter++ < MaxIt ) {
        
        // t=A*p
        mf->fct(mf->data, p, t);
        
        // alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
        temp2=fasp_blas_array_dotprod(m,t,p);
        alpha=temp1/temp2;
        
        // u_k=u_{k-1} + alpha_k*p_{k-1}
        fasp_blas_array_axpy(m,alpha,p,u->val);
        
        // r_k=r_{k-1} - alpha_k*A*p_{k-1}
        fasp_blas_array_axpy(m,-alpha,t,r);
        absres=fasp_blas_array_norm2(m,r);
        
        // compute reducation factor of residual ||r||
        factor=absres/absres0;
        
        // compute relative residual
        switch (stop_type) {
            case STOP_REL_PRECRES:
                // z = B(r)
                if (pc != NULL)
                    pc->fct(r,z,pc->data); /* Apply preconditioner */
                else
                    fasp_array_cp(m,r,z); /* No preconditioner */
                temp2=fasp_blas_array_dotprod(m,z,r);
                relres=sqrt(ABS(temp2))/normr0;
                break;
            case STOP_MOD_REL_RES:
                relres=absres/normu;
                break;
            default:
                relres=absres/normr0;
                break;
        }
        
        // output iteration information if needed
        print_itinfo(prtlvl,stop_type,iter,relres,absres,factor);
        
        // solution check, if soultion is too small, return ERROR_SOLVER_SOLSTAG.
        infnormu = fasp_blas_array_norminf(m, u->val);
        if ( infnormu <= sol_inf_tol ) {
            if ( prtlvl > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            break;
        }
        
        // compute relative difference
        normu=fasp_blas_dvec_norm2(u);
        reldiff=ABS(alpha)*fasp_blas_array_norm2(m,p)/normu;
        
        // stagnation check
        if ( (stag<=MaxStag) & (reldiff<maxdiff) ) {
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            mf->fct(mf->data, u->val, r);
            fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
            absres=fasp_blas_array_norm2(m,r);
            
            // relative residual
            switch (stop_type) {
                case STOP_REL_PRECRES:
                    // z = B(r)
                    if (pc != NULL)
                        pc->fct(r,z,pc->data); /* Apply preconditioner */
                    else
                        fasp_array_cp(m,r,z); /* No preconditioner */
                    temp2=fasp_blas_array_dotprod(m,z,r);
                    relres=sqrt(ABS(temp2))/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    relres=absres/normu;
                    break;
                default:
                    relres=absres/normr0;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( prtlvl > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    break;
                }
                fasp_array_set(m,p,0.0);
                ++stag;
                ++restart_step;
            }
        } // end of staggnation check!
        
        // safe-guard check
        if ( relres < tol ) {
            if ( prtlvl >= PRINT_MORE ) ITS_COMPRES(relres);
            
            mf->fct(mf->data, u->val, r);
            fasp_blas_array_axpby(m, 1.0, b->val, -1.0, r);
            
            // relative residual
            switch (stop_type) {
                case STOP_REL_PRECRES:
                    // z = B(r)
                    if (pc != NULL)
                        pc->fct(r,z,pc->data); /* Apply preconditioner */
                    else
                        fasp_array_cp(m,r,z); /* No preconditioner */
                    temp2=fasp_blas_array_dotprod(m,z,r);
                    relres=sqrt(ABS(temp2))/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres=fasp_blas_array_norm2(m,r);
                    relres=absres/normu;
                    break;
                default:
                    absres=fasp_blas_array_norm2(m,r);
                    relres=absres/normr0;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( more_step >= MaxRestartStep ) {
                if ( prtlvl > PRINT_MIN) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                break;
            }
            
            // prepare for restarting the method
            fasp_array_set(m,p,0.0);
            ++more_step;
            ++restart_step;
            
        } // end of safe-guard check!
        
        // update relative residual here
        absres0 = absres;
        
        // compute z_k = B(r_k)
        if ( stop_type != STOP_REL_PRECRES ) {
            if ( pc != NULL )
                pc->fct(r,z,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,r,z); /* No preconditioner, B=I */
        }
        
        // compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
        temp2=fasp_blas_array_dotprod(m,z,r);
        beta=temp2/temp1;
        temp1=temp2;
        
        // compute p_k = z_k + beta_k*p_{k-1}
        fasp_blas_array_axpby(m,1.0,z,beta,p);
        
    } // end of main PCG loop.
    
FINISHED:  // finish the iterative method
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE > 0
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
