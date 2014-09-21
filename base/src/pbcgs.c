/*! \file pbcgs.c
 *
 *  \brief Krylov subspace methods -- Preconditioned BiCGstab
 *
 *  Abstract algorithm
 *
 *  PBICGStab method to solve A*x=b is to generate {x_k} to approximate x
 *
 *  Note: We generate a series of {p_k} such that V_k=span{p_1,...,p_k}. 
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
 *  Convergence check: norm(r)/norm(b) < tol  
 *  
 *  Stagnation check:    
 *      - IF norm(alpha*p_k)/norm(x_{k+1}) < tol_stag 
 *          -# compute r=b-A*x_{k+1}; 
 *          -# convergence check; 
 *          -# IF ( not converged & restart_number < Max_Stag_Check ) restart;
 *      - END IF  
 *  
 *  Residual check:     
 *      - IF norm(r_{k+1})/norm(b) < tol             
 *          -# compute the real residual r = b-A*x_{k+1}; 
 *          -# convergence check; 
 *          -# IF ( not converged & restart_number < Max_Res_Check ) restart;
 *      - END IF 
 *
 *  \note Refer to Y. Saad 2003
 *        Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 *  \note See spbcgs.c for a safer version
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
 * \fn INT fasp_solver_dcsr_pbcgs (dCSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, 
 *                                 const SHORT stop_type, const SHORT print_level)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b 
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Chensong Zhang
 * \date   09/09/2009
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 * Modified by Feiteng Huang on 06/01/2012: fix restart param-init 
 * Modified by Chensong Zhang on 03/31/2013
 */
INT fasp_solver_dcsr_pbcgs (dCSRmat *A, 
                            dvector *b, 
                            dvector *u, 
                            precond *pc, 
                            const REAL tol,
                            const INT MaxIt, 
                            const SHORT stop_type,
                            const SHORT print_level)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         alpha, beta, omega, temp1, temp2;
    REAL         relres, normr0, reldiff;
    REAL         factor, absres, absres0, normd, normu, tempr, infnormu;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));    
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif

    // r = b-A*u
    fasp_array_cp(m,bval,r);
    fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_array_norm2(m,r);

    // compute initial relative residual
    switch (stop_type) {
        case STOP_REL_RES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_REL_PRECRES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_MOD_REL_RES:
            normu  = MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unrecognized stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if (relres<tol) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,iter,relres,absres0,0.0);

    // rho = r* := r
    fasp_array_cp(m,r,rho);    
    temp1 = fasp_blas_array_dotprod(m,r,rho);

    // p = r
    fasp_array_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dcsr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_array_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL ) {
            alpha = temp1/temp2;
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_array_cp(m,r,s);
        fasp_blas_array_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL ) 
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else 
            fasp_array_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dcsr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_array_dotprod(m,t,t);
        omega = fasp_blas_array_dotprod(m,s,t)/tempr;
        
        // delu = alpha pp + omega sp
        fasp_blas_array_axpby(m,alpha,pp,omega,sp);
        
        // u = u + delu
        fasp_blas_array_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_array_axpy(m,-omega,t,s);
        fasp_array_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_array_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL || ABS(omega) > SMALLREAL ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else { // Possible breakdown?
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_array_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_array_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_array_norm2(m,sp);
        normu   = fasp_blas_array_norm2(m,uval);
        reldiff = normd/normu;
        
        if ( normd < TOL_s ) { // Possible breakdown?
            ITS_SMALLSP; goto FINISHED;
        }

        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_array_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reducation factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(print_level,stop_type,iter,relres,absres,factor);
        
        // Check I: if soultion is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_array_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( print_level > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if staggenated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( print_level >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL ) 
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else 
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);    
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = absres/normr0; 
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( print_level > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            if ( print_level >= PRINT_MORE ) ITS_COMPRES(relres);
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL ) 
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else 
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);    
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;                    
                case STOP_REL_PRECRES:
                    if ( pc != NULL ) 
                        pc->fct(r,z,pc->data);
                    else 
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
             }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( more_step >= MaxRestartStep ) {
                if ( print_level > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( print_level > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end if safe guard
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish the iterative method
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if ( iter > MaxIt ) 
        return ERROR_SOLVER_MAXIT;
    else 
        return iter;
}

/**
 * \fn INT fasp_solver_dbsr_pbcgs (dBSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, 
 *                                 const SHORT stop_type, const SHORT print_level) 
 *
 * \brief Preconditioned BiCGstab method for solving Au=b 
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Chensong Zhang
 * \date   09/09/2009
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 * Modified by Feiteng Huang on 06/01/2012: fix restart param-init 
 * Modified by Chensong Zhang on 03/31/2013
 *
 */
INT fasp_solver_dbsr_pbcgs (dBSRmat *A,
                            dvector *b,
                            dvector *u,
                            precond *pc,
                            const REAL tol,
                            const INT MaxIt,
                            const SHORT stop_type,
                            const SHORT print_level) 
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         alpha, beta, omega, temp1, temp2;
    REAL         relres, normr0, reldiff;
    REAL         factor, absres, absres0, normd, normu, tempr, infnormu;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_array_norm2(m,r);
    
    // compute initial relative residual
    switch (stop_type) {
        case STOP_REL_RES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_REL_PRECRES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_MOD_REL_RES:
            normu  = MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unrecognized stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if (relres<tol) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,iter,relres,absres0,0.0);
    
    // rho = r* := r
    fasp_array_cp(m,r,rho);
    temp1 = fasp_blas_array_dotprod(m,r,rho);
    
    // p = r
    fasp_array_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dbsr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_array_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL ) {
            alpha = temp1/temp2;
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_array_cp(m,r,s);
        fasp_blas_array_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dbsr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_array_dotprod(m,t,t);
        omega = fasp_blas_array_dotprod(m,s,t)/tempr;
        
        // delu = alpha pp + omega sp
        fasp_blas_array_axpby(m,alpha,pp,omega,sp);
        
        // u = u + delu
        fasp_blas_array_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_array_axpy(m,-omega,t,s);
        fasp_array_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_array_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL || ABS(omega) > SMALLREAL ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_array_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_array_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_array_norm2(m,sp);
        normu   = fasp_blas_array_norm2(m,uval);
        reldiff = normd/normu;
        
        if ( normd < TOL_s ) {
            ITS_SMALLSP; goto FINISHED;
        }
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_array_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reducation factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(print_level,stop_type,iter,relres,absres,factor);
        
        // Check I: if soultion is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_array_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( print_level > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if staggenated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( print_level >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( print_level > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            if ( print_level >= PRINT_MORE ) ITS_COMPRES(relres);
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( more_step >= MaxRestartStep ) {
                if ( print_level > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( print_level > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end if safe guard
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish the iterative method
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if ( iter > MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/**
 * \fn INT fasp_solver_bdcsr_pbcgs (block_dCSRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                  const REAL tol, const INT MaxIt, 
 *                                  const SHORT stop_type, const SHORT print_level)
 *
 * \brief A preconditioned BiCGstab method for solving Au=b 
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Xiaozhe Hu
 * \date   05/24/2010
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 * Modified by Feiteng Huang on 06/01/2012: fix restart param-init 
 * Modified by Chensong Zhang on 03/31/2013
 *
 */
INT fasp_solver_bdcsr_pbcgs (block_dCSRmat *A, 
                             dvector *b, 
                             dvector *u, 
                             precond *pc, 
                             const REAL tol,
                             const INT MaxIt, 
                             const SHORT stop_type,
                             const SHORT print_level) 
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         alpha, beta, omega, temp1, temp2;
    REAL         relres, normr0, reldiff;
    REAL         factor, absres, absres0, normd, normu, tempr, infnormu;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    fasp_blas_bdcsr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_array_norm2(m,r);
    
    // compute initial relative residual
    switch (stop_type) {
        case STOP_REL_RES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_REL_PRECRES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_MOD_REL_RES:
            normu  = MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unrecognized stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if (relres<tol) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,iter,relres,absres0,0.0);
    
    // rho = r* := r
    fasp_array_cp(m,r,rho);
    temp1 = fasp_blas_array_dotprod(m,r,rho);
    
    // p = r
    fasp_array_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_bdcsr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_array_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL ) {
            alpha = temp1/temp2;
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_array_cp(m,r,s);
        fasp_blas_array_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_bdcsr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_array_dotprod(m,t,t);
        omega = fasp_blas_array_dotprod(m,s,t)/tempr;
        
        // delu = alpha pp + omega sp
        fasp_blas_array_axpby(m,alpha,pp,omega,sp);
        
        // u = u + delu
        fasp_blas_array_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_array_axpy(m,-omega,t,s);
        fasp_array_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_array_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL || ABS(omega) > SMALLREAL ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_array_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_array_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_array_norm2(m,sp);
        normu   = fasp_blas_array_norm2(m,uval);
        reldiff = normd/normu;
        
        if ( normd < TOL_s ) {
            ITS_SMALLSP; goto FINISHED;
        }
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_array_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reducation factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(print_level,stop_type,iter,relres,absres,factor);
        
        // Check I: if soultion is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_array_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( print_level > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if staggenated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( print_level >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_bdcsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( print_level > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            if ( print_level >= PRINT_MORE ) ITS_COMPRES(relres);
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_bdcsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( more_step >= MaxRestartStep ) {
                if ( print_level > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( print_level > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end if safe guard
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish the iterative method
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if ( iter > MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/**
 * \fn INT fasp_solver_dstr_pbcgs (dSTRmat *A, dvector *b, dvector *u, precond *pc, 
 *                                 const REAL tol, const INT MaxIt, 
 *                                 const SHORT stop_type, const SHORT print_level)
 * 
 * \brief Preconditioned BiCGstab method for solving Au=b 
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param print_level  How much information to print out
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou
 * \date   04/25/2010
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 * Modified by Feiteng Huang on 06/01/2012: fix restart param-init 
 * Modified by Chensong Zhang on 03/31/2013
 *
 */
INT fasp_solver_dstr_pbcgs (dSTRmat *A, 
                            dvector *b, 
                            dvector *u, 
                            precond *pc, 
                            const REAL tol,
                            const INT MaxIt, 
                            const SHORT stop_type,
                            const SHORT print_level) 
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         alpha, beta, omega, temp1, temp2;
    REAL         relres, normr0, reldiff;
    REAL         factor, absres, absres0, normd, normu, tempr, infnormu;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_array_norm2(m,r);
    
    // compute initial relative residual
    switch (stop_type) {
        case STOP_REL_RES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_REL_PRECRES:
            normr0 = MAX(SMALLREAL,absres0);
            relres = absres0/normr0;
            break;
        case STOP_MOD_REL_RES:
            normu  = MAX(SMALLREAL,fasp_blas_array_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unrecognized stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if (relres<tol) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(print_level,stop_type,iter,relres,absres0,0.0);
    
    // rho = r* := r
    fasp_array_cp(m,r,rho);
    temp1 = fasp_blas_array_dotprod(m,r,rho);
    
    // p = r
    fasp_array_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dstr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_array_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL ) {
            alpha = temp1/temp2;
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_array_cp(m,r,s);
        fasp_blas_array_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_array_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dstr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_array_dotprod(m,t,t);
        omega = fasp_blas_array_dotprod(m,s,t)/tempr;
        
        // delu = alpha pp + omega sp
        fasp_blas_array_axpby(m,alpha,pp,omega,sp);
        
        // u = u + delu
        fasp_blas_array_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_array_axpy(m,-omega,t,s);
        fasp_array_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_array_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL || ABS(omega) > SMALLREAL ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else {
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_array_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_array_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_array_norm2(m,sp);
        normu   = fasp_blas_array_norm2(m,uval);
        reldiff = normd/normu;
        
        if ( normd < TOL_s ) {
            ITS_SMALLSP; goto FINISHED;
        }
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_array_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_array_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reducation factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(print_level,stop_type,iter,relres,absres,factor);
        
        // Check I: if soultion is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_array_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( print_level > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if staggenated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( print_level >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( print_level > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            if ( print_level >= PRINT_MORE ) ITS_COMPRES(relres);
            
            // re-init iteration param
            fasp_array_cp(m,bval,r);
            fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_array_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_array_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_array_cp(m,r,rho);
            temp1 = fasp_blas_array_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_array_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_array_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_array_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( print_level >= PRINT_MORE ) ITS_REALRES(relres);
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( more_step >= MaxRestartStep ) {
                if ( print_level > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( print_level > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end if safe guard
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish the iterative method
    if ( print_level > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if ( iter > MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
