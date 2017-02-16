/*! \file  KryPbcgs.c
 *
 *  \brief Krylov subspace methods -- Preconditioned BiCGstab
 *
 *  \note  This file contains Level-3 (Kry) functions. It requires:
 *         AuxArray.c, AuxMemory.c, AuxMessage.c, BlaArray.c, BlaSpmvBLC.c,
 *         BlaSpmvBSR.c, BlaSpmvCSR.c, and BlaSpmvSTR.c
 *
 *  \note  See KrySPbcgs.c for a safer version
 *
 *  Reference:
 *         Y. Saad 2003
 *         Iterative methods for sparse linear systems (2nd Edition), SIAM
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  TODO: Use one single function for all! --Chensong
 */

#include <math.h>
#include <float.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/

#include "KryUtil.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dcsr_pbcgs (dCSRmat *A, dvector *b, dvector *u, precond *pc,
 *                                 const REAL tol, const INT MaxIt,
 *                                 const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b
 *
 * \param A            Pointer to coefficient matrix
 * \param b            Pointer to dvector of right hand side
 * \param u            Pointer to dvector of DOFs
 * \param pc           Pointer to precond: structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chensong Zhang
 * \date   09/09/2009
 *
 * Modified by Feiteng Huang on 06/01/2012: fix restart param-init
 * Modified by Chensong Zhang on 03/31/2013
 */
INT fasp_solver_dcsr_pbcgs (dCSRmat     *A,
                            dvector     *b,
                            dvector     *u,
                            precond     *pc,
                            const REAL   tol,
                            const INT    MaxIt,
                            const SHORT  stop_type,
                            const SHORT  prtlvl)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // stagnation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL         normd   = BIGREAL, normu  = BIGREAL, normr0 = BIGREAL;
    REAL         reldiff, factor, infnormu;
    REAL         alpha, beta, omega, temp1, temp2, tempr;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_darray_cp(m,bval,r);
    fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_darray_norm2(m,r);
    
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
            normu  = MAX(SMALLREAL,fasp_blas_darray_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unknown stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol || absres0 < 1e-3*tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(prtlvl,stop_type,iter,relres,absres0,0.0);
    
    // shadow residual rho = r* := r
    fasp_darray_cp(m,r,rho);
    temp1 = fasp_blas_darray_dotprod(m,r,rho);
    
    // p = r
    fasp_darray_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dcsr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_darray_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL2 ) {
            alpha = temp1/temp2;
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_darray_cp(m,r,s);
        fasp_blas_darray_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dcsr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_darray_dotprod(m,t,t);
        omega = fasp_blas_darray_dotprod(m,s,t)/tempr;
        
        // diffu = alpha pp + omega sp
        fasp_blas_darray_axpby(m,alpha,pp,omega,sp);
        
        // u = u + diffu
        fasp_blas_darray_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_darray_axpy(m,-omega,t,s);
        fasp_darray_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_darray_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL2 || ABS(omega) > SMALLREAL2 ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_darray_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_darray_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_darray_norm2(m,sp);
        if ( normd < TOL_s ) { // Possible breakdown?
            ITS_SMALLSP; goto FINISHED;
        }
        
        normu   = fasp_blas_darray_norm2(m,uval);
        reldiff = normd/normu;
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_darray_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reduction factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(prtlvl,stop_type,iter,relres,absres,factor);
        
        // Check I: if solution is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_darray_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( prtlvl > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if stagnated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_darray_cp(m,bval,r);
            fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( prtlvl > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            
            REAL computed_relres = relres;
            
            // compute current residual
            fasp_darray_cp(m,bval,r);
            fasp_blas_dcsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_COMPRES(computed_relres); ITS_REALRES(relres);
            }
            
            if ( more_step >= MaxRestartStep ) {
                if ( prtlvl > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( prtlvl > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end of false convergence check
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish iterative method
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE > 0
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
 *                                 const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b
 *
 * \param A            Pointer to coefficient matrix
 * \param b            Pointer to dvector of right hand side
 * \param u            Pointer to dvector of DOFs
 * \param pc           Pointer to precond: structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chensong Zhang
 * \date   09/09/2009
 *
 * Modified by Feiteng Huang on 06/01/2012: fix restart param-init
 * Modified by Chensong Zhang on 03/31/2013
 */
INT fasp_solver_dbsr_pbcgs (dBSRmat     *A,
                            dvector     *b,
                            dvector     *u,
                            precond     *pc,
                            const REAL   tol,
                            const INT    MaxIt,
                            const SHORT  stop_type,
                            const SHORT  prtlvl)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // stagnation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL         normd   = BIGREAL, normu  = BIGREAL, normr0 = BIGREAL;
    REAL         reldiff, factor, infnormu;
    REAL         alpha, beta, omega, temp1, temp2, tempr;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_darray_cp(m,bval,r);
    fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_darray_norm2(m,r);
    
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
            normu  = MAX(SMALLREAL,fasp_blas_darray_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unknown stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol || absres0 < 1e-3*tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(prtlvl,stop_type,iter,relres,absres0,0.0);
    
    // shadow residual rho = r* := r
    fasp_darray_cp(m,r,rho);
    temp1 = fasp_blas_darray_dotprod(m,r,rho);
    
    // p = r
    fasp_darray_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dbsr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_darray_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL2 ) {
            alpha = temp1/temp2;
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_darray_cp(m,r,s);
        fasp_blas_darray_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dbsr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_darray_dotprod(m,t,t);
        omega = fasp_blas_darray_dotprod(m,s,t)/tempr;
        
        // diffu = alpha pp + omega sp
        fasp_blas_darray_axpby(m,alpha,pp,omega,sp);
        
        // u = u + diffu
        fasp_blas_darray_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_darray_axpy(m,-omega,t,s);
        fasp_darray_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_darray_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL2 || ABS(omega) > SMALLREAL2 ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_darray_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_darray_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_darray_norm2(m,sp);
        if ( normd < TOL_s ) { // Possible breakdown?
            ITS_SMALLSP; goto FINISHED;
        }
        
        normu   = fasp_blas_darray_norm2(m,uval);
        reldiff = normd/normu;
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_darray_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reduction factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(prtlvl,stop_type,iter,relres,absres,factor);
        
        // Check I: if solution is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_darray_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( prtlvl > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if stagnated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_darray_cp(m,bval,r);
            fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( prtlvl > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            
            REAL computed_relres = relres;
            
            // compute current residual
            fasp_darray_cp(m,bval,r);
            fasp_blas_dbsr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_COMPRES(computed_relres); ITS_REALRES(relres);
            }
            
            if ( more_step >= MaxRestartStep ) {
                if ( prtlvl > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( prtlvl > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end of false convergence check
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish iterative method
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if ( iter > MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/**
 * \fn INT fasp_solver_dblc_pbcgs (dBLCmat *A, dvector *b, dvector *u, precond *pc,
 *                                 const REAL tol, const INT MaxIt,
 *                                 const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief A preconditioned BiCGstab method for solving Au=b
 *
 * \param A            Pointer to coefficient matrix
 * \param b            Pointer to dvector of right hand side
 * \param u            Pointer to dvector of DOFs
 * \param pc           Pointer to precond: structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   05/24/2010
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 */
INT fasp_solver_dblc_pbcgs (dBLCmat     *A,
                            dvector     *b,
                            dvector     *u,
                            precond     *pc,
                            const REAL   tol,
                            const INT    MaxIt,
                            const SHORT  stop_type,
                            const SHORT  prtlvl)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // stagnation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL         normd   = BIGREAL, normu  = BIGREAL, normr0 = BIGREAL;
    REAL         reldiff, factor, infnormu;
    REAL         alpha, beta, omega, temp1, temp2, tempr;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_darray_cp(m,bval,r);
    fasp_blas_dblc_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_darray_norm2(m,r);
    
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
            normu  = MAX(SMALLREAL,fasp_blas_darray_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unknown stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol || absres0 < 1e-3*tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(prtlvl,stop_type,iter,relres,absres0,0.0);
    
    // shadow residual rho = r* := r
    fasp_darray_cp(m,r,rho);
    temp1 = fasp_blas_darray_dotprod(m,r,rho);
    
    // p = r
    fasp_darray_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dblc_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_darray_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL2 ) {
            alpha = temp1/temp2;
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_darray_cp(m,r,s);
        fasp_blas_darray_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dblc_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_darray_dotprod(m,t,t);
        omega = fasp_blas_darray_dotprod(m,s,t)/tempr;
        
        // diffu = alpha pp + omega sp
        fasp_blas_darray_axpby(m,alpha,pp,omega,sp);
        
        // u = u + diffu
        fasp_blas_darray_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_darray_axpy(m,-omega,t,s);
        fasp_darray_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_darray_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL2 || ABS(omega) > SMALLREAL2 ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_darray_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_darray_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_darray_norm2(m,sp);
        if ( normd < TOL_s ) { // Possible breakdown?
            ITS_SMALLSP; goto FINISHED;
        }
        
        normu   = fasp_blas_darray_norm2(m,uval);
        reldiff = normd/normu;
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_darray_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reduction factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(prtlvl,stop_type,iter,relres,absres,factor);
        
        // Check I: if solution is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_darray_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( prtlvl > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if stagnated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_darray_cp(m,bval,r);
            fasp_blas_dblc_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( prtlvl > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            
            REAL computed_relres = relres;
            
            // compute current residual
            fasp_darray_cp(m,bval,r);
            fasp_blas_dblc_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_COMPRES(computed_relres); ITS_REALRES(relres);
            }
            
            if ( more_step >= MaxRestartStep ) {
                if ( prtlvl > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( prtlvl > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end of false convergence check
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish iterative method
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE > 0
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
 *                                 const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b
 *
 * \param A            Pointer to coefficient matrix
 * \param b            Pointer to dvector of right hand side
 * \param u            Pointer to dvector of DOFs
 * \param pc           Pointer to precond: structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Zhiyang Zhou
 * \date   04/25/2010
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 */
INT fasp_solver_dstr_pbcgs (dSTRmat     *A,
                            dvector     *b,
                            dvector     *u,
                            precond     *pc,
                            const REAL   tol,
                            const INT    MaxIt,
                            const SHORT  stop_type,
                            const SHORT  prtlvl)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // stagnation tolerance
    const REAL   sol_inf_tol = SMALLREAL; // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2; // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL         normd   = BIGREAL, normu  = BIGREAL, normr0 = BIGREAL;
    REAL         reldiff, factor, infnormu;
    REAL         alpha, beta, omega, temp1, temp2, tempr;
    REAL         *uval=u->val, *bval=b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_darray_cp(m,bval,r);
    fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
    absres0 = fasp_blas_darray_norm2(m,r);
    
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
            normu  = MAX(SMALLREAL,fasp_blas_darray_norm2(m,uval));
            relres = absres0/normu;
            break;
        default:
            printf("### ERROR: Unknown stopping type for %s!\n", __FUNCTION__);
            goto FINISHED;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol || absres0 < 1e-3*tol ) goto FINISHED;
    
    // output iteration information if needed
    print_itinfo(prtlvl,stop_type,iter,relres,absres0,0.0);
    
    // shadow residual rho = r* := r
    fasp_darray_cp(m,r,rho);
    temp1 = fasp_blas_darray_dotprod(m,r,rho);
    
    // p = r
    fasp_darray_cp(m,r,p);
    
    // main BiCGstab loop
    while ( iter++ < MaxIt ) {
        
        // pp = precond(p)
        if ( pc != NULL )
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,p,pp); /* No preconditioner */
        
        // z = A*pp
        fasp_blas_dstr_mxv(A,pp,z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2 = fasp_blas_darray_dotprod(m,z,rho);
        if ( ABS(temp2) > SMALLREAL2 ) {
            alpha = temp1/temp2;
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // s = r - alpha z
        fasp_darray_cp(m,r,s);
        fasp_blas_darray_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if ( pc != NULL )
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        fasp_blas_dstr_mxv(A,sp,t);
        
        // omega = (t,s)/(t,t)
        tempr = fasp_blas_darray_dotprod(m,t,t);
        omega = fasp_blas_darray_dotprod(m,s,t)/tempr;
        
        // diffu = alpha pp + omega sp
        fasp_blas_darray_axpby(m,alpha,pp,omega,sp);
        
        // u = u + diffu
        fasp_blas_darray_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_darray_axpy(m,-omega,t,s);
        fasp_darray_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2 = temp1;
        temp1 = fasp_blas_darray_dotprod(m,r,rho);
        
        if ( ABS(temp2) > SMALLREAL2 || ABS(omega) > SMALLREAL2 ) {
            beta = (temp1*alpha)/(temp2*omega);
        }
        else { // Possible breakdown
            ITS_DIVZERO; goto FINISHED;
        }
        
        // p = p - omega z
        fasp_blas_darray_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_darray_axpby(m,1.0,r,beta,p);
        
        // compute difference
        normd   = fasp_blas_darray_norm2(m,sp);
        if ( normd < TOL_s ) { // Possible breakdown?
            ITS_SMALLSP; goto FINISHED;
        }
        
        normu   = fasp_blas_darray_norm2(m,uval);
        reldiff = normd/normu;
        
        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normr0;
                break;
            case STOP_REL_PRECRES:
                if ( pc == NULL )
                    fasp_darray_cp(m,r,z);
                else
                    pc->fct(r,z,pc->data);
                absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                relres = absres/normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = fasp_blas_darray_norm2(m,r);
                relres = absres/normu;
                break;
        }
        
        // compute reduction factor of residual ||r||
        factor = absres/absres0;
        
        // output iteration information if needed
        print_itinfo(prtlvl,stop_type,iter,relres,absres,factor);
        
        // Check I: if solution is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = fasp_blas_darray_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( prtlvl > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // Check II: if stagnated, try to restart
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            fasp_darray_cp(m,bval,r);
            fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = absres/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            if ( relres < tol )
                break;
            else {
                if ( stag >= MaxStag ) {
                    if ( prtlvl > PRINT_MIN ) ITS_STAGGED;
                    iter = ERROR_SOLVER_STAG;
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // Check III: prevent false convergence
        if ( relres < tol ) {
            
            REAL computed_relres = relres;
            
            // compute current residual
            fasp_darray_cp(m,bval,r);
            fasp_blas_dstr_aAxpy(-1.0,A,uval,r);
            
            // pp = precond(p)
            fasp_darray_cp(m,r,p);
            if ( pc != NULL )
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1 = fasp_blas_darray_dotprod(m,r,rho);
            
            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normr0;
                    break;
                case STOP_REL_PRECRES:
                    if ( pc != NULL )
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    absres = sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres = tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    absres = fasp_blas_darray_norm2(m,r);
                    relres = absres/normu;
                    break;
            }
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_COMPRES(computed_relres); ITS_REALRES(relres);
            }
            
            if ( more_step >= MaxRestartStep ) {
                if ( prtlvl > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( prtlvl > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end of false convergence check
        
        absres0 = absres;
        
    } // end of main BiCGstab loop
    
FINISHED:  // finish iterative method
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    // clean up temp memory
    fasp_mem_free(work);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    if ( iter > MaxIt )
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/**
 * \fn INT fasp_solver_pbcgs (mxv_matfree *mf, dvector *b, dvector *u, precond *pc,
 *                            const REAL tol, const INT MaxIt,
 *                            const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b
 *
 * \param mf           Pointer to mxv_matfree: spmv operation
 * \param b            Pointer to dvector: right hand side
 * \param u            Pointer to dvector: unknowns
 * \param pc           Pointer to precond: structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chensong Zhang
 * \date   09/09/2009
 *
 * Rewritten by Chensong Zhang on 04/30/2012
 */
INT fasp_solver_pbcgs (mxv_matfree  *mf,
                       dvector      *b,
                       dvector      *u,
                       precond      *pc,
                       const REAL    tol,
                       const INT     MaxIt,
                       const SHORT   stop_type,
                       const SHORT   prtlvl)
{
    const SHORT  MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    const INT    m = b->row;
    const REAL   maxdiff = tol*STAG_RATIO; // staganation tolerance
    const REAL   sol_inf_tol = SMALLREAL;  // infinity norm tolerance
    const REAL   TOL_s = tol*1e-2;         // tolerance for norm(p)
    
    // local variables
    INT          iter = 0, stag = 1, more_step = 1, restart_step = 1;
    REAL         absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL         normd   = BIGREAL, normu  = BIGREAL, normr0 = BIGREAL;
    REAL         reldiff, factor, infnormu;
    REAL         alpha, beta, omega, temp1, temp2, tempr;
    REAL         *uval = u->val, *bval = b->val;
    
    // allocate temp memory (need 8*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(8*m,sizeof(REAL));
    REAL *p=work, *z=work+m, *r=z+m, *t=r+m;
    REAL *rho=t+m, *pp=rho+m, *s=pp+m, *sp=s+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // initialize counters
    stag=more_step=restart_step=1;
    
    // r = b-A*u
    mf->fct(mf->data, uval, r);
    fasp_blas_darray_axpby(m, 1.0, bval, -1.0, r);
    
    // pp=precond(p)
    fasp_darray_cp(m,r,p);
    if (pc != NULL)
        pc->fct(p,pp,pc->data); /* Apply preconditioner */
    else
        fasp_darray_cp(m,p,pp); /* No preconditioner */
    
    // compute initial relative residual
    switch (stop_type) {
        case STOP_REL_PRECRES:
            absres0=sqrt(ABS(fasp_blas_darray_dotprod(m,r,pp)));
            normr0=MAX(SMALLREAL,absres0);
            relres=absres0/normr0;
            break;
        case STOP_MOD_REL_RES:
            absres0=fasp_blas_darray_norm2(m,r);
            normu=MAX(SMALLREAL,fasp_blas_darray_norm2(m,uval));
            relres=absres0/normu;
            break;
        default:
            absres0=fasp_blas_darray_norm2(m,r);
            normr0=MAX(SMALLREAL,absres0);
            relres=absres0/normr0;
            break;
    }
    
    // if initial residual is small, no need to iterate!
    if ( relres < tol || absres0 < 1e-3*tol ) goto FINISHED;
    
    // rho = r* := r
    fasp_darray_cp(m,r,rho);
    temp1=fasp_blas_darray_dotprod(m,r,rho);
    
    while ( iter++ < MaxIt ) {
        
        // z = A*pp
        mf->fct(mf->data, pp, z);
        
        // alpha = (r,rho)/(A*p,rho)
        temp2=fasp_blas_darray_dotprod(m,z,rho);
        if (ABS(temp2)>SMALLREAL) {
            alpha=temp1/temp2;
        }
        else {
            ITS_DIVZERO;
            return ERROR_SOLVER_MISC;
        }
        
        // s = r - alpha z
        fasp_darray_cp(m,r,s);
        fasp_blas_darray_axpy(m,-alpha,z,s);
        
        // sp = precond(s)
        if (pc != NULL)
            pc->fct(s,sp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,s,sp); /* No preconditioner */
        
        // t = A*sp;
        mf->fct(mf->data, sp, t);
        
        // omega = (t,s)/(t,t)
        tempr=fasp_blas_darray_dotprod(m,t,t);
        
        if (ABS(tempr)>SMALLREAL) {
            omega=fasp_blas_darray_dotprod(m,s,t)/tempr;
        }
        else {
            if ( prtlvl >= PRINT_SOME ) ITS_DIVZERO;
            omega=0.0;
        }
        
        // delu = alpha pp + omega sp
        fasp_blas_darray_axpby(m,alpha,pp,omega,sp);
        
        // u = u + delu
        fasp_blas_darray_axpy(m,1.0,sp,uval);
        
        // r = s - omega t
        fasp_blas_darray_axpy(m,-omega,t,s);
        fasp_darray_cp(m,s,r);
        
        // beta = (r,rho)/(rp,rho)
        temp2=temp1;
        temp1=fasp_blas_darray_dotprod(m,r,rho);
        
        if (ABS(temp2)>SMALLREAL) {
            beta=(temp1*alpha)/(temp2*omega);
        }
        else {
            ITS_DIVZERO;
            return ERROR_SOLVER_MISC;
        }
        
        // p = p - omega z
        fasp_blas_darray_axpy(m,-omega,z,p);
        
        // p = r + beta p
        fasp_blas_darray_axpby(m,1.0,r,beta,p);
        
        // pp = precond(p)
        if (pc != NULL)
            pc->fct(p,pp,pc->data); /* Apply preconditioner */
        else
            fasp_darray_cp(m,p,pp); /* No preconditioner */
        
        // compute reducation factor of residual ||r||
        absres=fasp_blas_darray_norm2(m,r);
        factor=absres/absres0;
        
        // relative difference
        normd = fasp_blas_darray_norm2(m,sp);
        normu = fasp_blas_darray_norm2(m,uval);
        reldiff = normd/normu;
        
        if ( normd<TOL_s ) {
            ITS_SMALLSP; goto FINISHED;
        }
        
        // relative residual
        switch (stop_type) {
            case STOP_REL_PRECRES:
                if (pc == NULL) fasp_darray_cp(m,r,z);
                else pc->fct(r,z,pc->data);
                tempr=sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                relres=tempr/normr0;
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
        infnormu = fasp_blas_darray_norminf(m, uval);
        if ( infnormu <= sol_inf_tol ) {
            if ( prtlvl > PRINT_MIN ) ITS_ZEROSOL;
            iter = ERROR_SOLVER_SOLSTAG;
            goto FINISHED;
        }
        
        // stagnation check
        if ( (stag<=MaxStag) && (reldiff<maxdiff) ) {
            
            if ( prtlvl >= PRINT_MORE ) {
                ITS_DIFFRES(reldiff,relres);
                ITS_RESTART;
            }
            
            // re-init iteration param
            mf->fct(mf->data, uval, r);
            fasp_blas_darray_axpby(m, 1.0, bval, -1.0, r);
            
            // pp=precond(p)
            fasp_darray_cp(m,r,p);
            if (pc != NULL)
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1=fasp_blas_darray_dotprod(m,r,rho);
            absres=fasp_blas_darray_norm2(m,r);
            
            // relative residual
            switch (stop_type) {
                case STOP_REL_PRECRES:
                    if (pc != NULL)
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    tempr=sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres=tempr/normr0;
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
                    goto FINISHED;
                }
                ++stag;
                ++restart_step;
            }
            
        } // end of stagnation check!
        
        // safe-guard check
        if ( relres < tol ) {
            if ( prtlvl >= PRINT_MORE ) ITS_COMPRES(relres);
            
            // re-init iteration param
            mf->fct(mf->data, uval, r);
            fasp_blas_darray_axpby(m, 1.0, bval, -1.0, r);
            
            // pp=precond(p)
            fasp_darray_cp(m,r,p);
            if (pc != NULL)
                pc->fct(p,pp,pc->data); /* Apply preconditioner */
            else
                fasp_darray_cp(m,p,pp); /* No preconditioner */
            // rho = r* := r
            fasp_darray_cp(m,r,rho);
            temp1=fasp_blas_darray_dotprod(m,r,rho);
            absres=fasp_blas_darray_norm2(m,r);
            
            // relative residual
            switch (stop_type) {
                case STOP_REL_PRECRES:
                    if (pc != NULL)
                        pc->fct(r,z,pc->data);
                    else
                        fasp_darray_cp(m,r,z);
                    tempr=sqrt(ABS(fasp_blas_darray_dotprod(m,r,z)));
                    relres=tempr/normr0;
                    break;
                case STOP_MOD_REL_RES:
                    relres=absres/normu;
                    break;
                default:
                    relres=absres/normr0;
                    break;
            }
            
            if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
            
            // check convergence
            if ( relres < tol ) break;
            
            if ( more_step >= MaxRestartStep ) {
                if ( prtlvl > PRINT_MIN ) ITS_ZEROTOL;
                iter = ERROR_SOLVER_TOLSMALL;
                goto FINISHED;
            }
            else {
                if ( prtlvl > PRINT_NONE ) ITS_RESTART;
            }
            
            ++more_step;
            ++restart_step;
        } // end if safe guard
        
        absres0=absres;
    }
    
FINISHED:  // finish iterative method
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
