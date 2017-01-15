/*! \file pvbcgs.c
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
#include <float.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dcsr_pvbcgs (dCSRmat *A, dvector *b, dvector *u, precond *pc,
 *                                 const REAL tol, const INT MaxIt,
 *                                 const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b, Rewritten from Matlab 2011a
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chunsheng Feng
 * \date   03/04/2016
 *
 */
INT fasp_solver_dcsr_pvbcgs (dCSRmat     *A,
                             dvector     *b,
                             dvector     *u,
                             precond     *pc,
                             const REAL   tol,
                             const INT    MaxIt,
                             const SHORT  stop_type,
                             const SHORT  prtlvl)
{
    const INT    m = b->row;
    
    // local variables
    REAL     n2b,tolb;
    INT      iter=0, stag = 1, moresteps = 1, maxmsteps=1;
    INT      flag, maxstagsteps, half_step=0;
    REAL     absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL     alpha,beta,omega,rho,rho1,rtv,tt;
    REAL     normr,normr_act,normph,normx,imin;
    REAL     norm_sh,norm_xhalf,normrmin,factor;
    
    REAL     *x = u->val, *bval=b->val;
    
    // allocate temp memory (need 10*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(10*m,sizeof(REAL));
    REAL *r=work, *rt=r+m, *p=rt+m, *v=p+m;
    REAL *ph=v+m, *xhalf=ph+m, *s=xhalf+m, *sh=s+m;
    REAL *t = sh+m, *xmin = t+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    n2b = fasp_blas_array_norm2(m,r);
    
    flag = 1;
    fasp_array_cp(m,x,xmin);
    imin = 0;
    
    iter = 0;
    
    tolb = n2b*tol;
    
    fasp_blas_dcsr_aAxpy(-1.0, A, x, r);
    normr = fasp_blas_array_norm2(m,r);
    normr_act = normr;
    
    relres = normr/n2b;
    // if initial residual is small, no need to iterate!
    if (normr <= tolb) {
        flag =0;
        iter =0;
        goto FINISHED;
    }
    
    // output iteration information if needed
    
    print_itinfo(prtlvl,stop_type,iter,relres,n2b,0.0);
    
    // shadow residual rt = r* := r
    fasp_array_cp(m,r,rt);
    normrmin  = normr;
    
    rho = 1.0;
    omega = 1.0;
    stag = 0;
    alpha = 0.0;
    
    moresteps = 0;
    maxmsteps = 10;
    maxstagsteps = 3;
    
    // loop over maxit iterations (unless convergence or failure)
    for(iter=1;iter <= MaxIt;iter++){
        
        rho1 = rho;
        rho  = fasp_blas_array_dotprod(m,rt,r);
        
        if ((rho ==0.0 )|| (ABS(rho) >= DBL_MAX )) {
            flag = 4;
            goto FINISHED;
        }
        
        if (iter==1) {
            fasp_array_cp(m,r,p);
        }
        else  {
            beta = (rho/rho1)*(alpha/omega);
            
            if ((beta == 0)||( ABS(beta) > DBL_MAX )) {
                flag = 4;
                goto FINISHED;
            }
            
            // p = r + beta * (p - omega * v);
            fasp_blas_array_axpy(m,-omega,v,p);        //p=p - omega*v
            fasp_blas_array_axpby(m,1.0, r, beta, p);  //p = 1.0*r +beta*p
        }
        
        // pp = precond(p) ,ph
        if ( pc != NULL )
            pc->fct(p,ph,pc->data); /* Apply preconditioner */
        // if ph all is infinite then exit need add
        else
            fasp_array_cp(m,p,ph); /* No preconditioner */
        
        // v = A*ph
        fasp_blas_dcsr_mxv(A,ph,v);
        rtv = fasp_blas_array_dotprod(m,rt,v);
        
        if (( rtv==0.0 )||( ABS(rtv) > DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        alpha = rho/rtv;
        
        if ( ABS(alpha) > DBL_MAX ){
            flag = 4;
            ITS_DIVZERO;
            goto FINISHED;
        }
        
        normx =  fasp_blas_array_norm2(m,x);
        normph = fasp_blas_array_norm2(m,ph);
        if (ABS(alpha)*normph < DBL_EPSILON*normx )
            stag = stag + 1;
        else
            stag = 0;
        
        //       xhalf = x + alpha * ph;        // form the "half" iterate
        //       s = r - alpha * v;             // residual associated with xhalf
        fasp_blas_array_axpyz(m, alpha, ph, x , xhalf);  //       z= ax + y
        fasp_blas_array_axpyz(m, -alpha, v, r, s);
        normr = fasp_blas_array_norm2(m,s);  //       normr = norm(s);
        normr_act = normr;
        
        // compute reduction factor of residual ||r||
        absres = normr_act;
        factor = absres/absres0;
        print_itinfo(prtlvl,stop_type,iter,normr_act/n2b,absres,factor);
        
        //   check for convergence
        if ((normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            fasp_array_cp(m,bval,s);
            fasp_blas_dcsr_aAxpy(-1.0,A,xhalf,s);
            normr_act = fasp_blas_array_norm2(m,s);
            
            if (normr_act <= tolb){
                // x = xhalf;
                fasp_array_cp(m,xhalf,x);    // x = xhalf;
                flag = 0;
                imin = iter - 0.5;
                half_step++;
                if ( prtlvl >= PRINT_MORE )
                    printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                           flag,stag,imin,half_step);
                goto FINISHED;
            }
            else {
                if ((stag >= maxstagsteps) && (moresteps == 0))  stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps){
                    //     if ~warned
                    flag = 3;
                    fasp_array_cp(m,xhalf,x);
                    goto FINISHED;
                }
            }
        }
        
        if ( stag >= maxstagsteps){
            flag = 3;
            goto FINISHED;
        }
        
        if (normr_act < normrmin )      // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,xhalf,xmin);
            imin = iter - 0.5;
            half_step++;
            if ( prtlvl >= PRINT_MORE )
                printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                       flag,stag,imin,half_step);
        }
        
        // sh = precond(s)
        if ( pc != NULL ){
            pc->fct(s,sh,pc->data); /* Apply preconditioner */
            //if all is finite
        }
        else
            fasp_array_cp(m,s,sh); /* No preconditioner */
        
        // t = A*sh;
        fasp_blas_dcsr_mxv(A,sh,t);
        //tt = t' * t;
        tt = fasp_blas_array_dotprod(m,t,t);
        if ((tt == 0) ||( tt >= DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        //  omega = (t' * s) / tt;
        omega = fasp_blas_array_dotprod(m,s,t)/tt;
        if (ABS(omega) > DBL_MAX )
        {
            flag = 4;
            goto FINISHED;
        }
        
        norm_sh = fasp_blas_array_norm2(m,sh);
        norm_xhalf = fasp_blas_array_norm2(m,xhalf);
        
        if (ABS(omega)*norm_sh < DBL_EPSILON*norm_xhalf )
            stag = stag + 1;
        else
            stag = 0;
        
        fasp_blas_array_axpyz(m, omega,sh,xhalf, x);  //  x = xhalf + omega * sh;
        fasp_blas_array_axpyz(m, -omega, t, s, r);    //  r = s - omega * t;
        normr = fasp_blas_array_norm2(m,r);           //normr = norm(r);
        normr_act = normr;
        
        // % check for convergence
        if ( (normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            // normr_act = norm(r);
            fasp_array_cp(m,bval,r);
            fasp_blas_dcsr_aAxpy(-1.0,A,x,r);
            normr_act = fasp_blas_array_norm2(m,r);
            if (normr_act <= tolb)
            {
                flag = 0;
                goto FINISHED;
            }
            else
            {
                if ((stag >= maxstagsteps) && (moresteps == 0))       stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps)
                {
                    flag = 3;
                    goto FINISHED;
                }
            }
        }
        
        if (normr_act < normrmin)       // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,x,xmin);
            imin = iter;
        }
        
        if (stag >= maxstagsteps)
        {
            flag = 3;
            goto FINISHED;
        }
        
        if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
        
        absres0 = absres;
    }   // for iter = 1 : maxit
    
FINISHED:  // finish the iterative method
    // returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else {
        fasp_array_cp(m, bval,r);
        fasp_blas_dcsr_aAxpy(-1.0,A,xmin,r);
        normr = fasp_blas_array_norm2(m,r);
        
        if ( normr <= normr_act)  {
            fasp_array_cp(m, xmin,x);
            iter = imin;
            relres = normr/n2b;
        }else {
            iter = iter;
            relres = normr_act/n2b;
        }
    }
    
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);

    if ( prtlvl >= PRINT_MORE )
        printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
               flag,stag,imin,half_step);
    
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
 * \fn INT fasp_solver_dbsr_pvbcgs (dBSRmat *A, dvector *b, dvector *u, precond *pc,
 *                                 const REAL tol, const INT MaxIt,
 *                                 const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b, Rewritten from Matlab 2011a
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chunsheng Feng
 * \date   03/04/2016
 *
 */
INT fasp_solver_dbsr_pvbcgs (dBSRmat     *A,
                             dvector     *b,
                             dvector     *u,
                             precond     *pc,
                             const REAL   tol,
                             const INT    MaxIt,
                             const SHORT  stop_type,
                             const SHORT  prtlvl)
{
    const INT    m = b->row;
    
    // local variables
    REAL     n2b,tolb;
    INT      iter=0, stag = 1, moresteps = 1, maxmsteps=1;
    INT      flag, maxstagsteps, half_step=0;
    REAL     absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL     alpha,beta,omega,rho,rho1,rtv,tt;
    REAL     normr,normr_act,normph,normx,imin;
    REAL     norm_sh,norm_xhalf,normrmin,factor;
    
    REAL     *x = u->val, *bval=b->val;
    
    // allocate temp memory (need 10*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(10*m,sizeof(REAL));
    REAL *r=work, *rt=r+m, *p=rt+m, *v=p+m;
    REAL *ph=v+m, *xhalf=ph+m, *s=xhalf+m, *sh=s+m;
    REAL *t = sh+m, *xmin = t+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    n2b = fasp_blas_array_norm2(m,r);
    
    flag = 1;
    fasp_array_cp(m,x,xmin);
    imin = 0;
    
    iter = 0;
    
    tolb = n2b*tol;
    
    fasp_blas_dbsr_aAxpy(-1.0, A, x, r);
    normr = fasp_blas_array_norm2(m,r);
    normr_act = normr;
    
    relres = normr/n2b;
    // if initial residual is small, no need to iterate!
    if (normr <= tolb) {
        flag =0;
        iter =0;
        goto FINISHED;
    }
    
    // output iteration information if needed
    
    print_itinfo(prtlvl,stop_type,iter,relres,n2b,0.0);
    
    // shadow residual rt = r* := r
    fasp_array_cp(m,r,rt);
    normrmin  = normr;
    
    rho = 1.0;
    omega = 1.0;
    stag = 0;
    alpha =0.0;
    
    moresteps = 0;
    maxmsteps = 10;
    maxstagsteps = 3;
    
    // loop over maxit iterations (unless convergence or failure)
    for(iter=1;iter <= MaxIt;iter++){
        
        rho1 = rho;
        rho  = fasp_blas_array_dotprod(m,rt,r);
        
        if ((rho ==0.0 )|| (ABS(rho) >= DBL_MAX )) {
            flag = 4;
            goto FINISHED;
        }
        
        if (iter==1) {
            fasp_array_cp(m,r,p);
        }
        else  {
            beta = (rho/rho1)*(alpha/omega);
            
            if ((beta == 0)||( ABS(beta) > DBL_MAX )) {
                flag = 4;
                goto FINISHED;
            }
            
            // p = r + beta * (p - omega * v);
            fasp_blas_array_axpy(m,-omega,v,p);        //p=p - omega*v
            fasp_blas_array_axpby(m,1.0, r, beta, p);  //p = 1.0*r +beta*p
        }
        
        // pp = precond(p) ,ph
        if ( pc != NULL )
            pc->fct(p,ph,pc->data); /* Apply preconditioner */
        // if ph all is infinite then exit need add
        else
            fasp_array_cp(m,p,ph); /* No preconditioner */
        
        // v = A*ph
        fasp_blas_dbsr_mxv(A,ph,v);
        rtv = fasp_blas_array_dotprod(m,rt,v);
        
        if (( rtv==0.0 )||( ABS(rtv) > DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        alpha = rho/rtv;
        
        if ( ABS(alpha) > DBL_MAX ){
            flag = 4;
            ITS_DIVZERO;
            goto FINISHED;
        }
        
        normx =  fasp_blas_array_norm2(m,x);
        normph = fasp_blas_array_norm2(m,ph);
        if (ABS(alpha)*normph < DBL_EPSILON*normx )
            stag = stag + 1;
        else
            stag = 0;
        
        //       xhalf = x + alpha * ph;        // form the "half" iterate
        //       s = r - alpha * v;             // residual associated with xhalf
        fasp_blas_array_axpyz(m, alpha, ph, x , xhalf);  //       z= ax + y
        fasp_blas_array_axpyz(m, -alpha, v, r, s);
        normr = fasp_blas_array_norm2(m,s);  //       normr = norm(s);
        normr_act = normr;
        
        // compute reduction factor of residual ||r||
        absres = normr_act;
        factor = absres/absres0;
        print_itinfo(prtlvl,stop_type,iter,normr_act/n2b,absres,factor);
        
        //   check for convergence
        if ((normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            fasp_array_cp(m,bval,s);
            fasp_blas_dbsr_aAxpy(-1.0,A,xhalf,s);
            normr_act = fasp_blas_array_norm2(m,s);
            
            if (normr_act <= tolb){
                // x = xhalf;
                fasp_array_cp(m,xhalf,x);    // x = xhalf;
                flag = 0;
                imin = iter - 0.5;
                half_step++;
                if ( prtlvl >= PRINT_MORE )
                    printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                           flag,stag,imin,half_step);
                goto FINISHED;
            }
            else {
                if ((stag >= maxstagsteps) && (moresteps == 0))  stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps){
                    //     if ~warned
                    flag = 3;
                    fasp_array_cp(m,xhalf,x);
                    goto FINISHED;
                }
            }
        }
        
        if ( stag >= maxstagsteps){
            flag = 3;
            goto FINISHED;
        }
        
        if (normr_act < normrmin )      // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,xhalf,xmin);
            imin = iter - 0.5;
            half_step++;
            if ( prtlvl >= PRINT_MORE )
                printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                       flag,stag,imin,half_step);
        }
        
        // sh = precond(s)
        if ( pc != NULL ){
            pc->fct(s,sh,pc->data); /* Apply preconditioner */
            //if all is finite
        }
        else
            fasp_array_cp(m,s,sh); /* No preconditioner */
        
        // t = A*sh;
        fasp_blas_dbsr_mxv(A,sh,t);
        //tt = t' * t;
        tt = fasp_blas_array_dotprod(m,t,t);
        if ((tt == 0) ||( tt >= DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        //  omega = (t' * s) / tt;
        omega = fasp_blas_array_dotprod(m,s,t)/tt;
        if (ABS(omega) > DBL_MAX )
        {
            flag = 4;
            goto FINISHED;
        }
        
        norm_sh = fasp_blas_array_norm2(m,sh);
        norm_xhalf = fasp_blas_array_norm2(m,xhalf);
        
        if (ABS(omega)*norm_sh < DBL_EPSILON*norm_xhalf )
            stag = stag + 1;
        else
            stag = 0;
        
        fasp_blas_array_axpyz(m, omega,sh,xhalf, x);  //  x = xhalf + omega * sh;
        fasp_blas_array_axpyz(m, -omega, t, s, r);    //  r = s - omega * t;
        normr = fasp_blas_array_norm2(m,r);           //normr = norm(r);
        normr_act = normr;
        
        // % check for convergence
        if ( (normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            // normr_act = norm(r);
            fasp_array_cp(m,bval,r);
            fasp_blas_dbsr_aAxpy(-1.0,A,x,r);
            normr_act = fasp_blas_array_norm2(m,r);
            if (normr_act <= tolb)
            {
                flag = 0;
                goto FINISHED;
            }
            else
            {
                if ((stag >= maxstagsteps) && (moresteps == 0))       stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps)
                {
                    flag = 3;
                    goto FINISHED;
                }
            }
        }
        
        if (normr_act < normrmin)       // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,x,xmin);
            imin = iter;
        }
        
        if (stag >= maxstagsteps)
        {
            flag = 3;
            goto FINISHED;
        }
        
        if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
        
        absres0 = absres;
    }   // for iter = 1 : maxit
    
FINISHED:  // finish the iterative method
    // returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else {
        fasp_array_cp(m, bval,r);
        fasp_blas_dbsr_aAxpy(-1.0,A,xmin,r);
        normr = fasp_blas_array_norm2(m,r);
        
        if ( normr <= normr_act)  {
            fasp_array_cp(m, xmin,x);
            iter = imin;
            relres = normr/n2b;
        }else {
            iter = iter;
            relres = normr_act/n2b;
        }
    }
    
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    if ( prtlvl >= PRINT_MORE )
        printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
               flag,stag,imin,half_step);
    
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
 * \fn INT fasp_solver_dblc_pvbcgs (dBLCmat *A, dvector *b, dvector *u, precond *pc,
 *                                  const REAL tol, const INT MaxIt,
 *                                  const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b, Rewritten from Matlab 2011a
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chunsheng Feng
 * \date   03/04/2016
 *
 */
INT fasp_solver_dblc_pvbcgs (dBLCmat     *A,
                             dvector     *b,
                             dvector     *u,
                             precond     *pc,
                             const REAL   tol,
                             const INT    MaxIt,
                             const SHORT  stop_type,
                             const SHORT  prtlvl)
{
    const INT    m = b->row;
    
    // local variables
    REAL     n2b,tolb;
    INT      iter=0, stag = 1, moresteps = 1, maxmsteps=1;
    INT      flag, maxstagsteps, half_step=0;
    REAL     absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL     alpha,beta,omega,rho,rho1,rtv,tt;
    REAL     normr,normr_act,normph,normx,imin;
    REAL     norm_sh,norm_xhalf,normrmin,factor;
    
    REAL     *x = u->val, *bval=b->val;
    
    // allocate temp memory (need 10*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(10*m,sizeof(REAL));
    REAL *r=work, *rt=r+m, *p=rt+m, *v=p+m;
    REAL *ph=v+m, *xhalf=ph+m, *s=xhalf+m, *sh=s+m;
    REAL *t = sh+m, *xmin = t+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    n2b = fasp_blas_array_norm2(m,r);
    
    flag = 1;
    fasp_array_cp(m,x,xmin);
    imin = 0;
    
    iter = 0;
    
    tolb = n2b*tol;
    
    fasp_blas_dblc_aAxpy(-1.0, A, x, r);
    normr = fasp_blas_array_norm2(m,r);
    normr_act = normr;
    
    relres = normr/n2b;
    // if initial residual is small, no need to iterate!
    if (normr <= tolb) {
        flag =0;
        iter =0;
        goto FINISHED;
    }
    
    // output iteration information if needed
    
    print_itinfo(prtlvl,stop_type,iter,relres,n2b,0.0);
    
    // shadow residual rt = r* := r
    fasp_array_cp(m,r,rt);
    normrmin  = normr;
    
    rho = 1.0;
    omega = 1.0;
    stag = 0;
    alpha =0.0;
    
    moresteps = 0;
    maxmsteps = 10;
    maxstagsteps = 3;
    
    // loop over maxit iterations (unless convergence or failure)
    for(iter=1;iter <= MaxIt;iter++){
        
        rho1 = rho;
        rho  = fasp_blas_array_dotprod(m,rt,r);
        
        if ((rho ==0.0 )|| (ABS(rho) >= DBL_MAX )) {
            flag = 4;
            goto FINISHED;
        }
        
        if (iter==1) {
            fasp_array_cp(m,r,p);
        }
        else  {
            beta = (rho/rho1)*(alpha/omega);
            
            if ((beta == 0)||( ABS(beta) > DBL_MAX )) {
                flag = 4;
                goto FINISHED;
            }
            
            // p = r + beta * (p - omega * v);
            fasp_blas_array_axpy(m,-omega,v,p);        //p=p - omega*v
            fasp_blas_array_axpby(m,1.0, r, beta, p);  //p = 1.0*r +beta*p
        }
        
        // pp = precond(p) ,ph
        if ( pc != NULL )
            pc->fct(p,ph,pc->data); /* Apply preconditioner */
        // if ph all is infinite then exit need add
        else
            fasp_array_cp(m,p,ph); /* No preconditioner */
        
        // v = A*ph
        fasp_blas_dblc_mxv(A,ph,v);
        rtv = fasp_blas_array_dotprod(m,rt,v);
        
        if (( rtv==0.0 )||( ABS(rtv) > DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        alpha = rho/rtv;
        
        if ( ABS(alpha) > DBL_MAX ){
            flag = 4;
            ITS_DIVZERO;
            goto FINISHED;
        }
        
        normx =  fasp_blas_array_norm2(m,x);
        normph = fasp_blas_array_norm2(m,ph);
        if (ABS(alpha)*normph < DBL_EPSILON*normx )
            stag = stag + 1;
        else
            stag = 0;
        
        //       xhalf = x + alpha * ph;        // form the "half" iterate
        //       s = r - alpha * v;             // residual associated with xhalf
        fasp_blas_array_axpyz(m, alpha, ph, x , xhalf);  //       z= ax + y
        fasp_blas_array_axpyz(m, -alpha, v, r, s);
        normr = fasp_blas_array_norm2(m,s);  //       normr = norm(s);
        normr_act = normr;
        
        // compute reduction factor of residual ||r||
        absres = normr_act;
        factor = absres/absres0;
        print_itinfo(prtlvl,stop_type,iter,normr_act/n2b,absres,factor);
        
        //   check for convergence
        if ((normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            fasp_array_cp(m,bval,s);
            fasp_blas_dblc_aAxpy(-1.0,A,xhalf,s);
            normr_act = fasp_blas_array_norm2(m,s);
            
            if (normr_act <= tolb){
                // x = xhalf;
                fasp_array_cp(m,xhalf,x);    // x = xhalf;
                flag = 0;
                imin = iter - 0.5;
                half_step++;
                if ( prtlvl >= PRINT_MORE )
                    printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                           flag,stag,imin,half_step);
                goto FINISHED;
            }
            else {
                if ((stag >= maxstagsteps) && (moresteps == 0))  stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps){
                    //     if ~warned
                    flag = 3;
                    fasp_array_cp(m,xhalf,x);
                    goto FINISHED;
                }
            }
        }
        
        if ( stag >= maxstagsteps){
            flag = 3;
            goto FINISHED;
        }
        
        if (normr_act < normrmin )      // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,xhalf,xmin);
            imin = iter - 0.5;
            half_step++;
            if ( prtlvl >= PRINT_MORE )
                printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                       flag,stag,imin,half_step);
        }
        
        // sh = precond(s)
        if ( pc != NULL ){
            pc->fct(s,sh,pc->data); /* Apply preconditioner */
            //if all is finite
        }
        else
            fasp_array_cp(m,s,sh); /* No preconditioner */
        
        // t = A*sh;
        fasp_blas_dblc_mxv(A,sh,t);
        //tt = t' * t;
        tt = fasp_blas_array_dotprod(m,t,t);
        if ((tt == 0) ||( tt >= DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        //  omega = (t' * s) / tt;
        omega = fasp_blas_array_dotprod(m,s,t)/tt;
        if (ABS(omega) > DBL_MAX )
        {
            flag = 4;
            goto FINISHED;
        }
        
        norm_sh = fasp_blas_array_norm2(m,sh);
        norm_xhalf = fasp_blas_array_norm2(m,xhalf);
        
        if (ABS(omega)*norm_sh < DBL_EPSILON*norm_xhalf )
            stag = stag + 1;
        else
            stag = 0;
        
        fasp_blas_array_axpyz(m, omega,sh,xhalf, x);  //  x = xhalf + omega * sh;
        fasp_blas_array_axpyz(m, -omega, t, s, r);    //  r = s - omega * t;
        normr = fasp_blas_array_norm2(m,r);           //normr = norm(r);
        normr_act = normr;
        
        //  check for convergence
        if ( (normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            // normr_act = norm(r);
            fasp_array_cp(m,bval,r);
            fasp_blas_dblc_aAxpy(-1.0,A,x,r);
            normr_act = fasp_blas_array_norm2(m,r);
            if (normr_act <= tolb)
            {
                flag = 0;
                goto FINISHED;
            }
            else
            {
                if ((stag >= maxstagsteps) && (moresteps == 0))       stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps)
                {
                    flag = 3;
                    goto FINISHED;
                }
            }
        }
        
        if (normr_act < normrmin)       // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,x,xmin);
            imin = iter;
        }
        
        if (stag >= maxstagsteps)
        {
            flag = 3;
            goto FINISHED;
        }
        
        if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
        
        absres0 = absres;
    }   // for iter = 1 : maxit
    
FINISHED:  // finish the iterative method
    // returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else {
        fasp_array_cp(m, bval,r);
        fasp_blas_dblc_aAxpy(-1.0,A,xmin,r);
        normr = fasp_blas_array_norm2(m,r);
        
        if ( normr <= normr_act)  {
            fasp_array_cp(m, xmin,x);
            iter = imin;
            relres = normr/n2b;
        }else {
            iter = iter;
            relres = normr_act/n2b;
        }
    }
    
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    if ( prtlvl >= PRINT_MORE )
        printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
               flag,stag,imin,half_step);
    
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
 * \fn INT fasp_solver_dstr_pvbcgs (dSTRmat *A, dvector *b, dvector *u, precond *pc,
 *                                  const REAL tol, const INT MaxIt,
 *                                  const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b, Rewritten from Matlab 2011a
 *
 * \param A            Pointer to the coefficient matrix
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chunsheng Feng
 * \date   03/04/2016
 *
 */
INT fasp_solver_dstr_pvbcgs (dSTRmat     *A,
                             dvector     *b,
                             dvector     *u,
                             precond     *pc,
                             const REAL   tol,
                             const INT    MaxIt,
                             const SHORT  stop_type,
                             const SHORT  prtlvl)
{
    const INT    m = b->row;
    
    // local variables
    REAL     n2b,tolb;
    INT      iter=0, stag = 1, moresteps = 1, maxmsteps=1;
    INT      flag, maxstagsteps, half_step=0;
    REAL     absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL     alpha,beta,omega,rho,rho1,rtv,tt;
    REAL     normr,normr_act,normph,normx,imin;
    REAL     norm_sh,norm_xhalf,normrmin,factor;
    
    REAL     *x = u->val, *bval=b->val;
    
    // allocate temp memory (need 10*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(10*m,sizeof(REAL));
    REAL *r=work, *rt=r+m, *p=rt+m, *v=p+m;
    REAL *ph=v+m, *xhalf=ph+m, *s=xhalf+m, *sh=s+m;
    REAL *t = sh+m, *xmin = t+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    n2b = fasp_blas_array_norm2(m,r);
    
    flag = 1;
    fasp_array_cp(m,x,xmin);
    imin = 0;
    
    iter = 0;
    
    tolb = n2b*tol;
    
    fasp_blas_dstr_aAxpy(-1.0, A, x, r);
    normr = fasp_blas_array_norm2(m,r);
    normr_act = normr;
    
    relres = normr/n2b;
    // if initial residual is small, no need to iterate!
    if (normr <= tolb) {
        flag =0;
        iter =0;
        goto FINISHED;
    }
    
    // output iteration information if needed
    
    print_itinfo(prtlvl,stop_type,iter,relres,n2b,0.0);
    
    // shadow residual rt = r* := r
    fasp_array_cp(m,r,rt);
    normrmin  = normr;
    
    rho = 1.0;
    omega = 1.0;
    stag = 0;
    alpha =0.0;
    
    moresteps = 0;
    maxmsteps = 10;
    maxstagsteps = 3;
    
    // loop over maxit iterations (unless convergence or failure)
    for(iter=1;iter <= MaxIt;iter++){
        
        rho1 = rho;
        rho  = fasp_blas_array_dotprod(m,rt,r);
        
        if ((rho ==0.0 )|| (ABS(rho) >= DBL_MAX )) {
            flag = 4;
            goto FINISHED;
        }
        
        if (iter==1) {
            fasp_array_cp(m,r,p);
        }
        else  {
            beta = (rho/rho1)*(alpha/omega);
            
            if ((beta == 0)||( ABS(beta) > DBL_MAX )) {
                flag = 4;
                goto FINISHED;
            }
            
            // p = r + beta * (p - omega * v);
            fasp_blas_array_axpy(m,-omega,v,p);        //p=p - omega*v
            fasp_blas_array_axpby(m,1.0, r, beta, p);  //p = 1.0*r +beta*p
        }
        
        // pp = precond(p) ,ph
        if ( pc != NULL )
            pc->fct(p,ph,pc->data); /* Apply preconditioner */
        // if ph all is infinite then exit need add
        else
            fasp_array_cp(m,p,ph); /* No preconditioner */
        
        // v = A*ph
        fasp_blas_dstr_mxv(A,ph,v);
        rtv = fasp_blas_array_dotprod(m,rt,v);
        
        if (( rtv==0.0 )||( ABS(rtv) > DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        alpha = rho/rtv;
        
        if ( ABS(alpha) > DBL_MAX ){
            flag = 4;
            ITS_DIVZERO;
            goto FINISHED;
        }
        
        normx =  fasp_blas_array_norm2(m,x);
        normph = fasp_blas_array_norm2(m,ph);
        if (ABS(alpha)*normph < DBL_EPSILON*normx )
            stag = stag + 1;
        else
            stag = 0;
        
        //       xhalf = x + alpha * ph;        // form the "half" iterate
        //       s = r - alpha * v;             // residual associated with xhalf
        fasp_blas_array_axpyz(m, alpha, ph, x , xhalf);  //       z= ax + y
        fasp_blas_array_axpyz(m, -alpha, v, r, s);
        normr = fasp_blas_array_norm2(m,s);  //       normr = norm(s);
        normr_act = normr;
        
        // compute reduction factor of residual ||r||
        absres = normr_act;
        factor = absres/absres0;
        print_itinfo(prtlvl,stop_type,iter,normr_act/n2b,absres,factor);
        
        //   check for convergence
        if ((normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            fasp_array_cp(m,bval,s);
            fasp_blas_dstr_aAxpy(-1.0,A,xhalf,s);
            normr_act = fasp_blas_array_norm2(m,s);
            
            if (normr_act <= tolb){
                // x = xhalf;
                fasp_array_cp(m,xhalf,x);    // x = xhalf;
                flag = 0;
                imin = iter - 0.5;
                half_step++;
                if ( prtlvl >= PRINT_MORE )
                    printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                           flag,stag,imin,half_step);
                goto FINISHED;
            }
            else {
                if ((stag >= maxstagsteps) && (moresteps == 0))  stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps){
                    //     if ~warned
                    flag = 3;
                    fasp_array_cp(m,xhalf,x);
                    goto FINISHED;
                }
            } 
        }
        
        if ( stag >= maxstagsteps){
            flag = 3;
            goto FINISHED;
        }
        
        if (normr_act < normrmin )      // update minimal norm quantities
        {  
            normrmin = normr_act;
            fasp_array_cp(m,xhalf,xmin); 
            imin = iter - 0.5;
            half_step++;
            if ( prtlvl >= PRINT_MORE )
                printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                       flag,stag,imin,half_step);
        }
        
        // sh = precond(s)
        if ( pc != NULL ){
            pc->fct(s,sh,pc->data); /* Apply preconditioner */
            //if all is finite
        }
        else
            fasp_array_cp(m,s,sh); /* No preconditioner */
        
        // t = A*sh; 
        fasp_blas_dstr_mxv(A,sh,t);
        //tt = t' * t;
        tt = fasp_blas_array_dotprod(m,t,t);  
        if ((tt == 0) ||( tt >= DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        //  omega = (t' * s) / tt;
        omega = fasp_blas_array_dotprod(m,s,t)/tt;
        if (ABS(omega) > DBL_MAX ) 
        {
            flag = 4;
            goto FINISHED;
        }
        
        norm_sh = fasp_blas_array_norm2(m,sh);
        norm_xhalf = fasp_blas_array_norm2(m,xhalf);
        
        if (ABS(omega)*norm_sh < DBL_EPSILON*norm_xhalf )
            stag = stag + 1;
        else
            stag = 0;
        
        fasp_blas_array_axpyz(m, omega,sh,xhalf, x);  //  x = xhalf + omega * sh;
        fasp_blas_array_axpyz(m, -omega, t, s, r);    //  r = s - omega * t;
        normr = fasp_blas_array_norm2(m,r);           //normr = norm(r);
        normr_act = normr;
        
        // % check for convergence        
        if ( (normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            // normr_act = norm(r);
            fasp_array_cp(m,bval,r);
            fasp_blas_dstr_aAxpy(-1.0,A,x,r);
            normr_act = fasp_blas_array_norm2(m,r);
            if (normr_act <= tolb) 
            {
                flag = 0;
                goto FINISHED;
            }
            else     
            {
                if ((stag >= maxstagsteps) && (moresteps == 0))       stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps) 
                {
                    flag = 3;
                    goto FINISHED;
                }
            }        
        }
        
        if (normr_act < normrmin)       // update minimal norm quantities
        {  
            normrmin = normr_act;
            fasp_array_cp(m,x,xmin);
            imin = iter;
        }        
        
        if (stag >= maxstagsteps)
        {
            flag = 3;
            goto FINISHED;
        }
        
        if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
        
        absres0 = absres;
    }   // for iter = 1 : maxit
    
    
FINISHED:  // finish the iterative method
    // returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else {	  
        fasp_array_cp(m, bval,r);
        fasp_blas_dstr_aAxpy(-1.0,A,xmin,r);
        normr = fasp_blas_array_norm2(m,r);
        
        if ( normr <= normr_act)  {
            fasp_array_cp(m, xmin,x);
            iter = imin;
            relres = normr/n2b;
        }else {
            iter = iter;
            relres = normr_act/n2b;
        } 
    }
    
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    if ( prtlvl >= PRINT_MORE )
        printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
               flag,stag,imin,half_step);
    
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
 * \fn INT fasp_solver_pvbcgs (mxv_matfree *mf, dvector *b, dvector *u, precond *pc,
 *                             const REAL tol, const INT MaxIt,
 *                             const SHORT stop_type, const SHORT prtlvl)
 *
 * \brief Preconditioned BiCGstab method for solving Au=b, Rewritten from Matlab 2011a
 *
 * \param mf           Pointer to mxv_matfree: the spmv operation
 * \param b            Pointer to the dvector of right hand side
 * \param u            Pointer to the dvector of DOFs
 * \param pc           Pointer to precond: the structure of precondition
 * \param tol          Tolerance for stopping
 * \param MaxIt        Maximal number of iterations
 * \param stop_type    Stopping criteria type
 * \param prtlvl       How much information to print out
 *
 * \return             Iteration number if converges; ERROR otherwise.
 *
 * \author Chunsheng Feng
 * \date   03/04/2016
 */
INT fasp_solver_pvbcgs (mxv_matfree *mf,
                        dvector     *b,
                        dvector     *u,
                        precond     *pc,
                        const REAL   tol,
                        const INT    MaxIt,
                        const SHORT  stop_type,
                        const SHORT  prtlvl)
{
    const INT    m = b->row;
    
    // local variables
    REAL     n2b,tolb;
    INT      iter=0, stag = 1, moresteps = 1, maxmsteps=1;
    INT      flag, maxstagsteps, half_step=0;
    REAL     absres0 = BIGREAL, absres = BIGREAL, relres = BIGREAL;
    REAL     alpha,beta,omega,rho,rho1,rtv,tt;
    REAL     normr,normr_act,normph,normx,imin;
    REAL     norm_sh,norm_xhalf,normrmin,factor;
    
    REAL     *x = u->val, *bval=b->val;
    
    // allocate temp memory (need 10*m REAL)
    REAL *work=(REAL *)fasp_mem_calloc(10*m,sizeof(REAL));
    REAL *r=work, *rt=r+m, *p=rt+m, *v=p+m;
    REAL *ph=v+m, *xhalf=ph+m, *s=xhalf+m, *sh=s+m;
    REAL *t = sh+m, *xmin = t+m;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: maxit = %d, tol = %.4le\n", MaxIt, tol);
#endif
    
    // r = b-A*u
    fasp_array_cp(m,bval,r);
    n2b = fasp_blas_array_norm2(m,r);
    
    flag = 1;
    fasp_array_cp(m,x,xmin);
    imin = 0;
    
    iter = 0;
    
    tolb = n2b*tol;
    
    // r = b-A*x
    mf->fct(mf->data, x, r);
    fasp_blas_array_axpby(m, 1.0, bval, -1.0, r);
    normr = fasp_blas_array_norm2(m,r);
    normr_act = normr;
    
    relres = normr/n2b;
    // if initial residual is small, no need to iterate!
    if (normr <= tolb) {
        flag =0;
        iter =0;
        goto FINISHED;
    }
    
    // output iteration information if needed
    
    print_itinfo(prtlvl,stop_type,iter,relres,n2b,0.0);
    
    // shadow residual rt = r* := r
    fasp_array_cp(m,r,rt);
    normrmin  = normr;
    
    rho = 1.0;
    omega = 1.0;
    stag = 0;
    alpha =0.0;
    
    moresteps = 0;
    maxmsteps = 10;
    maxstagsteps = 3;
    
    // loop over maxit iterations (unless convergence or failure)
    for(iter=1;iter <= MaxIt;iter++){
        
        rho1 = rho;
        rho  = fasp_blas_array_dotprod(m,rt,r);
        
        if ((rho ==0.0 )|| (ABS(rho) >= DBL_MAX )) {
            flag = 4;
            goto FINISHED;
        }
        
        if (iter==1) {
            fasp_array_cp(m,r,p);
        }
        else  {
            beta = (rho/rho1)*(alpha/omega);
            
            if ((beta == 0)||( ABS(beta) > DBL_MAX )) {
                flag = 4;
                goto FINISHED;
            }
            
            // p = r + beta * (p - omega * v);
            fasp_blas_array_axpy(m,-omega,v,p);        //p=p - omega*v
            fasp_blas_array_axpby(m,1.0, r, beta, p);  //p = 1.0*r +beta*p
        }
        
        // pp = precond(p) ,ph
        if ( pc != NULL )
            pc->fct(p,ph,pc->data); /* Apply preconditioner */
        // if ph all is infinite then exit need add
        else
            fasp_array_cp(m,p,ph); /* No preconditioner */
        
        // v = A*ph
        mf->fct(mf->data, ph, v);
        rtv = fasp_blas_array_dotprod(m,rt,v);
        
        if (( rtv==0.0 )||( ABS(rtv) > DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        alpha = rho/rtv;
        
        if ( ABS(alpha) > DBL_MAX ){
            flag = 4;
            ITS_DIVZERO;
            goto FINISHED;
        }
        
        normx =  fasp_blas_array_norm2(m,x);
        normph = fasp_blas_array_norm2(m,ph);
        if (ABS(alpha)*normph < DBL_EPSILON*normx )
            stag = stag + 1;
        else
            stag = 0;
        
        //       xhalf = x + alpha * ph;        // form the "half" iterate
        //       s = r - alpha * v;             // residual associated with xhalf
        fasp_blas_array_axpyz(m, alpha, ph, x , xhalf);  //       z= ax + y
        fasp_blas_array_axpyz(m, -alpha, v, r, s);
        normr = fasp_blas_array_norm2(m,s);  //       normr = norm(s);
        normr_act = normr;
        
        // compute reduction factor of residual ||r||
        absres = normr_act;
        factor = absres/absres0;
        print_itinfo(prtlvl,stop_type,iter,normr_act/n2b,absres,factor);
        
        //   check for convergence
        if ((normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            // s = b-A*xhalf
            mf->fct(mf->data, xhalf, s);
            fasp_blas_array_axpby(m, 1.0, bval, -1.0, s);
            normr_act = fasp_blas_array_norm2(m,s);
            
            if (normr_act <= tolb){
                // x = xhalf;
                fasp_array_cp(m,xhalf,x);    // x = xhalf;
                flag = 0;
                imin = iter - 0.5;
                half_step++;
                if ( prtlvl >= PRINT_MORE )
                    printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                           flag,stag,imin,half_step);
                goto FINISHED;
            }
            else {
                if ((stag >= maxstagsteps) && (moresteps == 0))  stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps){
                    //     if ~warned
                    flag = 3;
                    fasp_array_cp(m,xhalf,x);
                    goto FINISHED;
                }
            }
        }
        
        if ( stag >= maxstagsteps){
            flag = 3;
            goto FINISHED;
        }
        
        if (normr_act < normrmin )      // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,xhalf,xmin);
            imin = iter - 0.5;
            half_step++;
            if ( prtlvl >= PRINT_MORE )
                printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
                       flag,stag,imin,half_step);
        }
        
        // sh = precond(s)
        if ( pc != NULL ){
            pc->fct(s,sh,pc->data); /* Apply preconditioner */
            //if all is finite
        }
        else
            fasp_array_cp(m,s,sh); /* No preconditioner */
        
        // t = A*sh;
        mf->fct(mf->data, sh, t);
        //tt = t' * t;
        tt = fasp_blas_array_dotprod(m,t,t);
        if ((tt == 0) ||( tt >= DBL_MAX )){
            flag = 4;
            goto FINISHED;
        }
        
        //  omega = (t' * s) / tt;
        omega = fasp_blas_array_dotprod(m,s,t)/tt;
        if (ABS(omega) > DBL_MAX )
        {
            flag = 4;
            goto FINISHED;
        }
        
        norm_sh = fasp_blas_array_norm2(m,sh);
        norm_xhalf = fasp_blas_array_norm2(m,xhalf);
        
        if (ABS(omega)*norm_sh < DBL_EPSILON*norm_xhalf )
            stag = stag + 1;
        else
            stag = 0;
        
        fasp_blas_array_axpyz(m, omega,sh,xhalf, x);  //  x = xhalf + omega * sh;
        fasp_blas_array_axpyz(m, -omega, t, s, r);    //  r = s - omega * t;
        normr = fasp_blas_array_norm2(m,r);           //normr = norm(r);
        normr_act = normr;
        
        // % check for convergence
        if ( (normr <= tolb)||(stag >= maxstagsteps)||moresteps)
        {
            // normr_act = norm(r);
            // r = b-A*x
            mf->fct(mf->data, x, r);
            fasp_blas_array_axpby(m, 1.0, bval, -1.0, r);
            normr_act = fasp_blas_array_norm2(m,r);
            if (normr_act <= tolb)
            {
                flag = 0;
                goto FINISHED;
            }
            else
            {
                if ((stag >= maxstagsteps) && (moresteps == 0))       stag = 0;
                
                moresteps = moresteps + 1;
                if (moresteps >= maxmsteps)
                {
                    flag = 3;
                    goto FINISHED;
                }
            }
        }
        
        if (normr_act < normrmin)       // update minimal norm quantities
        {
            normrmin = normr_act;
            fasp_array_cp(m,x,xmin);
            imin = iter;
        }
        
        if (stag >= maxstagsteps)
        {
            flag = 3;
            goto FINISHED;
        }
        
        if ( prtlvl >= PRINT_MORE ) ITS_REALRES(relres);
        
        absres0 = absres;
    }   // for iter = 1 : maxit
    
FINISHED:  // finish the iterative method
    // returned solution is first with minimal residual
    if (flag == 0)
        relres = normr_act / n2b;
    else {
        // r = b-A*xmin
        mf->fct(mf->data, xmin, r);
        fasp_blas_array_axpby(m, 1.0, bval, -1.0, r);
        normr = fasp_blas_array_norm2(m,r);
        
        if ( normr <= normr_act)  {
            fasp_array_cp(m, xmin,x);
            iter = imin;
            relres = normr/n2b;
        }else {
            iter = iter;
            relres = normr_act/n2b;
        }
    }
    
    if ( prtlvl > PRINT_NONE ) ITS_FINAL(iter,MaxIt,relres);
    
    if ( prtlvl >= PRINT_MORE )
        printf("Flag = %d Stag = %d Itermin = %.1f Half_step = %d\n",
               flag,stag,imin,half_step);
    
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
