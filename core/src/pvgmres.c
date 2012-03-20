/*! \file pvgmres.c
 *  \brief Krylov subspace methods -- Preconditioned Variable-Restarting GMRes.
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
 *      - prINT the result of k-th iteration; 
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
 *  Ref: Iterative methods for sparse linear systems (2nd Edition) 
 *  By Y. Saad, SIAM, 2003.
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
 * \fn INT fasp_solver_dcsr_pvgmres ( dCSRmat *A, dvector *b, dvector *x, const INT MaxIt, 
 *                                    const REAL tol, precond *pc, const SHORT print_level, 
 *                                    const SHORT stop_type, const SHORT restart )
 *
 * \brief Solve "Ax=b" using PGMRES(right preconditioned) iterative method in which the restart
 *        parameter can be adaptively modified during the iteration.
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param x	           Pointer to the dvector of DOFs
 * \param MaxIt        Maximal number of iterations
 * \param tol          Tolerance for stopping
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param print_level  How much information to print out
 * \param stop_type    Stopping criteria type
 * \param restart      Restarting steps
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou 
 * \date   2010/12/14
 *
 * Modified by Chensong Zhang on 2011/12/13
 */ 
INT fasp_solver_dcsr_pvgmres (dCSRmat *A, 
                              dvector *b, 
                              dvector *x, 
                              const INT MaxIt,
                              const REAL tol,
                              precond *pc, 
                              const SHORT print_level, 
                              const SHORT stop_type, 
                              const SHORT restart)
{
	const INT n                 = A->row;  
	const INT min_iter          = 0;

	//--------------------------------------------//
	//   Newly added parameters to monitor when   //
	//   to change the restart parameter          //
	//--------------------------------------------//	
	const REAL cr_max           = 0.99;    // = cos(8^o)  (experimental) 
	const REAL cr_min           = 0.174;   // = cos(80^o) (experimental)
	
    // local variables
	SHORT  converged            = 0; 
	INT    iter                 = 0;
	INT    status               = SUCCESS;	
	INT    restartplus1         = restart + 1;
	INT    i,j,k;
	
	REAL   epsmac               = SMALLREAL; 
	REAL   r_norm, b_norm, den_norm;
	REAL   epsilon, gamma, t, r_norm_0;   
	
	REAL  *c = NULL, *s = NULL, *rs = NULL; 
	REAL  *norms = NULL, *r = NULL, *w = NULL;
	REAL  **p = NULL, **hh = NULL;
	REAL  *work = NULL;
	
	REAL   cr          = 1.0;     // convergence rate
	REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
	INT    d           = 3;       // reduction for the restart parameter 
	INT    restart_max = restart; // upper bound for restart in each restart cycle 
	INT    restart_min = 3;	      // lower bound for restart in each restart cycle (should be small)
	INT    Restart;               // the real restart in some fixed restarted cycle
	
	/* allocate memory */
	work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));	
	p  = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));	
	hh = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
	
	if (print_level > PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
	
	r = work; w = r + n; rs = w + n; c = rs + restartplus1; s = c + restart;	
	for (i = 0; i < restartplus1; i ++) p[i] = s + restart + i*n;
	for (i = 0; i < restartplus1; i ++) hh[i] = p[restart] + n + i*restart;
	
	/* initialization */
	fasp_array_cp(n, b->val, p[0]);
	fasp_blas_dcsr_aAxpy(-1.0, A, x->val, p[0]);
	
	b_norm = fasp_blas_array_norm2(n, b->val);
	r_norm = fasp_blas_array_norm2(n, p[0]);
	r_norm_0 = r_norm;
	
	if ( print_level > PRINT_NONE)
	{
		norms[0] = r_norm;
		if ( print_level >= PRINT_SOME)
		{
            ITS_PUTNORM("right-hand side", b_norm);
            ITS_PUTNORM("residual", r_norm);
		}
	}
	
	if (b_norm > 0.0)  den_norm = b_norm;
	else               den_norm = r_norm;
	
	epsilon = tol*den_norm;
	
	/* outer iteration cycle */
	while (iter < MaxIt)
	{  
		rs[0] = r_norm;
		r_norm_old = r_norm;
		if (r_norm == 0.0)
		{
			fasp_mem_free(work); 
			fasp_mem_free(p); 
			fasp_mem_free(hh);
			fasp_mem_free(norms);
			return iter; 
		}
		
		//-----------------------------------//
		//   adjust the restart parameter    //
		//-----------------------------------//
		
		if (cr > cr_max || iter == 0)
		{
			Restart = restart_max;
		}
		else if (cr < cr_min)
		{
			Restart = Restart;
		}
		else
		{
			if (Restart - d > restart_min)
			{
				Restart -= d;
			}
			else
			{
				Restart = restart_max;
			}
		}
		
		if (r_norm <= epsilon && iter >= min_iter) 
		{
			fasp_array_cp(n, b->val, r);
			fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
			r_norm = fasp_blas_array_norm2(n, r);
			
			if (r_norm <= epsilon)
			{
				break;
			}
			else
			{
				if (print_level >= PRINT_SOME) ITS_FACONV;
			}
		}
		
		t = 1.0 / r_norm;
		for (j = 0; j < n; j ++) p[0][j] *= t;
		
		/* RESTART CYCLE (right-preconditioning) */
		i = 0;
		while (i < Restart && iter < MaxIt)
		{
			i ++;  iter ++;
			
			fasp_array_set(n, r, 0.0);
			
			/* apply the preconditioner */
			if (pc == NULL)
				fasp_array_cp(n, p[i-1], r);
			else
				pc->fct(p[i-1], r, pc->data);          
			
			fasp_blas_dcsr_mxv(A, r, p[i]);
			
			/* modified Gram_Schmidt */
			for (j = 0; j < i; j ++)
			{
				hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
				fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
			}
			t = fasp_blas_array_norm2(n, p[i]);
			hh[i][i-1] = t;	
			if (t != 0.0)
			{
				t = 1.0/t;
				for (j = 0; j < n; j ++) p[i][j] *= t;
			}
			
			for (j = 1; j < i; ++j)
			{
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
		for (k = i-2; k >= 0; k --)
		{
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
		
		if (r_norm  <= epsilon && iter >= min_iter) 
		{
			fasp_array_cp(n, b->val, r);
			fasp_blas_dcsr_aAxpy(-1.0, A, x->val, r);
			r_norm = fasp_blas_array_norm2(n, r);
			
			if (r_norm  <= epsilon)
			{
				converged = 1; break;
			}
			else
			{
				if (print_level >= PRINT_SOME) ITS_FACONV;
				fasp_array_cp(n, r, p[0]); i = 0;
			}
		} /* end of convergence check */
		
		
		/* compute residual vector and continue loop */
		for (j = i; j > 0; j--)
		{
			rs[j-1] = -s[j-1]*rs[j];
			rs[j] = c[j-1]*rs[j];
		}
		
		if (i) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
		
		for (j = i-1 ; j > 0; j --) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
		
		if (i)
		{
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
	
	if (iter>=MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/*!
 * \fn INT fasp_solver_dbsr_pvgmres( dBSRmat *A, dvector *b, dvector *x, const INT MaxIt, 
 *                                   const REAL tol, precond *pc, const SHORT print_level, 
 *                                   const SHORT stop_type, const SHORT restart )
 *
 * \brief Solve "Ax=b" using PGMRES(right preconditioned) iterative method in which the restart
 *        parameter can be adaptively modified during the iteration.
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param x	           Pointer to the dvector of DOFs
 * \param MaxIt        Maximal number of iterations
 * \param tol          Tolerance for stopping
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param print_level  How much information to print out
 * \param stop_type    Stopping criteria type
 * \param restart      Restarting steps
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou 
 * \date   2010/12/21
 */ 
INT fasp_solver_dbsr_pvgmres (dBSRmat *A, 
                              dvector *b, 
                              dvector *x, 
                              const INT MaxIt, 
                              const REAL tol,
                              precond *pc, 
                              const SHORT print_level, 
                              const SHORT stop_type, 
                              const SHORT restart )
{
	const INT n            = A->ROW*A->nb;  
	const INT min_iter     = 0;
	
	//--------------------------------------------//
	//   Newly added parameters to monitor when   //
	//   to change the restart parameter          //
	//--------------------------------------------//	
	const REAL cr_max      = 0.990;   // = cos(8^o)  (experimental) 
	const REAL cr_min      = 0.174;   // = cos(80^o) (experimental)
	
    // local variables
	SHORT  converged       = 0; 
	INT    iter            = 0;
	INT    status          = SUCCESS;	
	INT    restartplus1    = restart + 1;
	INT    i,j,k;
	
	REAL   epsmac          = SMALLREAL; 
	REAL   r_norm, b_norm, den_norm;
	REAL   epsilon, gamma, t, r_norm_0;   
	
	REAL  *c = NULL, *s = NULL, *rs = NULL; 
	REAL  *norms = NULL, *r = NULL, *w = NULL;
	REAL  **p = NULL, **hh = NULL;
	REAL  *work = NULL;
	
	REAL   cr          = 1.0;     // convergence rate
	REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
	INT    d           = 3;       // reduction for the restart parameter 
	INT    restart_max = restart; // upper bound for restart in each restart cycle 
	INT    restart_min = 3;	      // lower bound for restart in each restart cycle (should be small)
	INT    Restart;               // the real restart in some fixed restarted cycle
	
	/* allocate memory */
	work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));	
	p  = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));	
	hh = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
	
	if (print_level > PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
	
	r = work; w = r + n; rs = w + n; c = rs + restartplus1; s = c + restart;	
	for (i = 0; i < restartplus1; ++i) p[i] = s + restart + i*n;
	for (i = 0; i < restartplus1; ++i) hh[i] = p[restart] + n + i*restart;
	
	/* initialization */
	fasp_array_cp(n, b->val, p[0]);
	fasp_blas_dbsr_aAxpy(-1.0, A, x->val, p[0]);
	
	b_norm = fasp_blas_array_norm2(n, b->val);
	r_norm = fasp_blas_array_norm2(n, p[0]);
	r_norm_0 = r_norm;
	
	if ( print_level > PRINT_NONE)
	{
		norms[0] = r_norm;
		if ( print_level >= PRINT_SOME)
		{
            ITS_PUTNORM("right-hand side", b_norm);
            ITS_PUTNORM("residual", r_norm);
		}
	}
	
	if (b_norm > 0.0)  den_norm = b_norm;
	else               den_norm = r_norm;
	
	epsilon = tol*den_norm;
	
	/* outer iteration cycle */
	while (iter < MaxIt)
	{  
		rs[0] = r_norm;
		r_norm_old = r_norm;
		if (r_norm == 0.0)
		{
			fasp_mem_free(work); 
			fasp_mem_free(p); 
			fasp_mem_free(hh);
			fasp_mem_free(norms);
			return iter; 
		}
		
		//-----------------------------------//
		//   adjust the restart parameter    //
		//-----------------------------------//
		
		if (cr > cr_max || iter == 0)
		{
			Restart = restart_max;
		}
		else if (cr < cr_min)
		{
			Restart = Restart;
		}
		else
		{
			if (Restart - d > restart_min)
			{
				Restart -= d;
			}
			else
			{
				Restart = restart_max;
			}
		}
		
		if (r_norm <= epsilon && iter >= min_iter) 
		{
			fasp_array_cp(n, b->val, r);
			fasp_blas_dbsr_aAxpy(-1.0, A, x->val, r);
			r_norm = fasp_blas_array_norm2(n, r);
			
			if (r_norm <= epsilon)
			{
				break;
			}
			else
			{
				if (print_level >= PRINT_SOME) ITS_FACONV;
			}
		}
		
		t = 1.0 / r_norm;
		for (j = 0; j < n; ++j) p[0][j] *= t;
		
		/* RESTART CYCLE (right-preconditioning) */
		i = 0;
		while (i < Restart && iter < MaxIt)
		{
			++i;  ++iter;
			
			fasp_array_set(n, r, 0.0);
			
			/* apply the preconditioner */
			if (pc == NULL)
				fasp_array_cp(n, p[i-1], r);
			else
				pc->fct(p[i-1], r, pc->data);          
			
			fasp_blas_dbsr_mxv(A, r, p[i]);
			
			/* modified Gram_Schmidt */
			for (j = 0; j < i; ++j)
			{
				hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
				fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
			}
			t = fasp_blas_array_norm2(n, p[i]);
			hh[i][i-1] = t;	
			if (t != 0.0)
			{
				t = 1.0/t;
				for (j = 0; j < n; ++j) p[i][j] *= t;
			}
			
			for (j = 1; j < i; j++)
			{
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
			if (r_norm <= epsilon && iter >= min_iter)
			{
				break;
			}         
		} /* end of restart cycle */
		
		/* now compute solution, first solve upper triangular system */
		
		rs[i-1] = rs[i-1] / hh[i-1][i-1];
		for (k = i-2; k >= 0; --k)
		{
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
		
		if (r_norm  <= epsilon && iter >= min_iter) 
		{
			fasp_array_cp(n, b->val, r);
			fasp_blas_dbsr_aAxpy(-1.0, A, x->val, r);
			r_norm = fasp_blas_array_norm2(n, r);
			
			if (r_norm  <= epsilon)
			{
				converged = 1; break;
			}
			else
			{
				if (print_level >= PRINT_SOME) ITS_FACONV;
				fasp_array_cp(n, r, p[0]); i = 0;
			}
		} /* end of convergence check */
		
		/* compute residual vector and continue loop */
		for (j = i; j > 0; j--)
		{
			rs[j-1] = -s[j-1]*rs[j];
			rs[j] = c[j-1]*rs[j];
		}
		
		if (i) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
		
		for (j = i-1 ; j > 0; --j) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
		
		if (i)
		{
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
	
	if (iter>=MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/*!
 * \fn INT fasp_solver_dstr_pvgmres( dSTRmat *A, dvector *b, dvector *x, const INT MaxIt, 
 *                                   const REAL tol, precond *pc, const SHORT print_level, 
 *                                   const SHORT stop_type, const SHORT restart )
 *
 * \brief Solve "Ax=b" using PGMRES(right preconditioned) iterative method in which the restart
 *        parameter can be adaptively modified during the iteration.
 *
 * \param A	           Pointer to the coefficient matrix
 * \param b	           Pointer to the dvector of right hand side
 * \param x	           Pointer to the dvector of DOFs
 * \param MaxIt        Maximal number of iterations
 * \param tol          Tolerance for stopping
 * \param pc           Pointer to the structure of precondition (precond) 
 * \param print_level  How much information to print out
 * \param stop_type    Stopping criteria type
 * \param restart      Restarting steps
 *
 * \return             Number of iterations if converged, error message otherwise
 *
 * \author Zhiyang Zhou 
 * \date   2010/12/14
 */ 
INT fasp_solver_dstr_pvgmres (dSTRmat *A, 
                              dvector *b, 
                              dvector *x, 
                              const INT MaxIt, 
                              const REAL tol, 
                              precond *pc, 
                              const SHORT print_level, 
                              const SHORT stop_type, 
                              const SHORT restart )
{
	const INT n            = A->nc*A->ngrid;  
	const INT min_iter     = 0;
	
	//--------------------------------------------//
	//   Newly added parameters to monitor when   //
	//   to change the restart parameter          //
	//--------------------------------------------//	
	const REAL cr_max      = 0.990;   // = cos(8^o)  (experimental) 
	const REAL cr_min      = 0.174;   // = cos(80^o) (experimental)
	
    // local variables
	SHORT  converged       = 0; 
	INT    iter            = 0;
	INT    status          = SUCCESS;	
	INT    restartplus1    = restart + 1;
	INT    i,j,k;
	
	REAL   epsmac          = SMALLREAL; 
	REAL   r_norm, b_norm, den_norm;
	REAL   epsilon, gamma, t, r_norm_0;   
	
	REAL  *c = NULL, *s = NULL, *rs = NULL; 
	REAL  *norms = NULL, *r = NULL, *w = NULL;
	REAL **p = NULL, **hh = NULL;
	REAL  *work = NULL;
	
	REAL   cr          = 1.0;     // convergence rate
	REAL   r_norm_old  = 0.0;     // save the residual norm of the previous restart cycle
	INT    d           = 3;       // reduction for the restart parameter 
	INT    restart_max = restart; // upper bound for restart in each restart cycle 
	INT    restart_min = 3;	      // lower bound for restart in each restart cycle (should be small)
	INT    Restart;               // the real restart in some fixed restarted cycle
	
	/* allocate memory */
	work = (REAL *)fasp_mem_calloc((restart+4)*(restart+n)+1-n, sizeof(REAL));	
	p  = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *));	
	hh = (REAL **)fasp_mem_calloc(restartplus1, sizeof(REAL *)); 
	
	if (print_level > PRINT_NONE) norms = (REAL *)fasp_mem_calloc(MaxIt+1, sizeof(REAL)); 
	
	r = work; w = r + n; rs = w + n; c = rs + restartplus1; s = c + restart;	
	for (i = 0; i < restartplus1; ++i) p[i] = s + restart + i*n;
	for (i = 0; i < restartplus1; ++i) hh[i] = p[restart] + n + i*restart;
	
	/* initialization */
	fasp_array_cp(n, b->val, p[0]);
	fasp_blas_dstr_aAxpy(-1.0, A, x->val, p[0]);
	
	b_norm = fasp_blas_array_norm2(n, b->val);
	r_norm = fasp_blas_array_norm2(n, p[0]);
	r_norm_0 = r_norm;
	
	if ( print_level > PRINT_NONE )
	{
		norms[0] = r_norm;
		if ( print_level >= PRINT_SOME )
		{
            ITS_PUTNORM("right-hand side", b_norm);
            ITS_PUTNORM("residual", r_norm);
		}
	}
	
	if (b_norm > 0.0)  den_norm = b_norm;
	else               den_norm = r_norm;
	
	epsilon = tol*den_norm;
	
	/* outer iteration cycle */
	while (iter < MaxIt)
	{  
		rs[0] = r_norm;
		r_norm_old = r_norm;
		if (r_norm == 0.0)
		{
			fasp_mem_free(work); 
			fasp_mem_free(p); 
			fasp_mem_free(hh);
			fasp_mem_free(norms);
			return iter; 
		}
		
		//-----------------------------------//
		//   adjust the restart parameter    //
		//-----------------------------------//
		
		if (cr > cr_max || iter == 0)
		{
			Restart = restart_max;
		}
		else if (cr < cr_min)
		{
			Restart = Restart;
		}
		else
		{
			if (Restart - d > restart_min)
			{
				Restart -= d;
			}
			else
			{
				Restart = restart_max;
			}
		}
		
		if (r_norm <= epsilon && iter >= min_iter) 
		{
			fasp_array_cp(n, b->val, r);
			fasp_blas_dstr_aAxpy(-1.0, A, x->val, r);
			r_norm = fasp_blas_array_norm2(n, r);
			
			if (r_norm <= epsilon)
			{
				break;
			}
			else
			{
				if (print_level >= PRINT_SOME) ITS_FACONV;
			}
		}
		
		t = 1.0 / r_norm;
		for (j = 0; j < n; ++j) p[0][j] *= t;
		
		/* RESTART CYCLE (right-preconditioning) */
		i = 0;
		while (i < Restart && iter < MaxIt)
		{
			++i;  ++iter;
			
			fasp_array_set(n, r, 0.0);
			
			/* apply the preconditioner */
			if (pc == NULL)
				fasp_array_cp(n, p[i-1], r);
			else
				pc->fct(p[i-1], r, pc->data);          
			
			//fasp_blas_dbsr_mxv(A, r, p[i]);
			fasp_blas_dstr_aAxpy(1.0, A, r, p[i]); // we need spmxv_str here. --Zhiyang Zhou
			
			/* modified Gram_Schmidt */
			for (j = 0; j < i; ++j)
			{
				hh[j][i-1] = fasp_blas_array_dotprod(n, p[j], p[i]);
				fasp_blas_array_axpy(n, -hh[j][i-1], p[j], p[i]);
			}
			t = fasp_blas_array_norm2(n, p[i]);
			hh[i][i-1] = t;	
			if (t != 0.0)
			{
				t = 1.0/t;
				for (j = 0; j < n; ++j) p[i][j] *= t;
			}
			
			for (j = 1; j < i; ++j)
			{
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
			if (r_norm <= epsilon && iter >= min_iter)
			{
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
		
		if (r_norm  <= epsilon && iter >= min_iter) 
		{
			fasp_array_cp(n, b->val, r);
			fasp_blas_dstr_aAxpy(-1.0, A, x->val, r);
			r_norm = fasp_blas_array_norm2(n, r);
			
			if (r_norm  <= epsilon)
			{
				converged = 1; break;
			}
			else
			{
				if (print_level > PRINT_SOME) ITS_FACONV;
				fasp_array_cp(n, r, p[0]); i = 0;
			}
		} /* end of convergence check */
		
		/* compute residual vector and continue loop */
		for (j = i; j > 0; j--)
		{
			rs[j-1] = -s[j-1]*rs[j];
			rs[j] = c[j-1]*rs[j];
		}
		
		if (i) fasp_blas_array_axpy(n, rs[i]-1.0, p[i], p[i]);
		
		for (j = i-1 ; j > 0; j --) fasp_blas_array_axpy(n, rs[j], p[j], p[i]);
		
		if (i)
		{
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
	
	if (iter>=MaxIt) return ERROR_SOLVER_MAXIT;
	else if (status<0) return status;
	else return iter;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
