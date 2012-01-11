/*! \file wrapper.c
 *  \brief Wrappers for accessing functions by advanced users.
 *  
 *  \note Input variables shall not need fasp.h at all!
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_fwrapper_amg_(int *n, int *nnz, int *ia, int *ja, double *a, double *b, double *u, 
 *     						  double *tol, int *maxit, int *ptrlvl)
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG
 *
 * \param n      num of cols of A
 * \param nnz    num of nonzeros of A
 * \param ia     IA of A in CSR format
 * \param ja     JA of A in CSR format
 * \param a      VAL of A in CSR format
 * \param b      rhs vector
 * \param u      solution vector
 * \param tol    tolerance for iterative solvers
 * \param maxit  max num of iterations
 * \param ptrlvl print level for iterative solvers
 *
 * \author Chensong Zhang
 * \date 09/16/2010
 */
void fasp_fwrapper_amg_ (int *n, 
												 int *nnz, 
												 int *ia, 
												 int *ja, 
												 double *a, 
												 double *b, 
												 double *u, 
												 double *tol, 
												 int *maxit, 
												 int *ptrlvl)
{
	dCSRmat    mat;      // coefficient matrix
	dvector    rhs, sol; // right-hand-side, solution
	AMG_param  amgparam; // parameters for AMG
	
	fasp_param_amg_init(&amgparam);
	amgparam.tol         = *tol;
	amgparam.print_level = *ptrlvl;
	amgparam.maxit    = *maxit;
	
	mat.row = *n; mat.col = *n; mat.nnz = *nnz;
	mat.IA = ia;  mat.JA  = ja; mat.val = a;
	
	rhs.row = *n; rhs.val = b;
	sol.row = *n; sol.val = u;
	
	fasp_solver_amg(&mat, &rhs, &sol, &amgparam);
}

/**
 * \fn void fasp_fwrapper_amg_(int *n, int *nnz, int *ia, int *ja, double *a, double *b, double *u, 
 *     						  double *tol, int *maxit, int *ptrlvl)
 * \brief Solve Ax=b by Krylov method preconditioned by classic AMG
 *
 * \param n       num of cols of A
 * \param nnz     num of nonzeros of A
 * \param ia      IA of A in CSR format
 * \param ja      JA of A in CSR format
 * \param a       VAL of A in CSR format
 * \param b       rhs vector
 * \param u       solution vector
 * \param tol     tolerance for iterative solvers
 * \param maxit   max num of iterations
 * \param ptrlvl  print level for iterative solvers
 *
 * \author Chensong Zhang
 * \date 09/16/2010
 */
void fasp_fwrapper_krylov_amg_ (int *n, 
																int *nnz, 
																int *ia, 
																int *ja, 
																double *a, 
																double *b, 
																double *u, 
																double *tol, 
																int *maxit, 
																int *ptrlvl)
{
	dCSRmat    mat;      // coefficient matrix
	dvector    rhs, sol; // right-hand-side, solution	
	AMG_param  amgparam; // parameters for AMG
	itsolver_param  itparam;  // parameters for itsolver
	
	fasp_param_amg_init(&amgparam);	
	fasp_param_solver_init(&itparam);
	itparam.tol         = *tol;
	itparam.print_level = *ptrlvl;
	itparam.maxit       = *maxit;
	
	mat.row = *n; mat.col = *n; mat.nnz = *nnz;
	mat.IA = ia;  mat.JA  = ja; mat.val = a;
	
	rhs.row = *n; rhs.val = b;
	sol.row = *n; sol.val = u;
	
	fasp_solver_dcsr_krylov_amg(&mat, &rhs, &sol, &itparam, &amgparam);
}

/**
 * \fn void fasp_krylov_stokes_(int *nA, int *nnzA, int *ia, int *ja, double *aval, 
 * 												 int *nB, int *nnzB, int *ib, int *jb, double *bval,
 * 												 int *nM, int *nnzM, int *im, int *jm, double *mval,												 
 * 												 double *b, double *u, double *beta,
 * 												 double *tol, int *maxit, int *ptrlvl)
 * \brief Solve [A B;B' O]x=b by Krylov method with block diagonal preconditioner
 *
 * \param nA       num of cols of A
 * \param nnzA     num of nonzeros of A
 * \param ia       IA of A in CSR format
 * \param ja       JA of A in CSR format
 * \param aval     VAL of A in CSR format
 * \param nB       num of cols of B
 * \param nnzB     num of nonzeros of B
 * \param ib       IA of B in CSR format
 * \param jb       JA of B in CSR format
 * \param bval     VAL of B in CSR format
 * \param nM       num of cols of M
 * \param nnzM     num of nonzeros of M
 * \param im       IA of M in CSR format
 * \param jm       JA of M in CSR format
 * \param mval     VAL of M in CSR format 
 * \param nP       num of cols of P
 * \param nnzP     num of nonzeros of P
 * \param ip       IA of P in CSR format
 * \param jp       JA of P in CSR format
 * \param pval     VAL of P in CSR format 
 * \param b        rhs vector
 * \param u        solution vector
 * \param tol      tolerance for iterative solvers
 * \param maxit    max num of iterations
 * \param ptrlvl   print level for iterative solvers
 *
 * \author Chensong Zhang
 * \date 11/16/2010
 */
void fasp_fwrapper_krylov_stokes_ (int *nA, 
																	 int *nnzA, 
																	 int *ia,
																	 int *ja,
																	 double *aval, 
																	 int *nB,
																	 int *nnzB, 
																	 int *ib,
																	 int *jb,
																	 double *bval,
																	 int *nM, 
																	 int *nnzM,
																	 int *im, 
																	 int *jm, 
																	 double *mval,												 
																	 int *nP, 
																	 int *nnzP, 
																	 int *ip, 
																	 int *jp, 
																	 double *pval,												 
																	 double *b, 
																	 double *u, 
																	 double *beta,
																	 double *tol, 
																	 int *maxit, 
																	 int *ptrlvl)
{
	dCSRmat matA, matB, matBt, matM, matP, zero;      
	block_dCSRmat mat; // coefficient matrix
	dvector rhs, sol; // right-hand-side, solution	
	itsolver_param  itparam;  // parameters for itsolver
	precond_Stokes_param psparam; // parameters for Stokes precond
	precond_Stokes_data  psdata; // data for Stokes precond
	
	// initialize itsolver parameters
	fasp_param_solver_init(&itparam);
	itparam.itsolver_type = SOLVER_MinRes;
	itparam.tol         = *tol;
	itparam.print_level = *ptrlvl;
	itparam.maxit       = *maxit;
	
	// initialize precond parameters
	psparam.AMG_type = CLASSIC_AMG;
	psparam.print_level = 1;
	psparam.max_levels = 10;
	
	// initialize matrix	
	matA.row = *nA; matA.col = *nA; matA.nnz = *nnzA;
	matA.IA = ia;  matA.JA  = ja; matA.val = aval;
	
	matB.row = *nB; matB.col = *nB; matB.nnz = *nnzB;
	matB.IA = ib;  matB.JA  = jb; matB.val = bval;
	
	fasp_dcsr_init(&zero); // zero matrix
	
	fasp_dcsr_trans(&matB, &matBt);
	mat.brow = mat.bcol = 2;
	mat.blocks[0] = &matA; 
	mat.blocks[1] = &matB; 
	mat.blocks[2] = &matBt; 
	mat.blocks[3] = &zero; 
	
	matM.row = *nB; matM.col = *nB; matM.nnz = *nnzB;
	matM.IA = ib;  matM.JA  = jb; matM.val = bval;
	
	matP.row = *nP; matP.col = *nP; matP.nnz = *nnzP;
	matP.IA = ip;  matP.JA  = jp; matP.val = pval;
	
	rhs.row = *nA+*nB; rhs.val = b;
	sol.row = *nA+*nB; sol.val = u;
	
	// initialize precond data
	psdata.colA = *nA;
	psdata.colB = *nB;
	psdata.col = *nA+*nB;
	psdata.beta = *beta;
	psdata.M = &matM;
	psdata.P = &matP;
	
	fasp_solver_bdcsr_krylov_stokes(&mat, &rhs, &sol, &itparam, &psparam, &psdata);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
