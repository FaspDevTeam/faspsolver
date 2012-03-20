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
 * \fn void void fasp_fwrapper_amg_ (INT *n, INT *nnz, INT *ia, INT *ja, REAL *a, REAL *b, 
 *                                   REAL *u, REAL *tol, INT *maxit, INT *ptrlvl)
 *
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG
 *
 * \param n      Number of cols of A
 * \param nnz    Number of nonzeros of A
 * \param ia     IA of A in CSR format
 * \param ja     JA of A in CSR format
 * \param a      VAL of A in CSR format
 * \param b      RHS vector
 * \param u      Solution vector
 * \param tol    Tolerance for iterative solvers
 * \param maxit  Max num of iterations
 * \param ptrlvl Print level for iterative solvers
 *
 * \author Chensong Zhang
 * \date   09/16/2010
 */
void fasp_fwrapper_amg_ (INT *n, 
                         INT *nnz, 
                         INT *ia, 
                         INT *ja, 
                         REAL *a, 
                         REAL *b, 
                         REAL *u, 
                         REAL *tol, 
                         INT *maxit, 
                         INT *ptrlvl)
{
	dCSRmat    mat;      // coefficient matrix
	dvector    rhs, sol; // right-hand-side, solution
	AMG_param  amgparam; // parameters for AMG
	
	fasp_param_amg_init(&amgparam);
	amgparam.tol         = *tol;
	amgparam.print_level = *ptrlvl;
	amgparam.maxit    = *maxit;
	
	mat.row = *n; mat.col = *n; mat.nnz = *nnz;
	mat.IA  = ia; mat.JA  = ja; mat.val = a;
	
	rhs.row = *n; rhs.val = b;
	sol.row = *n; sol.val = u;
	
	fasp_solver_amg(&mat, &rhs, &sol, &amgparam);
}

/**
 * \fn void fasp_fwrapper_krylov_amg_ (INT *n, INT *nnz, INT *ia, INT *ja, REAL *a, REAL *b, REAL *u, 
 *                                     REAL *tol, INT *maxit, INT *ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned by classic AMG
 *
 * \param n       Number of cols of A
 * \param nnz     Number of nonzeros of A
 * \param ia      IA of A in CSR format
 * \param ja      JA of A in CSR format
 * \param a       VAL of A in CSR format
 * \param b       RHS vector
 * \param u       Solution vector
 * \param tol     Tolerance for iterative solvers
 * \param maxit   Max num of iterations
 * \param ptrlvl  Print level for iterative solvers
 *
 * \author Chensong Zhang
 * \date   09/16/2010
 */
void fasp_fwrapper_krylov_amg_ (INT *n, 
                                INT *nnz, 
                                INT *ia, 
                                INT *ja, 
                                REAL *a, 
                                REAL *b, 
                                REAL *u, 
                                REAL *tol, 
                                INT *maxit, 
                                INT *ptrlvl)
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


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
