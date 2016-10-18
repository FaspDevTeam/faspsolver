/*! \file wrapper.c
 *
 *  \brief Wrappers for accessing functions by advanced users
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_fwrapper_amg_ (INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
 *                              REAL *b, REAL *u, REAL *tol, INT *maxit,
 *                              INT *ptrlvl)
 *
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG
 *
 * \param n       Number of cols of A
 * \param nnz     Number of nonzeros of A
 * \param ia      IA of A in CSR format
 * \param ja      JA of A in CSR format
 * \param a       VAL of A in CSR format
 * \param b       RHS vector
 * \param u       Solution vector
 * \param tol     Tolerance for iterative solvers
 * \param maxit   Max number of iterations
 * \param ptrlvl  Print level for iterative solvers
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
 * \fn void fasp_fwrapper_krylov_amg_ (INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
 *                                     REAL *b, REAL *u, REAL *tol, INT *maxit,
 *                                     INT *ptrlvl)
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
 * \param maxit   Max number of iterations
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
    dCSRmat         mat;      // coefficient matrix
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    itsolver_param  itparam;  // parameters for itsolver
    
    fasp_param_amg_init(&amgparam);
    amgparam.AMG_type             = UA_AMG;
    amgparam.print_level          = *ptrlvl;
    amgparam.aggregation_type     = VMB;

    amgparam.coarse_dof           = 100;
    amgparam.presmooth_iter       = 1;
    amgparam.postsmooth_iter      = 1;
    
    amgparam.strong_coupled       = 0.00;
    amgparam.max_aggregation      = 100;
    
    fasp_param_solver_init(&itparam);
    itparam.tol                   = *tol;
    itparam.print_level           = *ptrlvl;
    itparam.maxit                 = *maxit;
    
    
    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA = ia;  mat.JA  = ja; mat.val = a;
    
    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;
    
    fasp_solver_dcsr_krylov_amg(&mat, &rhs, &sol, &itparam, &amgparam);
}

/**
 * \fn INT fasp_wrapper_dbsr_krylov_amg (INT n, INT nnz, INT nb, INT *ia, INT *ja,
 *                                       REAL *a, REAL *b, REAL *u, REAL tol,
 *                                       INT maxit, INT ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned by AMG (dcsr - > dbsr)
 *
 * \param n       Number of cols of A
 * \param nnz     Number of nonzeros of A
 * \param nb      Size of each small block
 * \param ia      IA of A in CSR format
 * \param ja      JA of A in CSR format
 * \param a       VAL of A in CSR format
 * \param b       RHS vector
 * \param u       Solution vector
 * \param tol     Tolerance for iterative solvers
 * \param maxit   Max number of iterations
 * \param ptrlvl  Print level for iterative solvers
 *
 * \return        Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   03/05/2013
 */
INT fasp_wrapper_dbsr_krylov_amg (INT n,
                                  INT nnz,
                                  INT nb,
                                  INT *ia,
                                  INT *ja,
                                  REAL *a,
                                  REAL *b,
                                  REAL *u,
                                  REAL tol,
                                  INT maxit,
                                  INT ptrlvl)
{
    dCSRmat         mat;      // coefficient matrix in CSR format
    dBSRmat         bsrmat;   // coefficient matrix in BSR format
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    itsolver_param  itparam;  // parameters for itsolver
    INT             status = FASP_SUCCESS; // return parameter
    
    // setup AMG parameters
    fasp_param_amg_init(&amgparam);
    amgparam.AMG_type             = UA_AMG;
    amgparam.print_level          = PRINT_NONE;
    
    amgparam.coarse_dof           = 100;
    amgparam.presmooth_iter       = 1;
    amgparam.postsmooth_iter      = 1;
    
    amgparam.strong_coupled       = 0.00;
    amgparam.max_aggregation      = 100;
    
    amgparam.ILU_type             = ILUk;
    amgparam.ILU_levels           = 1;
    amgparam.ILU_lfil             = 1;
    
    // setup Krylov method parameters
    fasp_param_solver_init(&itparam);
    itparam.tol         = tol;
    itparam.print_level = ptrlvl;
    itparam.maxit       = maxit;
    
    itparam.itsolver_type = SOLVER_VFGMRES;
    itparam.restart       = 30;
    
    mat.row = n; mat.col = n; mat.nnz = nnz;
    mat.IA = ia;  mat.JA  = ja; mat.val = a;
    
    // convert CSR to BSR format
    bsrmat = fasp_format_dcsr_dbsr(&mat, nb);
    
    rhs.row = n; rhs.val = b;
    sol.row = n; sol.val = u;
    
    // solve
    status = fasp_solver_dbsr_krylov_amg(&bsrmat, &rhs, &sol, &itparam, &amgparam);
    
    // clean up
    fasp_dbsr_free(&bsrmat);
    
    return status;
}

/**
 * \fn INT fasp_wrapper_dcoo_dbsr_krylov_amg (INT n, INT nnz, INT nb, INT *ia,
 *                                            INT *ja, REAL *a, REAL *b, REAL *u,
 *                                            REAL tol, INT maxit, INT ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned by AMG (dcoo - > dbsr)
 *
 * \param n       Number of cols of A
 * \param nnz     Number of nonzeros of A
 * \param nb      Size of each small block
 * \param ia      IA of A in COO format
 * \param ja      JA of A in COO format
 * \param a       VAL of A in COO format
 * \param b       RHS vector
 * \param u       Solution vector
 * \param tol     Tolerance for iterative solvers
 * \param maxit   Max number of iterations
 * \param ptrlvl  Print level for iterative solvers
 *
 * \return        Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   03/06/2013
 */
INT fasp_wrapper_dcoo_dbsr_krylov_amg (INT n,
                                       INT nnz,
                                       INT nb,
                                       INT *ia,
                                       INT *ja,
                                       REAL *a,
                                       REAL *b,
                                       REAL *u,
                                       REAL tol,
                                       INT maxit,
                                       INT ptrlvl)
{
    
    dCOOmat         coomat;   // coefficient matrix in COO format
    dCSRmat         csrmat;   // coefficient matrix in CSR format
    dBSRmat         bsrmat;   // coefficient matrix in BSR format
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    itsolver_param  itparam;  // parameters for itsolver
    INT             status = FASP_SUCCESS; // return parameter
    
    // setup AMG parameters
    fasp_param_amg_init(&amgparam);
    amgparam.AMG_type             = UA_AMG;
    amgparam.print_level          = ptrlvl;
    
    amgparam.coarse_dof           = 100;
    amgparam.presmooth_iter       = 1;
    amgparam.postsmooth_iter      = 1;
    
    amgparam.strong_coupled       = 0.00;
    amgparam.max_aggregation      = 100;
    
    amgparam.ILU_type             = ILUk;
    amgparam.ILU_levels           = 1;
    amgparam.ILU_lfil             = 1;
    
    // setup Krylov method parameters
    fasp_param_solver_init(&itparam);
    itparam.tol         = tol;
    itparam.print_level = ptrlvl;
    itparam.maxit       = maxit;
    
    itparam.itsolver_type = SOLVER_VFGMRES;
    itparam.restart       = 30;
    
    // COO format
    coomat.row = n; coomat.col = n; coomat.nnz = nnz;
    coomat.rowind = ia; coomat.colind = ja; coomat.val = a;
    
    // convert COO to CSR format
    fasp_format_dcoo_dcsr(&coomat,&csrmat);
    
    // convert CSR to BSR format
    bsrmat = fasp_format_dcsr_dbsr(&csrmat, nb);
    
    // clean up
    fasp_dcsr_free(&csrmat);
    
    rhs.row = n; rhs.val = b;
    sol.row = n; sol.val = u;
    
    // solve
    status = fasp_solver_dbsr_krylov_amg(&bsrmat, &rhs, &sol, &itparam, &amgparam);
    
    // clean up
    fasp_dbsr_free(&bsrmat);
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
