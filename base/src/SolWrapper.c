/*! \file  SolWrapper.c
 *
 *  \brief Wrappers for accessing functions by advanced users
 *
 *  \note  This file contains Level-5 (Sol) functions. It requires:
 *         AuxParam.c, BlaFormat.c, BlaSparseBSR.c, BlaSparseCSR.c, SolAMG.c,
 *         SolBSR.c, and SolCSR.c
 *
 *  \note  IMPORTANT: The wrappers DO NOT change the original matrix data. Users
 *         should shift the matrix indices in order to make the IA and JA to start
 *         from 0 instead of 1.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2018 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_fwrapper_dcsr_amg_ (INT *n, INT *nnz, INT *ia, INT *ja, REAL *a,
 *                                   REAL *b, REAL *u, REAL *tol, INT *maxit,
 *                                   INT *ptrlvl)
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
void fasp_fwrapper_dcsr_amg_ (INT  *n,
                              INT  *nnz,
                              INT  *ia,
                              INT  *ja,
                              REAL *a,
                              REAL *b,
                              REAL *u,
                              REAL *tol,
                              INT  *maxit,
                              INT  *ptrlvl)
{
    dCSRmat    mat;      // coefficient matrix
    dvector    rhs, sol; // right-hand-side, solution
    AMG_param  amgparam; // parameters for AMG
    
    // setup AMG parameters
    fasp_param_amg_init(&amgparam);
    
    amgparam.tol         = *tol;
    amgparam.print_level = *ptrlvl;
    amgparam.maxit       = *maxit;
    
    // set up coefficient matrix
    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA  = ia; mat.JA  = ja; mat.val = a;
    
    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;
    
    fasp_solver_amg(&mat, &rhs, &sol, &amgparam);
}

/**
 * \fn void fasp_fwrapper_dcsr_krylov_ilu_ (INT *n, INT *nnz, INT *ia, INT *ja,
 *                                          REAL *a, REAL *b, REAL *u, REAL *tol,
 *                                          INT *maxit, INT *ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned by ILUk
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
 * \date   03/24/2018
 */
void fasp_fwrapper_dcsr_krylov_ilu_ (INT  *n,
                                     INT  *nnz,
                                     INT  *ia,
                                     INT  *ja,
                                     REAL *a,
                                     REAL *b,
                                     REAL *u,
                                     REAL *tol,
                                     INT  *maxit,
                                     INT  *ptrlvl)
{
    dCSRmat    mat;      // coefficient matrix
    dvector    rhs, sol; // right-hand-side, solution
    ILU_param  iluparam; // parameters for ILU
    ITS_param  itsparam; // parameters for itsolver

    // setup ILU parameters
    fasp_param_ilu_init(&iluparam);

    iluparam.print_level   = *ptrlvl;

    // setup Krylov method parameters
    fasp_param_solver_init(&itsparam);

    itsparam.itsolver_type = SOLVER_VFGMRES;
    itsparam.tol           = *tol;
    itsparam.maxit         = *maxit;
    itsparam.print_level   = *ptrlvl;

    // set up coefficient matrix
    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA  = ia;  mat.JA = ja; mat.val = a;

    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;

    fasp_solver_dcsr_krylov_ilu(&mat, &rhs, &sol, &itsparam, &iluparam);
}

/**
 * \fn void fasp_fwrapper_dcsr_krylov_amg_ (INT *n, INT *nnz, INT *ia, INT *ja,
 *                                          REAL *a, REAL *b, REAL *u, REAL *tol,
 *                                          INT *maxit, INT *ptrlvl)
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
void fasp_fwrapper_dcsr_krylov_amg_ (INT  *n,
                                     INT  *nnz,
                                     INT  *ia,
                                     INT  *ja,
                                     REAL *a,
                                     REAL *b,
                                     REAL *u,
                                     REAL *tol,
                                     INT  *maxit,
                                     INT  *ptrlvl)
{
    dCSRmat    mat;      // coefficient matrix
    dvector    rhs, sol; // right-hand-side, solution
    AMG_param  amgparam; // parameters for AMG
    ITS_param  itsparam; // parameters for itsolver
    
    // setup AMG parameters
    fasp_param_amg_init(&amgparam);
    
    amgparam.AMG_type             = UA_AMG;
    amgparam.aggregation_type     = VMB;
    amgparam.coarse_dof           = 100;
    amgparam.max_aggregation      = 100;
    amgparam.presmooth_iter       = 1;
    amgparam.postsmooth_iter      = 1;
    amgparam.strong_coupled       = 0.00;
    amgparam.print_level          = *ptrlvl;

    // setup Krylov method parameters
    fasp_param_solver_init(&itsparam);
    
    itsparam.tol                  = *tol;
    itsparam.maxit                = *maxit;
    itsparam.print_level          = *ptrlvl;

    // set up coefficient matrix
    mat.row = *n; mat.col = *n; mat.nnz = *nnz;
    mat.IA  = ia;  mat.JA = ja; mat.val = a;
    
    rhs.row = *n; rhs.val = b;
    sol.row = *n; sol.val = u;
    
    fasp_solver_dcsr_krylov_amg(&mat, &rhs, &sol, &itsparam, &amgparam);
}

/**
 * \fn void fasp_wrapper_dbsr_krylov_ilu (INT n, INT nnz, INT nb, INT *ia, INT *ja,
 *                                        REAL *a, REAL *b, REAL *u, REAL tol,
 *                                        INT maxit, INT ptrlvl)
 *
 * \brief Solve Ax=b by Krylov method preconditioned by block ILU in BSR format
 *
 * \param n       Number of cols of A
 * \param nnz     Number of nonzeros of A
 * \param nb      Size of each small block
 * \param ia      IA of A in BSR format
 * \param ja      JA of A in BSR format
 * \param a       VAL of A in BSR format
 * \param b       RHS vector
 * \param u       Solution vector
 * \param tol     Tolerance for iterative solvers
 * \param maxit   Max number of iterations
 * \param ptrlvl  Print level for iterative solvers
 *
 * \return        Iteration number if converges; ERROR otherwise.
 *
 * \author Chensong Zhang
 * \date   03/25/2018
 */
void fasp_wrapper_dbsr_krylov_ilu (INT  *n,
                                   INT  *nnz,
                                   INT  *nb,
                                   INT  *ia,
                                   INT  *ja,
                                   REAL *a,
                                   REAL *b,
                                   REAL *u,
                                   REAL *tol,
                                   INT  *maxit,
                                   INT  *ptrlvl)
{
    dBSRmat    mat;      // coefficient matrix in BSR format
    dvector    rhs, sol; // right-hand-side, solution

    ILU_param  iluparam; // parameters for ILU
    ITS_param  itsparam; // parameters for itsolver

    // setup ILU parameters
    fasp_param_ilu_init(&iluparam);

    iluparam.print_level   = *ptrlvl;

    // setup Krylov method parameters
    fasp_param_solver_init(&itsparam);

    itsparam.itsolver_type = SOLVER_VFGMRES;
    itsparam.tol           = *tol;
    itsparam.maxit         = *maxit;
    itsparam.print_level   = *ptrlvl;

    // set up coefficient matrix
    mat.ROW = *n; mat.COL = *n; mat.NNZ = *nnz; mat.nb = *nb;
    mat.IA  = ia; mat.JA  = ja; mat.val = a;

    rhs.row = *n * *nb; rhs.val = b;
    sol.row = *n * *nb; sol.val = u;

    // solve
    fasp_solver_dbsr_krylov_ilu(&mat, &rhs, &sol, &itsparam, &iluparam);
}

#if 0 // Remove later --Chensong

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
INT fasp_wrapper_dbsr_krylov_amg (INT   n,
                                  INT   nnz,
                                  INT   nb,
                                  INT  *ia,
                                  INT  *ja,
                                  REAL *a,
                                  REAL *b,
                                  REAL *u,
                                  REAL  tol,
                                  INT   maxit,
                                  INT   ptrlvl)
{
    dCSRmat    mat;      // coefficient matrix in CSR format
    dBSRmat    bsrmat;   // coefficient matrix in BSR format
    dvector    rhs, sol; // right-hand-side, solution
    AMG_param  amgparam; // parameters for AMG
    ITS_param  itsparam; // parameters for itsolver
    INT        status = FASP_SUCCESS; // return parameter
    
    // setup AMG parameters
    fasp_param_amg_init(&amgparam);
    
    amgparam.AMG_type        = UA_AMG;
    amgparam.print_level     = PRINT_NONE;
    amgparam.coarse_dof      = 100;
    amgparam.max_aggregation = 100;
    amgparam.presmooth_iter  = 1;
    amgparam.postsmooth_iter = 1;
    amgparam.strong_coupled  = 0.00;
    amgparam.ILU_type        = ILUk;
    amgparam.ILU_levels      = 1;
    amgparam.ILU_lfil        = 1;
    
    // setup Krylov method parameters
    fasp_param_solver_init(&itsparam);
    
    itsparam.tol              = tol;
    itsparam.print_level      = ptrlvl;
    itsparam.maxit            = maxit;
    itsparam.itsolver_type    = SOLVER_VFGMRES;
    itsparam.restart          = 30;
    
    // set up coefficient matrix
    mat.row = n; mat.col = n;  mat.nnz = nnz;
    mat.IA = ia; mat.JA  = ja; mat.val = a;
    
    // convert CSR to BSR format
    bsrmat = fasp_format_dcsr_dbsr(&mat, nb);
    
    rhs.row = n; rhs.val = b;
    sol.row = n; sol.val = u;
    
    // solve
    status = fasp_solver_dbsr_krylov_amg(&bsrmat, &rhs, &sol, &itsparam, &amgparam);
    
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
INT fasp_wrapper_dcoo_dbsr_krylov_amg (INT   n,
                                       INT   nnz,
                                       INT   nb,
                                       INT  *ia,
                                       INT  *ja,
                                       REAL *a,
                                       REAL *b,
                                       REAL *u,
                                       REAL  tol,
                                       INT   maxit,
                                       INT   ptrlvl)
{
    
    dCOOmat         coomat;   // coefficient matrix in COO format
    dCSRmat         csrmat;   // coefficient matrix in CSR format
    dBSRmat         bsrmat;   // coefficient matrix in BSR format
    dvector         rhs, sol; // right-hand-side, solution
    AMG_param       amgparam; // parameters for AMG
    ITS_param       itsparam; // parameters for itsolver
    INT             status = FASP_SUCCESS; // return parameter
    
    // setup AMG parameters
    fasp_param_amg_init(&amgparam);
    
    amgparam.AMG_type        = UA_AMG;
    amgparam.print_level     = ptrlvl;
    amgparam.coarse_dof      = 100;
    amgparam.max_aggregation = 100;
    amgparam.presmooth_iter  = 1;
    amgparam.postsmooth_iter = 1;
    amgparam.strong_coupled  = 0.00;
    amgparam.ILU_type        = ILUk;
    amgparam.ILU_levels      = 1;
    amgparam.ILU_lfil        = 1;
    
    // setup Krylov method parameters
    fasp_param_solver_init(&itsparam);
    
    itsparam.tol              = tol;
    itsparam.print_level      = ptrlvl;
    itsparam.maxit            = maxit;
    itsparam.itsolver_type    = SOLVER_VFGMRES;
    itsparam.restart          = 30;
    
    // set up coefficient matrix
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
    status = fasp_solver_dbsr_krylov_amg(&bsrmat, &rhs, &sol, &itsparam, &amgparam);
    
    // clean up
    fasp_dbsr_free(&bsrmat);
    
    return status;
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
