/*! \file  KryUtil.inl
 *
 *  \brief Routines for iterative solvers
 *
 *  \note  This file contains Level-3 (Kry) functions, which are used in:
 *         KryPbcgs.c, KryPcg.c, KryPgcg.c, KryPgcr.c, KryPgmres.c
 *         KryPminres.c, KryPvbcgs.c, KryPvfgmres.c, KryPvgmres.c
 *         KrySPbcgs.c, KrySPcg.c, KrySPgmres.c, KrySPminres.c,
 *         KrySPvgmres.c, PreMGSolve.c, SolBLC.c, SolBSR.c, SolCSR.c
 *         SolMatFree.c, and SolSTR.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  \warning This file is also used in FASP4BLKOIL!!!
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

//! Warning for residual false convergence
#define ITS_FACONV  printf("### WARNING: False convergence!\n")

//! Warning for solution close to zero
#define ITS_ZEROSOL printf("### WARNING: Iteration stopped due to the solution is almost zero! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for iteration restarted
#define ITS_RESTART printf("### WARNING: Iteration restarted due to stagnation! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for stagged iteration
#define ITS_STAGGED printf("### WARNING: Iteration stopped due to staggnation! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for tolerance practically close to zero
#define ITS_ZEROTOL printf("### WARNING: The tolerence might be too small! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for divided by zero
#define ITS_DIVZERO printf("### WARNING: Divided by zero! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for actual relative residual
#define ITS_REALRES(relres) printf("### WARNING: The actual relative residual = %e!\n",(relres))

//! Warning for computed relative residual
#define ITS_COMPRES(relres) printf("### WARNING: The computed relative residual = %e!\n",(relres))

//! Warning for too small sp 
#define ITS_SMALLSP printf("### WARNING: sp is too small! %s : %d\n", __FUNCTION__, __LINE__)

//! Warning for restore previous iteration 
#define ITS_RESTORE(iter) printf("### WARNING: Discard current solution. Restore iteration:%d!\n",(iter));

//! Output relative difference and residual
#define ITS_DIFFRES(reldiff,relres) printf("||u-u'|| = %e and the comp. rel. res. = %e.\n",(reldiff),(relres));

//! Output L2 norm of some variable
#define ITS_PUTNORM(name,value) printf("L2 norm of %s = %e.\n",(name),(value));

/**
 * \fn inline static void ITS_CHECK (const INT MaxIt, const REAL tol)
 * \brief Safeguard checks to prevent unexpected error for iterative solvers
 *
 * \param MaxIt   Maximal number of iterations
 * \param tol     Tolerance for convergence check
 *
 * \author Chensong Zhang
 * \date   05/16/2012
 */
inline static void ITS_CHECK (const INT MaxIt, const REAL tol) 
{    
    if ( tol < SMALLREAL ) {
        printf("### WARNING: Convergence tolerance for iterative solver is too small!\n");
    }
    if ( MaxIt <= 0 ) {
        printf("### WARNING: Max number of iterations should be a POSITIVE integer!\n");
    }
}

/**
 * \fn inline static void ITS_FINAL (const INT iter, const INT MaxIt, const REAL relres) 
 * \brief Print out final status of an iterative method
 *
 * \param iter    Number of iterations
 * \param MaxIt   Maximal number of iterations
 * \param relres  Relative residual 
 *
 * \author Chensong Zhang
 * \date   01/11/2012
 */
inline static void ITS_FINAL (const INT iter, const INT MaxIt, const REAL relres) 
{
    if ( iter > MaxIt ) {
        printf("### WARNING: Max iter %d reached with rel. resid. %e.\n", MaxIt, relres);
    }
    else if ( iter >= 0 ) {
        printf("Number of iterations = %d with relative residual %e.\n", iter, relres);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
