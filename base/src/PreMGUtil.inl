/*! \file PreMGUtil.inl
 *
 *  \brief Routines for multigrid coarsest level solver
 *
 *  \note This file contains Level-4 (Pre) functions, which are used in
 *        PreBSR.c, PreCSR.c, PreMGCycle.c, PreMGCycleFull.c, PreMGRecur.c, 
 *        and PreMGRecurAMLI.c
 *
 *  \warning This file is also used in FASP4BLKOIL!!!
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void fasp_coarse_itsolver (dCSRmat *A, dvector *b, dvector *x,
 *                                       const REAL ctol, const SHORT prt_lvl)
 *
 * \brief Iterative on the coarset level
 *
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  ctol      tolerance for the coarsest level
 * \param  prt_lvl   level of output
 *
 * \author Chensong Zhang
 * \date   01/10/2012
 */
static void fasp_coarse_itsolver (dCSRmat *A,
                                  dvector *b,
                                  dvector *x,
                                  const REAL ctol,
                                  const SHORT prt_lvl)
{
    const INT n = A->row;
    const INT maxit = MAX(250,MIN(n*n, 1000)); // Should NOT be less!

    INT status = fasp_solver_dcsr_spcg(A, b, x, NULL, ctol, maxit, 1, 0);

    // If CG fails to converge, use GMRES as a safety net
    if ( status < 0 ) {
        status = fasp_solver_dcsr_spvgmres(A, b, x, NULL, ctol, maxit, 20, 1, 0);
    }

    if ( status < 0 && prt_lvl >= PRINT_MORE ) {
        printf("### WARNING: Coarse level solver failed to converge!\n");
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
