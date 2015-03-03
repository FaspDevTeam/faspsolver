/*! \file mg_util.inl
 *
 *  \brief Routines for algebraic multigrid cycles
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
    const INT csize  = A->row;
    const INT cmaxit = MAX(250,MIN(csize*csize, 1000)); // Should NOT be less!

    INT status = fasp_solver_dcsr_spcg (A, b, x, NULL, ctol, cmaxit, 1, PRINT_NONE);

    // If PCG fails to converge, use PGMRES as another safe net
    if ( status < 0 ) {
        status = fasp_solver_dcsr_spvgmres (A, b, x, NULL, ctol, cmaxit, 20, 1, PRINT_NONE);
    }

    if ( status < 0 && prt_lvl >= PRINT_MORE ) {
        printf("### WARNING: Coarse level solver failed to converge!\n");
    }
}

/**
 * \fn static void fasp_dcsr_presmoothing (const SHORT smoother, dCSRmat *A,
 *                                         dvector *b, dvector *x,
 *                                         const INT nsweeps, const INT istart,
 *                                         const INT iend, const INT istep,
 *                                         const REAL relax, const SHORT ndeg,
 *                                         const SHORT order, INT *ordering)
 *
 * \brief Multigrid presmoothing
 *
 * \param  smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 * \param  ndeg      degree of the polynomial smoother
 * \param  order     order for smoothing sweeps
 * \param  ordering  user defined ordering
 *
 * \author Chensong Zhang
 * \date   01/10/2012
 *
 * Modified by Xiaozhe on 06/04/2012: add ndeg as input
 * Modified by Chensong on 02/16/2013: GS -> SMOOTHER_GS, etc
 */
static void fasp_dcsr_presmoothing (const SHORT smoother,
                                    dCSRmat *A,
                                    dvector *b,
                                    dvector *x,
                                    const INT nsweeps,
                                    const INT istart,
                                    const INT iend,
                                    const INT istep,
                                    const REAL relax,
                                    const SHORT ndeg,
                                    const SHORT order,
                                    INT *ordering)
{
    switch (smoother) {

        case SMOOTHER_GS:
            if (order == NO_ORDER || ordering == NULL)
                fasp_smoother_dcsr_gs(x, istart, iend, istep, A, b, nsweeps);
            else if (order == CF_ORDER)
                fasp_smoother_dcsr_gs_cf(x, A, b, nsweeps, ordering, 1);
            break;

        case SMOOTHER_SGS:
            fasp_smoother_dcsr_sgs(x, A, b, nsweeps);
            break;

        case SMOOTHER_JACOBI:
            fasp_smoother_dcsr_jacobi(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SMOOTHER_L1DIAG:
            fasp_smoother_dcsr_L1diag(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SMOOTHER_POLY:
            fasp_smoother_dcsr_poly(A, b, x, iend+1, ndeg, nsweeps);
            break;

        case SMOOTHER_SOR:
            fasp_smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_SSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_GSOR:
            fasp_smoother_dcsr_gs (x, istart, iend, istep, A, b, nsweeps);
            fasp_smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_SGSOR:
            fasp_smoother_dcsr_gs (x, istart, iend, istep, A, b, nsweeps);
            fasp_smoother_dcsr_gs (x, iend, istart,-istep, A, b, nsweeps);
            fasp_smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_CG:
            fasp_solver_dcsr_pcg(A, b, x, NULL, 1e-3, nsweeps, 1, PRINT_NONE);
            break;

        default:
            printf("### ERROR: Wrong smoother type %d!\n", smoother);
            fasp_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

/**
 * \fn static void fasp_dcsr_postsmoothing (const SHORT smoother, dCSRmat *A,
 *                                          dvector *b, dvector *x,
 *                                          const INT nsweeps, const INT istart,
 *                                          const INT iend, const INT istep,
 *                                          const REAL relax, const SHORT ndeg,
 *                                          const SHORT order, INT *ordering)
 *
 * \brief Multigrid presmoothing
 *
 * \param  smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 * \param  ndeg      degree of the polynomial smoother
 * \param  order     order for smoothing sweeps
 * \param  ordering  user defined ordering
 *
 * \author Chensong Zhang
 * \date   01/10/2012
 *
 * Modified by Xiaozhe Hu on 06/04/2012: add ndeg as input
 * Modified by Chensong on 02/16/2013: GS -> SMOOTHER_GS, etc
 */
static void fasp_dcsr_postsmoothing (const SHORT smoother,
                                     dCSRmat *A,
                                     dvector *b,
                                     dvector *x,
                                     const INT nsweeps,
                                     const INT istart,
                                     const INT iend,
                                     const INT istep,
                                     const REAL relax,
                                     const SHORT ndeg,
                                     const SHORT order,
                                     INT *ordering)
{
    switch (smoother) {

        case SMOOTHER_GS:
            if (order == NO_ORDER || ordering == NULL) {
                fasp_smoother_dcsr_gs(x, iend, istart, istep, A, b, nsweeps);
            }
            else if (order == CF_ORDER)
                fasp_smoother_dcsr_gs_cf(x, A, b, nsweeps, ordering, -1);
            break;

        case SMOOTHER_SGS:
            fasp_smoother_dcsr_sgs(x, A, b, nsweeps);
            break;

        case SMOOTHER_JACOBI:
            fasp_smoother_dcsr_jacobi(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SMOOTHER_L1DIAG:
            fasp_smoother_dcsr_L1diag(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SMOOTHER_POLY:
            fasp_smoother_dcsr_poly(A, b, x, iend+1, ndeg, nsweeps);
            break;

        case SMOOTHER_SOR:
            fasp_smoother_dcsr_sor(x, iend, istart, istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_SSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,  istep, A, b, nsweeps, relax);
            break;

        case SMOOTHER_GSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_gs (x, iend, istart,  istep, A, b, nsweeps);
            break;

        case SMOOTHER_SGSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,  istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_gs (x, istart, iend, -istep, A, b, nsweeps);
            fasp_smoother_dcsr_gs (x, iend, istart,  istep, A, b, nsweeps);
            break;

        case SMOOTHER_CG:
            fasp_solver_dcsr_pcg(A, b, x, NULL, 1e-3, nsweeps, 1, PRINT_NONE);
            break;

        default:
            printf("### ERROR: Wrong smoother type %d!\n", smoother);
            fasp_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
