/*! \file mg_util.inl
 *  \brief Utilies for multigrid cycles
 */

#if With_DISOLVE
extern "C" {
    void DIRECT_MUMPS(const INT *n, const INT *nnz, INT *ia, INT *ja, 
                      REAL *a, REAL *b, REAL *x);
}
#endif 

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
 * \date 01/10/2012
 */
static void fasp_coarse_itsolver (dCSRmat *A,
                                  dvector *b,
                                  dvector *x,
                                  const REAL ctol,
                                  const SHORT prt_lvl)
{
    const INT csize  = A->row;
    const INT cmaxit = MAX(500,MIN(csize*csize, 2000)); // coarse level iteration number

   
    INT Status = SUCCESS;
    
    // fasp_smoother_dcsr_sgs(x, A, b, 10);

    INT status = fasp_solver_dcsr_pcg (A, b, x, NULL, ctol, cmaxit, 1, PRINT_NONE);
    
    if (status < 0) { // If PCG does not converge, use BiCGstab as a saft net.
        status = fasp_solver_dcsr_pvgmres (A, b, x, NULL, ctol, cmaxit, 25, 1, PRINT_NONE);
    }
    
    if ( status < 0 && prt_lvl > PRINT_MIN ) {
        printf("### WARNING: coarse level solver does not converge in %d steps!\n", cmaxit);
    }
}

/**
 * \fn static void fasp_dcsr_presmoothing (const SHORT smoother, dCSRmat *A, dvector *b, dvector *x,
 *                                         const INT nsweeps, const INT istart, const INT iend,
 *                                         const INT istep, const REAL relax, const SHORT ndeg, 
 *                                         const SHORT order, INT *ordering)
 *
 * \brief Multigrid presmoothing
 *
 * \param   smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 * \param  ndeg      degree of the polynomial smoother
 * \param order  order for smoothing sweeps
 * \param ordering user defined ordering
 *
 * \author Chensong Zhang
 * \date 01/10/2012
 *
 * \note Modified by Xiaozhe on 06/04/2012: add ndeg as input
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
    //const SHORT ndeg = 6; // for polynomial smoothers
    
    switch (smoother) {

        case GS:
            if (order == NO_ORDER || ordering == NULL) 
                fasp_smoother_dcsr_gs(x, istart, iend, istep, A, b, nsweeps);
            else if (order == CF_ORDER) 
                fasp_smoother_dcsr_gs_cf(x, A, b, nsweeps, ordering, 1);
            break;
            
        case SGS:
            fasp_smoother_dcsr_sgs(x, A, b, nsweeps);
            break;
        
        case JACOBI: 
            fasp_smoother_dcsr_jacobi(x, istart, iend, istep, A, b, nsweeps);
            break;
        
        case L1_DIAG: 
            fasp_smoother_dcsr_L1diag(x, istart, iend, istep, A, b, nsweeps);
            break;
            
        case POLY:
            fasp_smoother_dcsr_poly(A, b, x, iend+1, ndeg, nsweeps); 
            break;
            
        case SOR:
            fasp_smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            break;
        
        case SSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;
        
        case GSOR:
            fasp_smoother_dcsr_gs (x, istart, iend, istep, A, b, nsweeps);
            fasp_smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;
        
        case SGSOR:
            fasp_smoother_dcsr_gs (x, istart, iend, istep, A, b, nsweeps);
            fasp_smoother_dcsr_gs (x, iend, istart,-istep, A, b, nsweeps);
            fasp_smoother_dcsr_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,-istep, A, b, nsweeps, relax);
            break;
        
        case CG:
            fasp_solver_dcsr_pcg(A, b, x, NULL, 1e-3, nsweeps, 1, PRINT_NONE);
            break;
        
        default:
            printf("### ERROR: Wrong smoother type %d!!!\n", smoother); 
            exit(ERROR_INPUT_PAR);
    }
}

/**
 * \fn static void fasp_dcsr_postsmoothing (const SHORT smoother, dCSRmat *A, dvector *b, dvector *x,
 *                                          const INT nsweeps, const INT istart, const INT iend,
 *                                          const INT istep, const REAL relax, const SHORT ndeg,
 *                                          const SHORT order, INT *ordering)
 *
 * \brief Multigrid presmoothing
 *
 * \param   smoother  type of smoother
 * \param  A         pointer to matrix data
 * \param  b         pointer to rhs data
 * \param  x         pointer to sol data
 * \param  nsweeps   number of smoothing sweeps
 * \param  istart    starting index
 * \param  iend      ending index
 * \param  istep     step size
 * \param  relax     relaxation parameter for SOR-type smoothers
 * \param  ndeg      degree of the polynomial smoother 
 * \param order  order for smoothing sweeps
 * \param ordering user defined ordering
 *
 * \author Chensong Zhang
 * \date 01/10/2012
 *
 * \note: modified by Xiaozhe Hu on 06/04/2012: add ndeg as input
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
    //const SHORT ndeg = 6; // for polynomial smoothers

    switch (smoother) {
        
        case GS:
            if (order == NO_ORDER || ordering == NULL)
                fasp_smoother_dcsr_gs(x, iend, istart, istep, A, b, nsweeps);
            else if (order == CF_ORDER)
                fasp_smoother_dcsr_gs_cf(x, A, b, nsweeps, ordering, -1);
            break;
        
        case SGS:
            fasp_smoother_dcsr_sgs(x, A, b, nsweeps);
            break;
            
        case JACOBI:
            fasp_smoother_dcsr_jacobi(x, iend, istart, istep, A, b, nsweeps);
            break;
            
        case L1_DIAG:
            fasp_smoother_dcsr_L1diag(x, iend, istart, istep, A, b, nsweeps);
            break;            
       
        case POLY:
            fasp_smoother_dcsr_poly(A, b, x, iend+1, ndeg, nsweeps); 
            break;
            
        case SOR:
            fasp_smoother_dcsr_sor(x, iend, istart, istep, A, b, nsweeps, relax);
            break;
        
        case SSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,  istep, A, b, nsweeps, relax);
            break;
        
        case GSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_gs (x, iend, istart,  istep, A, b, nsweeps);
            break;
        
        case SGSOR:
            fasp_smoother_dcsr_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_sor(x, iend, istart,  istep, A, b, nsweeps, relax);
            fasp_smoother_dcsr_gs (x, istart, iend, -istep, A, b, nsweeps);
            fasp_smoother_dcsr_gs (x, iend, istart,  istep, A, b, nsweeps);
            break;
        
        case CG:
            fasp_solver_dcsr_pcg(A, b, x, NULL, 1e-3, nsweeps, 1, PRINT_NONE);
            break;
        
        default:
            printf("### ERROR: Wrong smoother type %d!!!\n", smoother); 
            exit(ERROR_INPUT_PAR);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
