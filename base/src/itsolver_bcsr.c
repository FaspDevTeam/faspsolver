/*! \file itsolver_bcsr.c
 *
 *  \brief Iterative solvers for block_dCSRmat matrices
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_bdcsr_itsolver (block_dCSRmat *A, dvector *b, dvector *x,
 *                                     precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A        Pointer to the coeff matrix in block_dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Iteration number if converges; ERROR otherwise.
 *
 *
 * \author Chensong Zhang
 * \date   11/25/2010
 */
INT fasp_solver_bdcsr_itsolver (block_dCSRmat *A,
                                dvector *b,
                                dvector *x,
                                precond *pc,
                                itsolver_param *itparam)
{
    const SHORT prtlvl = itparam->print_level;
    const SHORT itsolver_type = itparam->itsolver_type;
    const SHORT stop_type = itparam->stop_type;
    const SHORT restart = itparam->restart;
    const INT   MaxIt = itparam->maxit;
    const REAL  tol = itparam->tol;
    
    REAL  solver_start, solver_end, solver_duration;
    INT   iter = ERROR_SOLVER_TYPE;
    
    fasp_gettime(&solver_start);

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);
#endif

    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );
    
    switch (itsolver_type) {
            
        case SOLVER_BiCGstab:
            if ( prtlvl > PRINT_NONE ) printf("\nCalling BiCGstab solver (Block CSR) ...\n");
            iter=fasp_solver_bdcsr_pbcgs(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            break;
            
        case SOLVER_MinRes:
            if ( prtlvl > PRINT_NONE ) printf("\nCalling MinRes solver (Block CSR) ...\n");
            iter=fasp_solver_bdcsr_pminres(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            break;
            
        case SOLVER_GMRES:
            if ( prtlvl > PRINT_NONE ) printf("\nCalling GMRES solver (Block CSR) ...\n");
            iter=fasp_solver_bdcsr_pgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;
            
        case SOLVER_VGMRES:
            if ( prtlvl > PRINT_NONE ) printf("Calling vGMRES solver (Block CSR) ...\n");
            iter=fasp_solver_bdcsr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;
            
        case SOLVER_VFGMRES:
            if ( prtlvl > PRINT_NONE ) printf("Calling FGMRES solver (Block CSR) ...\n");
            iter=fasp_solver_bdcsr_pvfgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, prtlvl);
            break;
            
        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
            
    }
    
    if ( (prtlvl >= PRINT_MIN) && (iter >= 0) ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return iter;
}

/**
 * \fn INT fasp_solver_bdcsr_krylov (block_dCSRmat *A, dvector *b, dvector *x,
 *                                   itsolver_param *itparam)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   07/18/2010
 */
INT fasp_solver_bdcsr_krylov (block_dCSRmat *A,
                              dvector *b,
                              dvector *x,
                              itsolver_param *itparam)
{
    const SHORT prtlvl = itparam->print_level;
    
    INT status = FASP_SUCCESS;
    REAL solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    // solver part
    fasp_gettime(&solver_start);
    
    status = fasp_solver_bdcsr_itsolver(A,b,x,NULL,itparam);
    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( prtlvl >= PRINT_MIN )
        print_cputime("Krylov method totally", solver_duration);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn INT fasp_solver_bdcsr_krylov_block_3 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   07/10/2014
 *
 * \note only works for 3by3 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT fasp_solver_bdcsr_krylov_block_3 (block_dCSRmat *A,
                                      dvector *b,
                                      dvector *x,
                                      itsolver_param *itparam,
                                      AMG_param *amgparam,
                                      dCSRmat *A_diag)
{
    const SHORT prtlvl = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    
    INT status = FASP_SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;
    
    const SHORT max_levels = amgparam->max_levels;
    INT m, n, nnz, i;
    
    AMG_data **mgl = NULL;
    
#if WITH_UMFPACK
    void **LU_diag = (void **)fasp_mem_calloc(3, sizeof(void *));
#endif
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    /* setup preconditioner */
    fasp_gettime(&setup_start);
    
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
        
#if WITH_UMFPACK
        // Need to sort the diagonal blocks for UMFPACK format
        dCSRmat A_tran;
        
        for (i=0; i<3; i++){
            
            A_tran = fasp_dcsr_create(A_diag[i].row, A_diag[i].col, A_diag[i].nnz);
            fasp_dcsr_transz(&A_diag[i], NULL, &A_tran);
            fasp_dcsr_cp(&A_tran, &A_diag[i]);
            
            printf("Factorization for %d-th diagnol: \n", i);
            LU_diag[i] = fasp_umfpack_factorize(&A_diag[i], prtlvl);
            
        }
        
        fasp_dcsr_free(&A_tran);
#endif
    
    }
    
    /* diagonal blocks are solved by AMG */
    else if ( precond_type > 30 && precond_type < 40 ) {
     
        mgl = (AMG_data **)fasp_mem_calloc(3, sizeof(AMG_data *));
        
        for (i=0; i<3; i++){
            
            mgl[i] = fasp_amg_data_create(max_levels);
            m = A_diag[i].row; n = A_diag[i].col; nnz = A_diag[i].nnz;
            mgl[i][0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(&A_diag[i],&mgl[i][0].A);
            mgl[i][0].b=fasp_dvec_create(n); mgl[i][0].x=fasp_dvec_create(n);
            
            switch (amgparam->AMG_type) {
                case SA_AMG: // Smoothed Aggregation AMG
                    status = fasp_amg_setup_sa(mgl[i], amgparam); break;
                case UA_AMG: // Unsmoothed Aggregation AMG
                    status = fasp_amg_setup_ua(mgl[i], amgparam); break;
                default: // Classical AMG
                    status = fasp_amg_setup_rs(mgl[i], amgparam); break;
            }

        }
        
    }
    
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
    precond_block_data precdata;
    precdata.Abcsr = A;
    precdata.A_diag = A_diag;
    precdata.r = fasp_dvec_create(b->row);
    
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
#if WITH_UMFPACK
        precdata.LU_diag = LU_diag;
#endif
    }
    /* diagonal blocks are solved by AMG */
    else if ( precond_type > 30 && precond_type < 40 ) {
        precdata.amgparam = amgparam;
        precdata.mgl = mgl;
    }
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
    precond prec; prec.data = &precdata;
    
    switch (precond_type)
    {
        case 21:
            prec.fct = fasp_precond_block_diag_3;
            break;
            
        case 22:
            prec.fct = fasp_precond_block_lower_3;
            break;
            
        case 23:
            prec.fct = fasp_precond_block_upper_3;
            break;
            
        case 24:
            prec.fct = fasp_precond_block_SGS_3;
            break;
            
        case 31:
            prec.fct = fasp_precond_block_diag_3_amg;
            break;
            
        case 32:
            prec.fct = fasp_precond_block_lower_3_amg;
            break;
            
        case 33:
            prec.fct = fasp_precond_block_upper_3_amg;
            break;
        
        case 34:
            prec.fct = fasp_precond_block_SGS_3_amg;
            break;
            
        default:
            fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
            break;
    }
    
    if ( prtlvl >= PRINT_MIN ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }
    
    
    // solver part
    fasp_gettime(&solver_start);
    
    status=fasp_solver_bdcsr_itsolver(A,b,x, &prec,itparam);
    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( prtlvl >= PRINT_MIN )
        print_cputime("Krylov method totally", solver_duration);
    
    // clean
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
#if WITH_UMFPACK
        for (i=0; i<3; i++) fasp_umfpack_free_numeric(LU_diag[i]);
#endif
    }
    /* diagonal blocks are solved by AMG */
    else if ( precond_type > 30 && precond_type < 40 ) {
        for (i=0; i<3; i++) fasp_amg_data_free(mgl[i], amgparam);
        if (mgl) fasp_mem_free(mgl);
    }
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn INT fasp_solver_bdcsr_krylov_block_4 (block_dCSRmat *A, dvector *b, dvector *x,
 *                                           itsolver_param *itparam,
 *                                           AMG_param *amgparam, dCSRmat *A_diag)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A         Pointer to the coeff matrix in block_dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG solvers
 * \param A_diag    Digonal blocks of A
 *
 * \return          Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   07/06/2014
 *
 * \note only works for 4 by 4 block dCSRmat problems!! -- Xiaozhe Hu
 */
INT fasp_solver_bdcsr_krylov_block_4 (block_dCSRmat *A,
                                      dvector *b,
                                      dvector *x,
                                      itsolver_param *itparam,
                                      AMG_param *amgparam,
                                      dCSRmat *A_diag)
{
    const SHORT prtlvl = itparam->print_level;
    const SHORT precond_type = itparam->precond_type;
    
    INT status = FASP_SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solver_start, solver_end, solver_duration;

#if WITH_UMFPACK
    void **LU_diag = (void **)fasp_mem_calloc(4, sizeof(void *));
    INT i;
#endif

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    /* setup preconditioner */
    fasp_gettime(&setup_start);
    
    /* diagonal blocks are solved exactly */
    if ( precond_type > 20 && precond_type < 30 ) {
    
#if WITH_UMFPACK
    // Need to sort the matrices local_A for UMFPACK format
    dCSRmat A_tran;
    
    for (i=0; i<4; i++){
        
        A_tran = fasp_dcsr_create(A_diag[i].row, A_diag[i].col, A_diag[i].nnz);
        fasp_dcsr_transz(&A_diag[i], NULL, &A_tran);
        fasp_dcsr_cp(&A_tran, &A_diag[i]);
        
        printf("Factorization for %d-th diagnol: \n", i);
        LU_diag[i] = fasp_umfpack_factorize(&A_diag[i], prtlvl);
        
    }
    
    fasp_dcsr_free(&A_tran);
#endif
        
    }
    else {
        fasp_chkerr(ERROR_SOLVER_PRECTYPE, __FUNCTION__);
    }
    
    precond_block_data precdata;
    
    precdata.Abcsr = A;
    precdata.A_diag = A_diag;
#if WITH_UMFPACK
    precdata.LU_diag = LU_diag;
#endif
    precdata.r = fasp_dvec_create(b->row);
    
    precond prec; prec.data = &precdata;
    
    switch (precond_type)
    {
        case 21:
            prec.fct = fasp_precond_block_diag_4;
            break;
            
        case 22:
            prec.fct = fasp_precond_block_lower_4;
            break;
    }
    
    if ( prtlvl >= PRINT_MIN ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }

    // solver part
    fasp_gettime(&solver_start);
    
    status=fasp_solver_bdcsr_itsolver(A,b,x, &prec,itparam);
    
    fasp_gettime(&solver_end);
    
    solver_duration = solver_end - solver_start;
    
    if ( prtlvl >= PRINT_MIN )
        print_cputime("Krylov method totally", solver_duration);
    
    // clean
#if WITH_UMFPACK
    for (i=0; i<4; i++) fasp_umfpack_free_numeric(LU_diag[i]);
#endif
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/**
 * \fn INT fasp_solver_bdcsr_krylov_sweeping (block_dCSRmat *A, dvector *b,
 *                                            dvector *x, itsolver_param *itparam,
 *                                            INT NumLayers, block_dCSRmat *Ai,
 *                                            dCSRmat *local_A, ivector *local_index)
 *
 * \brief Solve Ax = b by standard Krylov methods
 *
 * \param A             Pointer to the coeff matrix in block_dCSRmat format
 * \param b             Pointer to the right hand side in dvector format
 * \param x             Pointer to the approx solution in dvector format
 * \param itparam       Pointer to parameters for iterative solvers
 * \param NumLayers     Number of layers used for sweeping preconditioner
 * \param Ai            Pointer to the coeff matrix for the preconditioner in block_dCSRmat format
 * \param local_A       Pointer to the local coeff matrices in the dCSRmat format
 * \param local_index   Pointer to the local index in ivector format
 *
 * \return              Iteration number if converges; ERROR otherwise.
 *
 * \author Xiaozhe Hu
 * \date   05/01/2014
 */
INT fasp_solver_bdcsr_krylov_sweeping (block_dCSRmat *A,
                                       dvector *b,
                                       dvector *x,
                                       itsolver_param *itparam,
                                       INT NumLayers,
                                       block_dCSRmat *Ai,
                                       dCSRmat *local_A,
                                       ivector *local_index)
{
    const SHORT prtlvl = itparam->print_level;
    
    INT status = FASP_SUCCESS;
    REAL setup_start, setup_end, setup_duration;
    REAL solve_start, solve_end, solve_duration;
    
    /* setup preconditioner */
    fasp_gettime(&setup_start);
    
    void **local_LU = NULL;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
#if WITH_UMFPACK
    // Need to sort the matrices local_A for UMFPACK format
    INT l;
    dCSRmat A_tran;
    local_LU = (void **)fasp_mem_calloc(NumLayers, sizeof(void *));
    
    for ( l=0; l<NumLayers; l++ ) {
        
        A_tran = fasp_dcsr_create(local_A[l].row, local_A[l].col, local_A[l].nnz);
        fasp_dcsr_transz(&local_A[l], NULL, &A_tran);
        fasp_dcsr_cp(&A_tran, &local_A[l]);
        
        printf("Factorization for layer %d: \n", l);
        local_LU[l] = fasp_umfpack_factorize(&local_A[l], prtlvl);
        
    }
    
    fasp_dcsr_free(&A_tran);
#endif
    
    precond_sweeping_data precdata;
    precdata.NumLayers = NumLayers;
    precdata.A = A;
    precdata.Ai = Ai;
    precdata.local_A = local_A;
    precdata.local_LU = local_LU;
    precdata.local_index = local_index;
    precdata.r = fasp_dvec_create(b->row);
    precdata.w = (REAL *)fasp_mem_calloc(10*b->row,sizeof(REAL));
    
    precond prec; prec.data = &precdata;
    prec.fct = fasp_precond_sweeping;
    
    if ( prtlvl >= PRINT_MIN ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_cputime("Setup totally", setup_duration);
    }
    
    /* solver part */
    fasp_gettime(&solve_start);
    
    status = fasp_solver_bdcsr_itsolver(A,b,x, &prec,itparam);
    
    fasp_gettime(&solve_end);
    
    solve_duration = solve_end - solve_start;
    
    if ( prtlvl >= PRINT_MIN )
        print_cputime("Krylov method totally", setup_duration+solve_duration);
    
    // clean
#if WITH_UMFPACK
    for (l=0; l<NumLayers; l++) fasp_umfpack_free_numeric(local_LU[l]);
#endif
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
