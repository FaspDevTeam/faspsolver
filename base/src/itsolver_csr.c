/*! \file itsolver_csr.c
 *
 *  \brief Iterative solvers for dCSRmat matrices
 */

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"
#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dcsr_itsolver (dCSRmat *A, dvector *b, dvector *x, 
 *                                    precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods for CSR matrices
 *
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date   09/25/2009 
 *
 * \note This is an abstract interface for iterative methods.
 */
INT fasp_solver_dcsr_itsolver (dCSRmat *A, 
                               dvector *b, 
                               dvector *x, 
                               precond *pc, 
                               itsolver_param *itparam)
{
    const SHORT  print_level   = itparam->print_level;    
    const SHORT  itsolver_type = itparam->itsolver_type;
    const SHORT  stop_type     = itparam->stop_type;
    const SHORT  restart       = itparam->restart;
    const INT    MaxIt         = itparam->maxit;
    const REAL   tol           = itparam->tol; 
    
    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT iter;
    
    fasp_gettime(&solver_start);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);
#endif
    
    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );
    
    /* Choose a desirable Krylov iterative solver */
    switch ( itsolver_type ) {
        case SOLVER_CG:
            if (print_level>PRINT_NONE) printf("\nCalling PCG solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pcg(A, b, x, pc, tol, MaxIt, stop_type, print_level);
            break;
            
        case SOLVER_BiCGstab:
            if (print_level>PRINT_NONE) printf("\nCalling PBiCGstab solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pbcgs(A, b, x, pc, tol, MaxIt, stop_type, print_level);
            break;
            
        case SOLVER_MinRes:
            if (print_level>PRINT_NONE) printf("\nCalling PMinRes solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pminres(A, b, x, pc, tol, MaxIt, stop_type, print_level);
            break;
            
        case SOLVER_GMRES:
            if (print_level>PRINT_NONE) printf("\nCalling PGMRes solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_VGMRES: 
            if (print_level>PRINT_NONE) printf("\nCalling PVGMRes solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pvgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_VFGMRES: 
            if (print_level>PRINT_NONE) printf("\nCalling PVFGMRes solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pvfgmres(A, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_GCG:
            if (print_level>PRINT_NONE) printf("\nCalling PGCG solver (CSR) ...\n");
            iter = fasp_solver_dcsr_pgcg(A, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
            return ERROR_SOLVER_TYPE;
            
    } 
    
    if ( (print_level >= PRINT_SOME) && (iter >= 0) ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return iter;
}    

/**
 * \fn INT fasp_solver_dcsr_krylov (dCSRmat *A, dvector *b, dvector *x, 
 *                                  itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods for CSR matrices 
 *
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date   09/25/2009 
 */
INT fasp_solver_dcsr_krylov (dCSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             itsolver_param *itparam)
{
    const SHORT print_level = itparam->print_level;
    
    /* Local Variables */
    INT      status = FASP_SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);
    
    status = fasp_solver_dcsr_itsolver(A,b,x,NULL,itparam);
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Krylov method totally", solver_duration);
    }    
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_diag (dCSRmat *A, dvector *b, dvector *x, 
 *                                       itsolver_param *itparam)
 *
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 * 
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date   09/25/2009 
 */
INT fasp_solver_dcsr_krylov_diag (dCSRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam)
{
    const INT print_level = itparam->print_level;    
    
    /* Local Variables */
    INT       status = FASP_SUCCESS;
    REAL      solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_diag ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);
    
    // setup preconditioner
    dvector diag; fasp_dcsr_getdiag(0,A,&diag);    
    
    precond pc; 
    pc.data = &diag; 
    pc.fct  = fasp_precond_diag;
    
    // call iterative solver
    status = fasp_solver_dcsr_itsolver(A,b,x,&pc,itparam);
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Diag_Krylov method totally", solver_duration);
    }
    
    fasp_dvec_free(&diag);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_diag ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_schwarz (dCSRmat *A, dvector *b, dvector *x, 
 *                                          itsolver_param *itparam, 
 *                                          Schwarz_param *schparam)
 *
 * \brief Solve Ax=b by overlapping schwarz Krylov methods 
 * 
 * \param A        Pointer to the coeff matrix in dCSRmat format
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for iterative solvers
 * \param schparam Pointer to parameters for Schwarz methods
 *
 * \return         Number of iterations
 *
 * \author Xiaozhe Hu
 * \date   03/21/2011
 *
 * Modified by Chensong on 07/02/2012: change interface
 */
INT fasp_solver_dcsr_krylov_schwarz (dCSRmat *A, 
                                     dvector *b, 
                                     dvector *x, 
                                     itsolver_param *itparam,
                                     Schwarz_param *schparam)
{
    const INT print_level    = itparam->print_level;	
    const INT schwarz_mmsize = schparam->schwarz_mmsize;
    const INT schwarz_maxlvl = schparam->schwarz_maxlvl;
    const INT schwarz_type   = schparam->schwarz_type;
	
    REAL setup_start, setup_end, setup_duration, solver_start, solver_end, solver_duration;
    INT status = FASP_SUCCESS;
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_schwarz ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
    
	fasp_gettime(&setup_start);
    
	// setup preconditioner
	Schwarz_data schwarz_data;
	
	// symmetrize the matrix (for now, we have to do this. We will get rid of this later)
	schwarz_data.A=fasp_dcsr_sympat(A);
	
	// construct schwarz precondtioner
	fasp_dcsr_shift (&schwarz_data.A, 1);
	fasp_schwarz_setup(&schwarz_data, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
	
	fasp_gettime(&setup_end);
	setup_duration = setup_end - setup_start;
        printf("Schwarz_Krylov method setup costs %f seconds.\n", setup_duration);
	
	precond prec;
	prec.data = &schwarz_data; 
	prec.fct = fasp_precond_schwarz;
	
	fasp_gettime(&solver_start);
	
	// solver part
	status=fasp_solver_dcsr_itsolver(A,b,x,&prec,itparam);
	
	if (print_level>PRINT_NONE) {
            fasp_gettime(&solver_end);
            solver_duration = solver_end - solver_start;
            printf("Schwarz_Krylov method totally costs %f seconds.\n", solver_duration);
	}
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_schwarz ...... [Finish]\n");
#endif
	
	fasp_schwarz_data_free(&schwarz_data);
	
	return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_amg (dCSRmat *A, dvector *b, dvector *x, 
 *                                      itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve Ax=b by AMG preconditioned Krylov methods
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 *
 * \return          Number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date   09/25/2009  
 */
INT fasp_solver_dcsr_krylov_amg (dCSRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 AMG_param *amgparam)
{
    const INT print_level = itparam->print_level;
    const INT max_levels = amgparam->max_levels;
    const INT nnz=A->nnz, m=A->row, n=A->col;    
    
    /* Local Variables */
    INT      status = FASP_SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_amg ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);
    
    // initialize A, b, x for mgl[0]    
    AMG_data *mgl=fasp_amg_data_create(max_levels);
    mgl[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
    mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n); 
    
    // setup preconditioner  
    switch (amgparam->AMG_type) {
            
        case SA_AMG: // Smoothed Aggregation AMG
            status = fasp_amg_setup_sa(mgl, amgparam); break;
        
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = fasp_amg_setup_ua(mgl, amgparam); break;
        
        default: // Classical AMG
#ifdef _OPENMP
            status = fasp_amg_setup_rs_omp(mgl, amgparam); break;
#else
            status = fasp_amg_setup_rs(mgl, amgparam); break;
#endif
            
    }
    
#if CHMEM_MODE    
    fasp_mem_usage();
#endif
    
    if (status < 0) goto FINISHED;
    
    // setup preconditioner
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,amgparam);
    pcdata.max_levels = mgl[0].num_levels;
    pcdata.mgl_data = mgl;
    
    precond pc; pc.data = &pcdata;
    
    if (itparam->precond_type == PREC_FMG) {
        pc.fct = fasp_precond_famg; // Full AMG
    }
    else {
        switch (amgparam->cycle_type) {
            case AMLI_CYCLE: // AMLI cycle
                pc.fct = fasp_precond_amli; break;
            case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
                pc.fct = fasp_precond_nl_amli; break;
            default: // V,W-Cycle AMG
                pc.fct = fasp_precond_amg; break;
        }
    }
    
    // call iterative solver
    status = fasp_solver_dcsr_itsolver(A, b, x, &pc, itparam);    
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_Krylov method totally", solver_duration);
    }
    
FINISHED:
    fasp_amg_data_free(mgl, amgparam);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_amg ...... [Finish]\n");
#endif
    return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_ilu (dCSRmat *A, dvector *b, dvector *x, 
 *                                      itsolver_param *itparam, ILU_param *iluparam)
 *
 * \brief Solve Ax=b by ILUs preconditioned Krylov methods
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param iluparam  Pointer to parameters for ILU
 *
 * \return          Number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date   09/25/2009 
 */
INT fasp_solver_dcsr_krylov_ilu (dCSRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 ILU_param *iluparam)
{
    const INT print_level = itparam->print_level;
    
    /* Local Variables */
    INT      status = FASP_SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_ilu ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);
    
    // ILU setup for whole matrix
    ILU_data LU; 
    if ( (status = fasp_ilu_dcsr_setup(A,&LU,iluparam))<0 ) goto FINISHED;
    
    // check iludata
    if ( (status = fasp_mem_iludata_check(&LU))<0 ) goto FINISHED;
    
    // set preconditioner 
    precond pc; 
    pc.data = &LU; 
    pc.fct  = fasp_precond_ilu;
    
    // call iterative solver
    status = fasp_solver_dcsr_itsolver(A,b,x,&pc,itparam);
    
    if (print_level>=PRINT_MIN) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        
        switch (iluparam->ILU_type) {
            case ILUk:
                print_cputime("ILUk_Krylov method totally", solver_duration);
                break;
            case ILUt:
                print_cputime("ILUt_Krylov method totally", solver_duration);
                break;
            case ILUtp:
                print_cputime("ILUtp_Krylov method totally", solver_duration);
                break;
            default:
                print_cputime("ILUs_Krylov method totally", solver_duration);
                break;
        }
    }
    
FINISHED: 
    fasp_ilu_data_free(&LU);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_ilu ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_ilu_M (dCSRmat *A, dvector *b, dvector *x, 
 *                                        itsolver_param *itparam, ILU_param *iluparam, 
 *                                        dCSRmat *M)
 *
 * \brief Solve Ax=b by ILUs preconditioned Krylov methods: ILU of M as preconditioner
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param iluparam  Pointer to parameters for ILU
 * \param M         Pointer to the preconditioning matrix in dCSRmat format 
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date 09/25/2009 
 *
 * \note This function is specially designed for reservoir simulation.
 *       Have not been tested in any other places. 
 */
INT fasp_solver_dcsr_krylov_ilu_M (dCSRmat *A, 
                                   dvector *b, 
                                   dvector *x, 
                                   itsolver_param *itparam, 
                                   ILU_param *iluparam, 
                                   dCSRmat *M)
{
    const INT print_level = itparam->print_level;
    
    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT status = FASP_SUCCESS;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_ilu_M ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);

    // ILU setup for M
    ILU_data LU; 
    if ( (status = fasp_ilu_dcsr_setup(M,&LU,iluparam))<0 ) goto FINISHED;
    
    // check iludata
    if ( (status = fasp_mem_iludata_check(&LU))<0 ) goto FINISHED;
    
    // set precondtioner
    precond pc; 
    pc.data = &LU; 
    pc.fct  = fasp_precond_ilu;
    
    // call iterative solver
    status = fasp_solver_dcsr_itsolver(A,b,x,&pc,itparam);
    
    if (print_level>=PRINT_MIN) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        
        switch (iluparam->ILU_type) {
            case ILUt:
                print_cputime("ILUt_Krylov method totally", solver_duration);
                break;
            case ILUtp:
                print_cputime("ILUtp_Krylov method totally", solver_duration);
                break;
            default:
                print_cputime("ILUk_Krylov method totally", solver_duration);
                break;
        }
    }    
    
FINISHED:    
    fasp_ilu_data_free(&LU);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_ilu_M ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_amg_nk (dCSRmat *A, dvector *b, dvector *x,
 *                                      itsolver_param *itparam, AMG_param *amgparam,
                                        dCSRmat *A_nk, dCSRmat *P_nk, dCSRmat *R_nk)
 *
 * \brief Solve Ax=b by AMG preconditioned Krylov methods with an extra near kernel solve
 *
 * \param A         Pointer to the coeff matrix in dCSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param amgparam  Pointer to parameters for AMG methods
 * \param A_nk      Pointer to the coeff matrix of near kernel space in dCSRmat format
 * \param P_nk      Pointer to the prolongation of near kernel space in dCSRmat format
 * \param R_nk      Pointer to the restriction of near kernel space in dCSRmat format
 *
 * \return          Number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
INT fasp_solver_dcsr_krylov_amg_nk (dCSRmat *A,
                                 dvector *b,
                                 dvector *x,
                                 itsolver_param *itparam,
                                 AMG_param *amgparam,
                                 dCSRmat *A_nk,
                                 dCSRmat *P_nk,
                                 dCSRmat *R_nk)
{
    const INT print_level = itparam->print_level;
    const INT max_levels = amgparam->max_levels;
    const INT nnz=A->nnz, m=A->row, n=A->col;
    
    /* Local Variables */
    INT      status = FASP_SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_amg ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);
#endif
    
    fasp_gettime(&solver_start);
    
    // initialize A, b, x for mgl[0]
    AMG_data *mgl=fasp_amg_data_create(max_levels);
    mgl[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
    mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n);
    
    // setup preconditioner
    switch (amgparam->AMG_type) {
            
        case SA_AMG: // Smoothed Aggregation AMG
            status = fasp_amg_setup_sa(mgl, amgparam); break;
            
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = fasp_amg_setup_ua(mgl, amgparam); break;
            
        default: // Classical AMG
#ifdef _OPENMP
            status = fasp_amg_setup_rs_omp(mgl, amgparam); break;
#else
            status = fasp_amg_setup_rs(mgl, amgparam); break;
#endif
            
    }
    
#if CHMEM_MODE
    fasp_mem_usage();
#endif
    
    if (status < 0) goto FINISHED;
    
    // setup preconditioner
    precond_data pcdata;
    fasp_param_amg_to_prec(&pcdata,amgparam);
    pcdata.max_levels = mgl[0].num_levels;
    pcdata.mgl_data = mgl;
    
    // near kernel space
#if WITH_UMFPACK // use UMFPACK directly
    dCSRmat A_tran;
    fasp_dcsr_trans(A_nk, &A_tran);
    fasp_dcsr_sort(&A_tran);
    pcdata.A_nk = &A_tran;
#else
    pcdata.A_nk = A_nk;
#endif
    pcdata.P_nk = P_nk;
    pcdata.R_nk = R_nk;
    
    precond pc; pc.data = &pcdata;
    
    pc.fct = fasp_precond_amg_nk;
    
    // call iterative solver
    status = fasp_solver_dcsr_itsolver(A, b, x, &pc, itparam);
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("AMG_NK_Krylov method totally", solver_duration);
    }
    
FINISHED:
    fasp_amg_data_free(mgl, amgparam);
#if WITH_UMFPACK // use UMFPACK directly
    fasp_dcsr_free(&A_tran);
#endif
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_dcsr_krylov_amg ...... [Finish]\n");
#endif
    return status;
}
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
