/*! \file itsolver.c
 *  \brief Iterative solvers (main file)
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_dcsr_itsolver (dCSRmat *A, dvector *b, dvector *x, 
 *                                    precond *prec, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods -- this is an abstract interface
 *
 * \param A        pointer to the coeff matrix in dCSRmat format
 * \param b        pointer to the right hand side in dvector format
 * \param x        pointer to the approx solution in dvector format
 * \param prec     pointer to the preconditioning action
 * \param itparam  pointer to parameters for iterative solvers
 *
 * \return          return the number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date 09/25/2009 
 */
INT fasp_solver_dcsr_itsolver (dCSRmat *A, 
                               dvector *b, 
                               dvector *x, 
                               precond *prec, 
                               itsolver_param *itparam)
{
	const INT   print_level   = itparam->print_level;	
	const INT   itsolver_type = itparam->itsolver_type;
	const INT   stop_type     = itparam->stop_type;
	const INT   MaxIt         = itparam->maxit;
	const INT   restart       = itparam->restart;
	const REAL  tol           = itparam->tol; 
	
    /* Local Variables */
	INT iter;
	clock_t solver_start = clock();
    
#if DEBUG_MODE
	printf("### DEBUG: itsolver ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
    
    /* Choose a desirable Krylov iterative solver */
	switch (itsolver_type) {
			
		case SOLVER_CG:
			if (print_level>0) printf("Calling PCG solver ...\n");
			iter = fasp_solver_dcsr_pcg(A, b, x, MaxIt, tol, prec, print_level, stop_type); 
            break;
			
		case SOLVER_BiCGstab:
			if (print_level>0) printf("Calling BiCGstab solver ...\n");
			iter = fasp_solver_dcsr_pbcgs(A, b, x, MaxIt, tol, prec, print_level, stop_type); 
            break;
			
		case SOLVER_MinRes:
			if (print_level>0) printf("Calling MinRes solver ...\n");		
			iter = fasp_solver_dcsr_pminres(A, b, x, MaxIt, tol, prec, print_level, stop_type); 
            break;
			
		case SOLVER_GMRES:
			if (print_level>0) printf("Calling GMRes solver ...\n");		
			iter = fasp_solver_dcsr_pgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);
            break;
			
		case SOLVER_VGMRES: 
			if (print_level>0) printf("Calling vGMRes solver ...\n");		
			iter = fasp_solver_dcsr_pvgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	
            break;
            
        case SOLVER_GCG:
			if (print_level>0) printf("Calling GCG solver ...\n");
			iter = fasp_solver_dcsr_pgcg(A, b, x, MaxIt, tol, prec, print_level, stop_type); 
            break;
            
        case SOLVER_VFGMRES: 
			if (print_level>0) printf("Calling vFGMRes solver ...\n");		
			iter = fasp_solver_dcsr_pvfgmres(A, b, x, MaxIt, tol, prec, print_level, stop_type, restart);	
            break;
			
		default:
			printf("### ERROR: Itertive solver type %d not supported!\n", itsolver_type);
			return ERROR_SOLVER_TYPE;
			
	} 
	
	if ((print_level>=PRINT_SOME) && (iter >= 0)) {
		clock_t solver_end = clock();	
		REAL solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Iterative solver costs %f seconds.\n", solver_duration);
	}
    
#if DEBUG_MODE
	printf("### DEBUG: itsolver ...... [Finish]\n");
#endif
	
	return iter;
}	

/**
 * \fn INT fasp_solver_dcsr_krylov (dCSRmat *A, dvector *b, dvector *x, 
 *                                  itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods -- no preconditioning 
 *
 * \param A        pointer to the coeff matrix in dCSRmat format
 * \param b        pointer to the right hand side in dvector format
 * \param x        pointer to the approx solution in dvector format
 * \param itparam  pointer to parameters for iterative solvers
 *
 * \return          return the number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date 09/25/2009 
 */
INT fasp_solver_dcsr_krylov (dCSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;	

    /* Local Variables */
	clock_t  solver_start, solver_end;
	INT      status = SUCCESS;
	REAL     solver_duration;
	
#if DEBUG_MODE
	printf("### DEBUG: krylov ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
	
	solver_start = clock();
	status = fasp_solver_dcsr_itsolver(A,b,x,NULL,itparam);
	
	if (print_level>=PRINT_MIN) {
		solver_end = clock();	
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Krylov method totally costs %f seconds.\n", solver_duration);
	}	
	
#if DEBUG_MODE
	printf("### DEBUG: krylov ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_diag (dCSRmat *A, dvector *b, dvector *x, itsolver_param *itparam)
 * \brief Solve Ax=b by diagonal preconditioned Krylov methods 
 * 
 * \param A        pointer to the coeff matrix in dCSRmat format
 * \param b        pointer to the right hand side in dvector format
 * \param x        pointer to the approx solution in dvector format
 * \param itparam  pointer to parameters for iterative solvers
 * \return          return the number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date 09/25/2009 
 */
INT fasp_solver_dcsr_krylov_diag (dCSRmat *A, 
                                  dvector *b, 
                                  dvector *x, 
                                  itsolver_param *itparam)
{
	const INT print_level = itparam->print_level;	

    /* Local Variables */
	clock_t   solver_start, solver_end;
	INT       status = SUCCESS;
	REAL      solver_duration;
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_diag ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
	
	solver_start = clock();

	// setup preconditioner
	dvector diag; fasp_dcsr_getdiag(0,A,&diag);	
	
	precond prec; 
    prec.data = &diag; 
    prec.fct  = fasp_precond_diag;
		
	// call iterative solver
	status = fasp_solver_dcsr_itsolver(A,b,x,&prec,itparam);
	
	if (print_level>=PRINT_MIN) {
		solver_end = clock();	
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("Diag_Krylov method totally costs %f seconds.\n", solver_duration);
	}
	
	fasp_dvec_free(&diag);

#if DEBUG_MODE
	printf("### DEBUG: krylov_diag ...... [Finish]\n");
#endif
		
	return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_amg (dCSRmat *A, dvector *b, dvector *x, 
 *                                      itsolver_param *itparam, AMG_param *amgparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods with AMG as precondition
 *
 * \param A        pointer to the coeff matrix in dCSRmat format
 * \param b        pointer to the right hand side in dvector format
 * \param x        pointer to the approx solution in dvector format
 * \param itparam  pointer to parameters for iterative solvers
 * \param amgparam AMG parameters
 *
 * \return         the number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date 09/25/2009  
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
	clock_t  solver_start, solver_end;
	INT      status = SUCCESS;
	REAL     solver_duration;
	    
#if DEBUG_MODE
	printf("### DEBUG: krylov_amg ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
	
	solver_start = clock();
	
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
			status = fasp_amg_setup_rs(mgl, amgparam); break;
	}
	
#if CHMEM_MODE	
	fasp_mem_usage();
#endif
	
	if (status < 0) goto FINISHED;
	
	// setup preconditioner
	precond_data precdata;
    fasp_param_amg_to_prec(&precdata,amgparam);
	precdata.max_levels = mgl[0].num_levels;
	precdata.mgl_data = mgl;
	
	precond prec; prec.data = &precdata;
    
	if (itparam->precond_type == PREC_FMG) {
		prec.fct = fasp_precond_famg; // Full AMG
	}
	else {
        switch (amgparam->cycle_type) {
            case AMLI_CYCLE: // AMLI cycle
                prec.fct = fasp_precond_amli; break;
            case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
                prec.fct = fasp_precond_nl_amli; break;
            default: // V,W-Cycle AMG
                prec.fct = fasp_precond_amg; break;
        }
	}
	
	// call iterative solver
	status = fasp_solver_dcsr_itsolver(A,b,x,&prec,itparam);	
    
	if (print_level>=PRINT_MIN) {
		solver_end = clock();	
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		printf("AMG_Krylov method totally costs %f seconds.\n", solver_duration);
	}
	
FINISHED:
	fasp_amg_data_free(mgl);
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_amg ...... [Finish]\n");
#endif
	
    fasp_chkerr(status, "fasp_solver_dcsr_krylov_amg");

	return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_ilu (dCSRmat *A, dvector *b, dvector *x, 
 *                                      itsolver_param *itparam, ILU_param *iluparam)
 *
 * \brief Solve Ax=b by ILUs preconditioned Krylov methods
 *
 * \param A        pointer to the coeff matrix in dCSRmat format
 * \param b        pointer to the right hand side in dvector format
 * \param x        pointer to the approx solution in dvector format
 * \param itparam  pointer to parameters for iterative solvers
 * \param iluparam pointer to parameters for ILU
 *
 * \return          return the number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date 09/25/2009 
 */
INT fasp_solver_dcsr_krylov_ilu (dCSRmat *A, 
                                 dvector *b, 
                                 dvector *x, 
                                 itsolver_param *itparam, 
                                 ILU_param *iluparam)
{
	const INT print_level = itparam->print_level;

    /* Local Variables */
	clock_t  solver_start, solver_end;
	INT      status = SUCCESS;
	REAL     solver_duration;
		
#if DEBUG_MODE
	printf("### DEBUG: krylov_ilu ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
	
	solver_start = clock();
	
	// ILU setup for whole matrix
	ILU_data LU; 
	if ( (status = fasp_ilu_dcsr_setup(A,&LU,iluparam))<0 ) goto FINISHED;
	
	// check iludata
	if ( (status = fasp_mem_iludata_check(&LU))<0 ) goto FINISHED;
	
	// set preconditioner 
	precond prec; 
	prec.data = &LU; 
    prec.fct  = fasp_precond_ilu;
	
	// call iterative solver
	status = fasp_solver_dcsr_itsolver(A,b,x,&prec,itparam);
	
	if (print_level>=PRINT_MIN) {
		solver_end = clock();	
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);
		
		switch (iluparam->ILU_type) {
			case ILUk:
				printf("ILUk_Krylov method totally costs %f seconds.\n", solver_duration);
				break;
			case ILUt:
				printf("ILUt_Krylov method totally costs %f seconds.\n", solver_duration);
				break;
			case ILUtp:
				printf("ILUtp_Krylov method totally costs %f seconds.\n", solver_duration);
				break;
			default:
				printf("ILUs_Krylov method totally costs %f seconds.\n", solver_duration);
				break;
		}
	}
	
FINISHED: 
	fasp_ilu_data_free(&LU);
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_ilu ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn INT fasp_solver_dcsr_krylov_ilu_M (dCSRmat *A, dvector *b, dvector *x, 
 *                                        itsolver_param *itparam, ILU_param *iluparam, dCSRmat *M)
 *
 * \brief Solve Ax=b by ILUs preconditioned Krylov methods: ILU of M plays as preconditioner
 *
 * \param A         pointer to the coeff matrix in dCSRmat format
 * \param b         pointer to the right hand side in dvector format
 * \param x         pointer to the approx solution in dvector format
 * \param itparam   pointer to parameters for iterative solvers
 * \param iluparam  pointer to parameters for ILU
 * \param M	        pointer to the preconditioning matrix in dCSRmat format 
 *
 * \return           return the number of iterations if succeed
 *
 * \author Xiaozhe Hu
 * \date 09/25/2009 
 *
 * \note  This function is specially designed for reservoir simulation.
 *        Have not been tested in any other places. 
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
	clock_t solver_start, solver_end;
	REAL solver_duration;
	INT status = SUCCESS;
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_ilu_M ...... [Start]\n");
	printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
	printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);		
#endif
	
	solver_start=clock();
	
	// ILU setup for M
    ILU_data LU; 
	if ( (status = fasp_ilu_dcsr_setup(M,&LU,iluparam))<0 ) goto FINISHED;
	
	// check iludata
	if ( (status = fasp_mem_iludata_check(&LU))<0 ) goto FINISHED;
	
	// set precondtioner
	precond prec; 
    prec.data = &LU; 
    prec.fct  = fasp_precond_ilu;
	
    // call iterative solver
	status = fasp_solver_dcsr_itsolver(A,b,x,&prec,itparam);
	
	if (print_level>=PRINT_MIN) {
		solver_end=clock();	
		solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);

        switch (iluparam->ILU_type) {
            case ILUk:
                printf("ILUk_Krylov method totally costs %f seconds.\n", solver_duration);
                break;
            case ILUt:
                printf("ILUt_Krylov method totally costs %f seconds.\n", solver_duration);
                break;
            case ILUtp:
                printf("ILUtp_Krylov method totally costs %f seconds.\n", solver_duration);
                break;
            default:
                printf("ILUs_Krylov method totally costs %f seconds.\n", solver_duration);
                break;		
        }
	}	
    	
FINISHED:	
	fasp_ilu_data_free(&LU);
	
#if DEBUG_MODE
	printf("### DEBUG: krylov_ilu_M ...... [Finish]\n");
#endif
	
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
