/*! \file testomp.c
 *  \brief The main test function for FASP OpenMP solvers
 * 
 *  \note  This file becomes obsolete as test.c can do the same job now!
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for a few simple OMP tests.
 *
 * \author Chensong Zhang
 * \date   10/06/2011
 *
 * Modified by Chesong Zhang on 09/29/2013: revise according to test.c
 */
int main (int argc, const char * argv[]) 
{
	dCSRmat A;
	dvector b, x;
	int status=SUCCESS;
	
    //------------------------//
	// Step 0. Set parameters //
    //------------------------//
	input_param     inparam;  // parameters from input files
	itsolver_param  itparam;  // parameters for itsolver
	AMG_param       amgparam; // parameters for AMG
	ILU_param       iluparam; // parameters for ILU
    
    // Read input parameters from a disk file
	char *inputfile = "ini/openmp.dat";
	fasp_param_init(inputfile,&inparam,&itparam,&amgparam,&iluparam,NULL);
	
    // Set local parameters
	const int print_level   = inparam.print_level;
	const int problem_num   = inparam.problem_num;
	const int solver_type   = inparam.solver_type;
	const int precond_type  = inparam.precond_type;
	const int output_type   = inparam.output_type;
	
	if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}
    	
    //----------------------------------------------------//
	// Step 1. Input stiffness matrix and right-hand side //
    //----------------------------------------------------//
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
    
	printf("Test Problem %d\n", problem_num);

	// Read A and b -- P1 FE discretization for Poisson.
	if (problem_num == 10) {				
		datafile1="csrmat_FE.dat";
		strcat(filename1,datafile1);
		datafile2="rhs_FE.dat";
		strcat(filename2,datafile2);        
		fasp_dcsrvec2_read(filename1, filename2, &A, &b);
	}		
 
	else {
		printf("### ERROR: Unrecognized problem number %d\n", problem_num);
		return ERROR_INPUT_PAR;
	}
	
	fasp_mem_usage();
    
    // Print problem size
    printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
    printf("b: n = %d\n", b.row);

	// Solve the system
	if (print_level > PRINT_NONE) {
		printf("Max it num = %d\n", inparam.itsolver_maxit);
		printf("Tolerance  = %e\n", inparam.itsolver_tol);
	}
	   
    //--------------------------//
	// Step 2. Solve the system //
    //--------------------------//

	// Initial guess
    fasp_dvec_alloc(A.row, &x); 
    fasp_dvec_set(A.row,&x,0.0);
	
    // Preconditioned Krylov methods
    if ( solver_type >= 1 && solver_type <= 20) {
        
		// Using no preconditioner for Krylov iterative methods
		if (precond_type == PREC_NULL) {
			status = fasp_solver_dcsr_krylov(&A, &b, &x, &itparam);
		}
        
		// Using diag(A) as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_DIAG) {
			status = fasp_solver_dcsr_krylov_diag(&A, &b, &x, &itparam);
		}
        
		// Using AMG as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level > PRINT_NONE) fasp_param_amg_print(&amgparam);
			status = fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
		}
        
		// Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
		else if (precond_type == PREC_ILU) {
            if (print_level > PRINT_NONE) fasp_param_ilu_print(&iluparam);
			status = fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &itparam, &iluparam);
		}
        
		else {
			printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);
			status = ERROR_SOLVER_PRECTYPE;
		}
        
	}
    
    // AMG as the iterative solver
	else if (solver_type == SOLVER_AMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgparam);
		fasp_solver_amg(&A, &b, &x, &amgparam);
        
	}
    
    // Full AMG as the iterative solver
    else if (solver_type == SOLVER_FMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgparam);
        fasp_solver_famg(&A, &b, &x, &amgparam);
    }

	else {
		printf("### ERROR: Wrong solver type %d!!!\n", solver_type);		
		exit(ERROR_SOLVER_TYPE);
	}
	
	fasp_mem_usage();
	
	// Output solution to a diskfile
	if (status<0) {
		printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status);
	}
	else {
		printf("\nSolver succeed!\n\n");
	}
    
    if (output_type) fclose (stdout);

    // Clean up memory
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
	fasp_dvec_free(&x);
		
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
