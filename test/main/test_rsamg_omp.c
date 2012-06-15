/**
 *------------------------------------------------------
 *
 *  Modified by Zheng Lee on 05/25/2012
 *  Based on testomp.c created by Chensong Zhang
 *
 *------------------------------------------------------
 */

#include "fasp.h"
#include "fasp_functs.h"

int
main(int argc, const char * argv[])
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
    Schwarz_param   schparam; // parameters for Shcwarz method

    // Read input parameters from a disk file
	fasp_param_init("ini/openmp.dat",&inparam,&itparam,&amgparam,&iluparam,&schparam);
    
    // Set local parameters
	const int print_level   = inparam.print_level;
	const int problem_num   = inparam.problem_num;
	const int itsolver_type = inparam.solver_type;
	const int precond_type  = inparam.precond_type;
	const int output_type   = inparam.output_type;
    
    // Set output device
    if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}
    
    printf("Test Problem %d\n", problem_num);
    
    //----------------------------------------------------//
	// Step 1. Input stiffness matrix and right-hand side //
    //----------------------------------------------------//
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
    
	// Read A and b -- P1 F1 discretization for Poisson.
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

	// Print problem size
	if (print_level>PRINT_NONE) {
		printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
		printf("b: n = %d\n", b.row);
		fasp_mem_usage();
	}

	// Print out solver parameters
	if (print_level>PRINT_NONE) fasp_param_solver_print(&itparam);

    //--------------------------//
	// Step 2. Solve the system //
    //--------------------------//

	// Set initial guess
	fasp_dvec_alloc(A.row, &x); 
	fasp_dvec_set(A.row, &x, 0.0);

    // AMG as the iterative solver
	if (itsolver_type == 21) { 	
		 fasp_solver_amg(&A, &b, &x, &amgparam); 
	}
	else {
		printf("### ERROR: Wrong solver type %d!!!\n", itsolver_type);
		status = ERROR_SOLVER_TYPE;
		goto FINISHED;
	}

	fasp_mem_usage();

	if (status<0) {
		printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status);
	}
	else {
		printf("\nSolver finished successfully!\n\n", status);
	}

	if (output_type) fclose (stdout);

FINISHED:
	// Clean up memory
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
	fasp_dvec_free(&x);

	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
