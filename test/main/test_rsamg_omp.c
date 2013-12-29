/**
 *------------------------------------------------------
 *
 *  Modified by Zheng Lee on 05/25/2012
 *  Based on testomp.c created by Chensong Zhang
 *
 *------------------------------------------------------
 */

/*! \file test_rsamg_omp.c
 *  \brief The main test function for FASP solvers: OMP version
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
	input_param     inpar;  // parameters from input files
	itsolver_param  itpar;  // parameters for itsolver
	AMG_param       amgpar; // parameters for AMG
	ILU_param       ilupar; // parameters for ILU

    // Set solver parameters: use ./ini/openmp.dat
    fasp_param_set(argc, argv, &inpar);
    fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);

    
    // Set local parameters
	const int print_level   = inpar.print_level;
	const int problem_num   = inpar.problem_num;
	const int itsolver_type = inpar.solver_type;
	const int output_type   = inpar.output_type;
    
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
	
	strncpy(filename1,inpar.workdir,128);
	strncpy(filename2,inpar.workdir,128);


	// Read A and b -- P1 F1 discretization for Poisson.
	if (problem_num == 10) {
		datafile1="csrmat_FE.dat";
		strcat(filename1,datafile1);
		datafile2="rhs_FE.dat";
		strcat(filename2,datafile2);
		fasp_dcsrvec2_read(filename1, filename2, &A, &b);
	}
	else if(problem_num == 13) {
		datafile1="csrmat_49.dat";
		strcat(filename1,datafile1);
		datafile2="rhs_49.dat";
		strcat(filename2,datafile2);
		fasp_dcsrvec2_read(filename1, filename2, &A, &b);
	}
	else if (problem_num == 11) {
		datafile1="csrmat_FD.dat";
		strcat(filename1,datafile1);
		datafile2="rhs_FD.dat";
		strcat(filename2,datafile2);
		fasp_dcsrvec2_read(filename1, filename2, &A, &b);
	}
	else if (problem_num == 100) {
		datafile1="mat_csr_100X100.dat";
		strcat(filename1,datafile1);
		datafile2="rhs_100X100.dat";
		strcat(filename2,datafile2);
		fasp_dcsrvec2_read(filename1, filename2, &A, &b);
	}
	else if (problem_num == 512) {
		datafile1="mat_csr_512X512.dat";
		strcat(filename1,datafile1);
		datafile2="rhs_512X512.dat";
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
	if (print_level>PRINT_NONE) fasp_param_solver_print(&itpar);

    //--------------------------//
	// Step 2. Solve the system //
    //--------------------------//

	// Set initial guess
	fasp_dvec_alloc(A.row, &x); 
	fasp_dvec_set(A.row, &x, 0.0);

    // AMG as the iterative solver
	if (itsolver_type == 21) { 	
		 fasp_solver_amg(&A, &b, &x, &amgpar); 
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
		printf("\nSolver finished successfully!\n\n");
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
