/**
 *		Test FASP solvers with a few simple problems. 
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 10/06/2011.
 *
 *------------------------------------------------------
 *
 */

/*! \file testomp.c
 *  \brief The main test function for FASP OpenMP solvers
 */

#include "fasp.h"
#include "fasp_functs.h"

// Declare P1 assembling function
extern int assemble(dCSRmat *ptr_A, dvector *ptr_b, int levelNum);

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for a few simple tests.
 */
int main (int argc, const char * argv[]) 
{
	dCSRmat A;
	dvector b, uh;
	int status=SUCCESS;
	
	//! Step 0. Set parameters
	input_param     inparam;  // parameters from input files
	itsolver_param  itparam;  // parameters for itsolver
	AMG_param       amgparam; // parameters for AMG
	ILU_param       iluparam; // parameters for ILU
    
    // Read input parameters from a disk file
	char *inputfile = "ini/testomp.dat";
	fasp_param_init(inputfile,&inparam,&itparam,&amgparam,&iluparam);
	
    // Set local parameters
	const int print_level   = inparam.print_level;
	const int problem_num   = inparam.problem_num;
	const int itsolver_type = inparam.itsolver_type;
	const int precond_type  = inparam.precond_type;
	const int output_type   = inparam.output_type;
	
	if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}
    	
	//! Step 1. Assemble or read matrix and right-hand side
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
    
	printf("Test Problem %d\n", problem_num);

	// Read A and b -- P1 FE discretization for Poisson.
	if (problem_num == 10) {				
		datafile1="matP1.dat";
		strcat(filename1,datafile1);
		datafile2="rhsP1.dat";
		strcat(filename2,datafile2);
		fasp_dcoo_read(filename1, &A);
		fasp_dvecind_read(filename2, &b);
	}		
 
	// Assemble A and b -- P1 FE discretization for Poisson.
	else if (problem_num == 19) {	
	  	assemble(&A,&b,9);
	}

	else {
		printf("### ERROR: Unrecognized problem number %d\n", problem_num);
		return ERROR_INPUT_PAR;
	}
	
	fasp_mem_usage();
    
    // Print problem size
    printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
    printf("b: n = %d\n", b.row);

	//! Step 2. Solve the system
	if (print_level>0) {
		printf("Max it num = %d\n", inparam.itsolver_maxit);
		printf("Tolerance  = %e\n", inparam.itsolver_tol);
	}
	    
	// initial guess
    fasp_dvec_alloc(A.row, &uh); 
    fasp_dvec_set(A.row,&uh,0.0);
	
    // AMG as the iterative solver
	if (itsolver_type == 0) { 	
		status = fasp_solver_amg(&A, &b, &uh, &amgparam); 
	}

    // Full AMG as the iterative solver 
	else if (itsolver_type == 8) {
		status = fasp_solver_famg(&A, &b, &uh, &amgparam);
	}
	
#if FASP_USE_OPENMP
    // OMP version AMG as the iterative solver
	else if( itsolver_type == 110) {        
        int nts = 1;
 		printf("omp test itsolver _ type = %d amgparam.max_iter = %d, amgparam.tol = %lf\n",
                itsolver_type, amgparam.max_iter, amgparam.tol);
		omp_set_num_threads(nts);
		status = fasp_solver_amg_omp(&A, &b, &uh, &amgparam, nts, 1000);
	}
#endif	

	else if (itsolver_type >= 1 && itsolver_type <= 5) {
        
		// Using no preconditioner for Krylov iterative methods
		if (precond_type == PREC_NULL) {
			status = fasp_solver_dcsr_krylov(&A, &b, &uh, &itparam);
		}	
        
		// Using diag(A) as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_DIAG) {
			status = fasp_solver_dcsr_krylov_diag(&A, &b, &uh, &itparam);
		}
        
		// Using AMG as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
			status = fasp_solver_dcsr_krylov_amg(&A, &b, &uh, &itparam, &amgparam);
		}
        
		// Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
		else if (precond_type == PREC_ILU) {
			status = fasp_solver_dcsr_krylov_ilu(&A, &b, &uh, &itparam, &iluparam);
		}
        
		else {
			printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);		
			exit(ERROR_SOLVER_PRECTYPE);
		}
		
	}

	else {
		printf("### ERROR: Wrong solver type %d!!!\n", itsolver_type);		
		exit(ERROR_SOLVER_TYPE);
	}
	
	fasp_mem_usage();
	
	// output solution to a diskfile
	if (status<0) {
		printf("\n### WARNING Solver failed! Exit status = %d.\n\n", status);
	}
	else {
		printf("\nSolver succeed! Exit status = %d.\n\n", status);
	}
    
    if (output_type) fclose (stdout);

	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
		
	return SUCCESS;
}
