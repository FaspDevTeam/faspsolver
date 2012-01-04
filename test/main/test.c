/**
 *		Test FASP solvers with a few simple problems. 
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/31/2009.
 *      Modified by Chensong Zhang on 09/09/2011.
 *
 *------------------------------------------------------
 *
 */

/*! \file test.c
 *  \brief The main test function for FASP solvers
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
	fasp_param_init("ini/input.dat",&inparam,&itparam,&amgparam,&iluparam);
    
    // Set local parameters
	const int print_level   = inparam.print_level;
	const int problem_num   = inparam.problem_num;
	const int itsolver_type = inparam.itsolver_type;
	const int precond_type  = inparam.precond_type;
	const int output_type   = inparam.output_type;
    
    // Set output device
    if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}
    
    printf("Test Problem %d\n", problem_num);
    
	//! Step 1. Assemble stiffness matrix and right-hand side
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
    
	// Read A and b -- P1 FE discretization for Poisson.
	if (problem_num == 10) {				
		datafile1="matP1.dat";
		strcat(filename1,datafile1);
		datafile2="rhsP1.dat";
		strcat(filename2,datafile2);
		fasp_dcoo_read(filename1, &A);
		fasp_dvecind_read(filename2, &b);
	}	
    
	// Read A and b -- P1 FE discretization for Poisson, large    
    else if (problem_num == 11) {
		datafile1="cooA_1046529.dat";
		strcat(filename1,datafile1);
		fasp_dcoo_read(filename1, &A);
        
        dvector sol = fasp_dvec_create(A.row);
        fasp_dvec_rand(A.row, &sol);
        
        // Form the right-hand-side b = A*sol
        b = fasp_dvec_create(A.row);
        fasp_blas_dcsr_mxv(&A, sol.val, b.val);     
        
        fasp_dvec_free(&sol);
	}	
    
	// Assemble A and b -- P1 FE discretization for Poisson.
	else if (problem_num == 19) {	
        assemble(&A,&b, 10);
        
        fasp_dcsr_compress_inplace(&A, SMALLREAL);
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
    
	//! Step 2. Solve the system
    
    // Print out solver parameters
    if (print_level>PRINT_NONE) fasp_param_solver_print(&itparam);
    
    // Set initial guess
    fasp_dvec_alloc(A.row, &uh); 
    fasp_dvec_set(A.row,&uh,0.0);
	
    // AMG as the iterative solver
	if (itsolver_type == 0) {
        
        if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
		status = fasp_solver_amg(&A, &b, &uh, &amgparam); 
    
	}
    
    // Preconditioned Krylov methods
	else if ( (itsolver_type >= 1 && itsolver_type <= 5) || (itsolver_type == 9) || (itsolver_type == 10)) {
        
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
            if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
			status = fasp_solver_dcsr_krylov_amg(&A, &b, &uh, &itparam, &amgparam);
		}
        
		// Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
		else if (precond_type == PREC_ILU) {
            if (print_level>PRINT_NONE) fasp_param_ilu_print(&iluparam);
			status = fasp_solver_dcsr_krylov_ilu(&A, &b, &uh, &itparam, &iluparam);
		}
        
		else {
			printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);		
			exit(ERROR_SOLVER_PRECTYPE);
		}
        
	}
    
#if With_SuperLU // use SuperLU directly
	else if (itsolver_type == 6) {
		status = superlu(&A, &b, &uh, print_level);	 
	}
#endif	 
	
#if With_UMFPACK // use UMFPACK directly
	else if (itsolver_type == 7) {
		status = umfpack(&A, &b, &uh, print_level);	 
	}
#endif	 
    
    // Full AMG as the iterative solver 
    else if (itsolver_type == 8) {
        if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
        status = fasp_solver_famg(&A, &b, &uh, &amgparam);
    }
    
	else {
		printf("### ERROR: Wrong solver type %d!!!\n", itsolver_type);		
		exit(ERROR_SOLVER_TYPE);
	}
	
	fasp_mem_usage();
	
	if (status<0) {
		printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status);
	}
	else {
		printf("\nSolver finished!\n\n", status);
	}
    
    if (output_type) fclose (stdout);
    
    // Clean up memory
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
    
	return SUCCESS;
}
