/*! \file testbsr.c
 *  \brief The main test function for FASP solvers -- BSR format
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for a few simple tests.
 *
 * \author Chensong Zhang
 * \date   03/31/2009
 * 
 * Modified by Chensong Zhang on 09/09/2011
 */
int main (int argc, const char * argv[]) 
{
    dBSRmat Absr;
	dvector b, uh;
    
	int status=SUCCESS;
	
	// Step 0. Set parameters
	input_param     inparam;  // parameters from input files
	itsolver_param  itparam;  // parameters for itsolver
	AMG_param       amgparam; // parameters for AMG
	ILU_param       iluparam; // parameters for ILU
    Schwarz_param   schparam; // parameters for Shcwarz method
    
    // Read input parameters from a disk file
	fasp_param_init("ini/bsr.dat",&inparam,&itparam,&amgparam,&iluparam,&schparam);
    
    // Set local parameters
	const int print_level   = inparam.print_level;
	const int problem_num   = inparam.problem_num;
	const int itsolver_type = inparam.solver_type;
	const int precond_type  = inparam.precond_type;
	const int output_type   = inparam.output_type;
    
    // Set output devices
    if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}
    
    printf("Test Problem %d\n", problem_num);
    
	// Step 1. Input stiffness matrix and right-hand side
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
    
    // Default test problem from black-oil benchmark: SPE01
	if (problem_num == 10) {				
        // Read the stiffness matrix from bsrmat_SPE01.dat
        strncpy(filename1,inparam.workdir,128);    
        datafile1="bsrmat_SPE01.dat"; strcat(filename1,datafile1);
        fasp_dbsr_read(filename1, &Absr);
        
        // Read the RHS from rhs_SPE01.dat
        strncpy(filename2,inparam.workdir,128);
        datafile2="rhs_SPE01.dat"; strcat(filename2,datafile2);
        fasp_dvec_read(filename2, &b);
    }
    
    else {
		printf("### ERROR: Unrecognized problem number %d\n", problem_num);
		return ERROR_INPUT_PAR;
	}
    
    // Print problem size
	if (print_level>PRINT_NONE) {
        printf("A: m = %d, n = %d, nnz = %d\n", Absr.ROW, Absr.COL, Absr.NNZ);
        printf("b: n = %d\n", b.row);
	}
    
	// Step 2. Solve the system
    
    // Print out solver parameters
    if (print_level>PRINT_NONE) fasp_param_solver_print(&itparam);
    
    // Set initial guess
    fasp_dvec_alloc(b.row, &uh); 
    fasp_dvec_set(b.row, &uh, 0.0);
    
    // Preconditioned Krylov methods
    if ( itsolver_type > 0 && itsolver_type < 20 ) {
        
		// Using no preconditioner for Krylov iterative methods
		if (precond_type == PREC_NULL) {
            status = fasp_solver_dbsr_krylov(&Absr, &b, &uh, &itparam);
		}	
        
		// Using diag(A) as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_DIAG) {
            status = fasp_solver_dbsr_krylov_diag(&Absr, &b, &uh, &itparam);
		}
        
		// Using AMG as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
            status = fasp_solver_dbsr_krylov_amg(&Absr, &b, &uh, &itparam, &amgparam); 
		}
        
		// Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
		else if (precond_type == PREC_ILU) {
            if (print_level>PRINT_NONE) fasp_param_ilu_print(&iluparam);
            status = fasp_solver_dbsr_krylov_ilu(&Absr, &b, &uh, &itparam, &iluparam);
		}
        
		else {
			printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);		
			exit(ERROR_SOLVER_PRECTYPE);
		}
        
	}
    
	else {
		printf("### ERROR: Wrong solver type %d!!!\n", itsolver_type);		
		status = ERROR_SOLVER_TYPE;
        goto FINISHED;
	}
    	
	if (status<0) {
		printf("\n### WARNING: Solver failed! Exit status = %d.\n\n", status);
	}
	else {
		printf("\nSolver finished successfully!\n\n");
	}
    
    if (output_type) fclose (stdout);
    
 FINISHED:
    // Clean up memory
    fasp_dbsr_free(&Absr);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
    
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
