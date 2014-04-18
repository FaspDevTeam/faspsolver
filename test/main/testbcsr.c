/*! \file testbcsr.c
 *  \brief The main test function for FASP solvers -- BCSR format
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for a few simple tests.
 *
 * \author Xiaozhe Hu
 * \date   04/07/2014
 * 
 * Modified by Chensong Zhang on 09/09/2011
 */
int main (int argc, const char * argv[]) 
{
    block_dCSRmat Abcsr;
	dvector b, uh;
    
    dBSRmat Absr;
    
	int status=SUCCESS;
	
	// Step 0. Set parameters
	input_param     inpar;  // parameters from input files
	itsolver_param  itpar;  // parameters for itsolver
	AMG_param       amgpar; // parameters for AMG
	ILU_param       ilupar; // parameters for ILU
    
    // Set solver parameters: use ./ini/bsr.dat
    fasp_param_set(argc, argv, &inpar);
    fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);
    
    // Set local parameters
	const int print_level   = inpar.print_level;
	const int problem_num   = inpar.problem_num;
	const int itsolver_type = inpar.solver_type;
	const int precond_type  = inpar.precond_type;
	const int output_type   = inpar.output_type;
    
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
	
	strncpy(filename1,inpar.workdir,128);
	strncpy(filename2,inpar.workdir,128);
    
    // Default test problem from black-oil benchmark: SPE01
	if (problem_num == 10) {
        
        strncpy(filename1,inpar.workdir,128);
        datafile1="/pnp-data/A.dat";
        strcat(filename1,datafile1);

        dCSRmat A;
        fasp_dcoo_shift_read(filename1, &A);
        
        // form block CSR matrix
        INT i;
        INT row = A.row/3;
		ivector phi_idx;
        ivector n_idx;
        ivector p_idx;
        
        fasp_ivec_alloc(row, &phi_idx);
        fasp_ivec_alloc(row, &n_idx);
        fasp_ivec_alloc(row, &p_idx);
        
        printf("row = %d\n",row);
        for (i=0; i<row; i++){
            phi_idx.val[i] = 3*i;
            n_idx.val[i] = 3*i+1;
            p_idx.val[i] = 3*i+2;
        }

        // Assemble the matrix in block dCSR format
		Abcsr.brow = 3; Abcsr.bcol = 3;
		Abcsr.blocks = (dCSRmat **)calloc(9, sizeof(dCSRmat *));
        for (i=0; i<9 ;i++) {
            Abcsr.blocks[i] = (dCSRmat *)fasp_mem_calloc(1, sizeof(dCSRmat));
        }
		      
		// A11
        fasp_dcsr_getblk(&A, phi_idx.val, phi_idx.val, row, row, Abcsr.blocks[0]);
        // A12
        fasp_dcsr_getblk(&A, phi_idx.val, n_idx.val, row, row, Abcsr.blocks[1]);
        // A13
        fasp_dcsr_getblk(&A, phi_idx.val, p_idx.val, row, row, Abcsr.blocks[2]);
        // A21
        fasp_dcsr_getblk(&A, n_idx.val, phi_idx.val, row, row, Abcsr.blocks[3]);
        // A22
        fasp_dcsr_getblk(&A, n_idx.val, n_idx.val, row, row, Abcsr.blocks[4]);
        // A23
        fasp_dcsr_getblk(&A, n_idx.val, p_idx.val, row, row, Abcsr.blocks[5]);
        // A31
        fasp_dcsr_getblk(&A, p_idx.val, phi_idx.val, row, row, Abcsr.blocks[6]);
        // A32
        fasp_dcsr_getblk(&A, p_idx.val, n_idx.val, row, row, Abcsr.blocks[7]);
        // A33
        fasp_dcsr_getblk(&A, p_idx.val, p_idx.val, row, row, Abcsr.blocks[8]);
        
        // form right hand side
        dvector b_temp;
        strncpy(filename2,inpar.workdir,128);
        datafile2="/pnp-data/rhs.dat"; strcat(filename2,datafile2);
        fasp_dvec_read(filename2, &b_temp);
        
        
        fasp_dvec_alloc(b_temp.row, &b);
        for (i=0; i<row; i++){
            
            b.val[i]        = b_temp.val[3*i];
            b.val[row+i]    = b_temp.val[3*i+1];
            b.val[2*row+i]  = b_temp.val[3*i+2];
            
        }
        
        //Absr = fasp_format_dcsr_dbsr(&A, 3);
        
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
    if (print_level>PRINT_NONE) fasp_param_solver_print(&itpar);
    
    // Set initial guess
    fasp_dvec_alloc(b.row, &uh); 
    fasp_dvec_set(b.row, &uh, 0.0);
    
    // Preconditioned Krylov methods
    if ( itsolver_type > 0 && itsolver_type < 20 ) {
        
		// Using no preconditioner for Krylov iterative methods
		if (precond_type == PREC_NULL) {
            status = fasp_solver_bdcsr_krylov(&Abcsr, &b, &uh, &itpar);
		}	
        
		// Using diag(A) as preconditioner for Krylov iterative methods
		else if (precond_type > 20 &&  precond_type < 30) {
            status = fasp_solver_bdcsr_krylov_block(&Abcsr, &b, &uh, &itpar, &amgpar);
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
    fasp_bdcsr_free(&Abcsr);
    fasp_dbsr_free(&Absr);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
    
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
