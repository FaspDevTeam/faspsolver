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
	dCSRmat A;
	dvector b, uh;
    
    dBSRmat Absr;
    
    dCSRmat A11;

	int status=SUCCESS;
	
	// Step 0. Set parameters
	input_param     inparam;  // parameters from input files
	itsolver_param  itparam;  // parameters for itsolver
	AMG_param       amgparam; // parameters for AMG
	ILU_param       iluparam; // parameters for ILU
    
    // Read input parameters from a disk file
	fasp_param_init("ini/bsr.dat",&inparam,&itparam,&amgparam,&iluparam);
    
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
    
	// Step 1. Input stiffness matrix and right-hand side
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;
	
	strncpy(filename1,inparam.workdir,128);
	strncpy(filename2,inparam.workdir,128);
    
	// Read A and b -- P1 FE discretization for Poisson.
	if (problem_num == 11) {				
		datafile1="PNNL/test1/A.dat";
		strcat(filename1,datafile1);
		datafile2="PNNL/test1/b.dat";
		strcat(filename2,datafile2);
        
		fasp_dcoo_read(filename1, &A);
		fasp_dvec_read(filename2, &b);
        
        fasp_format_dcsr_dbsr(&Absr, 6, &A);        
    }	
    
    else if (problem_num == 12) {				
		datafile1="PNNL/test2/A.dat";
		strcat(filename1,datafile1);
		datafile2="PNNL/test3/b.dat";
		strcat(filename2,datafile2);
        
		fasp_dcoo_read(filename1, &A);
		fasp_dvec_read(filename2, &b);
        
        fasp_format_dcsr_dbsr(&Absr, 6, &A);
        
        // Form the right-hand-side b = A*sol
        //dvector sol = fasp_dvec_create(A.row);
        //fasp_dvec_set(A.row,&sol,1);
        //fasp_dvec_rand(A.row, &sol);
        
        //b = fasp_dvec_create(A.row);
        //fasp_blas_dcsr_mxv(&A, sol.val, b.val);     
        
        //fasp_dvec_free(&sol);
    }
    
    else if (problem_num == 13) {				
		datafile1="PNNL/test3/A.dat";
		strcat(filename1,datafile1);
		datafile2="PNNL/test3/b.dat";
		strcat(filename2,datafile2);
        
		fasp_dcoo_read(filename1, &A);
		fasp_dvec_read(filename2, &b);
        
        fasp_format_dcsr_dbsr(&Absr, 6, &A);
    }	
    
    else if (problem_num == 20) {
		datafile1="Melanie/A.dat";
		strcat(filename1,datafile1);
		fasp_dcoo_read(filename1, &A);
        
        dvector sol = fasp_dvec_create(A.row);
        fasp_dvec_rand(A.row, &sol);
        
        // Form the right-hand-side b = A*sol
        b = fasp_dvec_create(A.row);
        fasp_blas_dcsr_mxv(&A, sol.val, b.val);     
        
        /*
          int *Is, *Js;
          int nRow=A.row/5;//we must ensure this is a integer
          //printf("nRow = %d\n", nRow);
          int nCol=A.col/5;
          //printf("nCol = %d\n", nCol);
          Is=(int *)fasp_mem_calloc(nRow, sizeof(int));
          Js=(int *)fasp_mem_calloc(nCol, sizeof(int));
        
          int k;
          for(k=0;k<nRow;k++){
          Is[k]=k*5;
          Js[k]=k*5;
          }
          dCSRmat tmpMat;
          status=fasp_dcsr_getblk(&A,Is,Js,nRow,nCol,&tmpMat);
          fasp_dcsr_write("./A11.dat", &tmpMat);
          //printf("nnz = %d\n", tmpMat.nnz);
          */
        
        //fasp_dcsr_write("Amelane.dat", &A);
        //fasp_format_dcsr_dbsr(&Absr, 5, &A);
        
        fasp_dvec_free(&sol);
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
    
	// Step 2. Solve the system
    
    // Print out solver parameters
    if (print_level>PRINT_NONE) fasp_param_solver_print(&itparam);
    
    // Set initial guess
    fasp_dvec_alloc(A.row, &uh); 
    fasp_dvec_set(A.row,&uh,0.0);

    // Preconditioned Krylov methods
    if ( itsolver_type >= 1 && itsolver_type <= 20) {
        
		// Using no preconditioner for Krylov iterative methods
		if (precond_type == PREC_NULL) {
            //status = fasp_solver_dbsr_krylov(&Absr, &b, &uh, &itparam);
            status = fasp_solver_dcsr_krylov(&A, &b, &uh, &itparam);

		}	
        
		// Using diag(A) as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_DIAG) {
            status = fasp_solver_dbsr_krylov_diag(&Absr, &b, &uh, &itparam);
		}
        
		// Using AMG as preconditioner for Krylov iterative methods
		else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
			//status = fasp_solver_dcsr_krylov_amg(&A11, &b, &uh, &itparam, &amgparam);
            status = fasp_solver_dbsr_krylov_amg(&Absr, &b, &uh, &itparam,&amgparam); 
		}
        
		// Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
		else if (precond_type == PREC_ILU) {
            if (print_level>PRINT_NONE) fasp_param_ilu_print(&iluparam);
            status = fasp_solver_dbsr_krylov_ilu(&Absr, &b, &uh, &itparam,&iluparam);

		}
        
		else {
			printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);		
			exit(ERROR_SOLVER_PRECTYPE);
		}
        
	}
    
    // AMG as the iterative solver
    /*
      else if (itsolver_type == SOLVER_AMG) {
      if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
      fasp_solver_amg(&A, &b, &uh, &amgparam); 
        
      }
    */

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
    //fasp_dbsr_free(&Absr);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
    
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
