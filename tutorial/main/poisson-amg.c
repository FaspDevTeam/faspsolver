/*! \file poisson-amg.c
 *  \brief The first test example for FASP: using AMG to solve 
 *         the discrete Poisson equation from P1 finite element
 *
 *	\note  AMG example for FASP: C version
 *
 *  Solving the Poisson equation (P1 FEM) with AMG
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for the first example.
 *
 * \author Chensong Zhang
 * \date   12/21/2011
 */
int main (int argc, const char * argv[]) 
{
    input_param     inparam;  // parameters from input files
    AMG_param       amgparam; // parameters for AMG
    
    printf("\n========================================");
    printf("\n||   FASP: AMG example -- C version   ||");
    printf("\n========================================\n\n");

    // Step 0. Set parameters
    // Read input and AMG parameters from a disk file
    // In this example, we read everything from a disk file:
    //          "./ini/amg.dat"
    // See the reference manual for details of the parameters. 
    fasp_param_init("ini/amg.dat",&inparam,NULL,&amgparam,NULL);
    
    // Set local parameters using the input values
    const int print_level   = inparam.print_level;
    const int problem_num   = inparam.problem_num;
    const int itsolver_type = inparam.itsolver_type;
    const int precond_type  = inparam.precond_type;
    const int output_type   = inparam.output_type;
    
    // Step 1. Get stiffness matrix and right-hand side
    // Read A and b -- P1 FE discretization for Poisson.
    // The location of the data files are given in "amg.dat".
    dCSRmat A;
    dvector b, x;
    char filename1[512], *datafile1;
    char filename2[512], *datafile2;
	
    // Read the stiffness matrix from matP1.dat
    strncpy(filename1,inparam.workdir,128);    
    datafile1="matP1.dat"; strcat(filename1,datafile1);
    fasp_dcoo_read(filename1, &A);

    // Read the RHS from rhsP1.dat
    strncpy(filename2,inparam.workdir,128);
    datafile2="rhsP1.dat"; strcat(filename2,datafile2);
    fasp_dvecind_read(filename2, &b);    
    
    // Step 2. Print problem size and AMG parameters
    if (print_level>PRINT_NONE) {
        printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
        printf("b: n = %d\n", b.row);
        fasp_param_amg_print(&amgparam);
    }
    
    // Step 3. Solve the system with AMG as an iterative solver
    // Set the initial guess to be zero and then solve it
    // with AMG method as an iterative procedure
    fasp_dvec_alloc(A.row, &x); fasp_dvec_set(A.row,&x,0.0);
    fasp_solver_amg(&A, &b, &x, &amgparam);
    
    // Step 4. Clean up memory
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);
    
    return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
