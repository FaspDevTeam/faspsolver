/*! \file poisson-its.c
 *  \brief The second test example for FASP: using ITS to solve 
 *         the discrete Poisson equation from P1 finite element
 *
 *  \note  ITS example for FASP: C version
 *
 *  Solving the Poisson equation (P1 FEM) with iterative methods
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for the second example.
 *
 * \author Feiteng Huang
 * \date   04/13/2012
 */
int main (int argc, const char * argv[]) 
{
    input_param          inparam;  // parameters from input files
    itsolver_param       itparam;  // parameters for itsolver

    printf("\n========================================");
    printf("\n||   FASP: ITS example -- C version   ||");
    printf("\n========================================\n\n");
    
    // Step 0. Set parameters
    // Read input and AMG parameters from a disk file
    // In this example, we read everything from a disk file:
    //          "./ini/its.dat"
    // See the reference manual for details of the parameters. 
    fasp_param_init("ini/its.dat",&inparam,&itparam,NULL,NULL);

    // Set local parameters
    const int print_level   = inparam.print_level;
    
    // Step 1. Get stiffness matrix and right-hand side
    // Read A and b -- P1 FE discretization for Poisson.
    // The location of the data files are given in "its.dat".
    dCSRmat A;
    dvector b, x;
    char filename1[512], *datafile1;
    char filename2[512], *datafile2;

    // Stiffness matrix from matP1.dat
    strncpy(filename1,inparam.workdir,128);    
    datafile1="matP1.dat"; strcat(filename1,datafile1);
    fasp_dcoo_read(filename1, &A);

    // RHS from rhsP1.dat
    strncpy(filename2,inparam.workdir,128);
    datafile2="rhsP1.dat"; strcat(filename2,datafile2);
    fasp_dvecind_read(filename2, &b);    
    
    // Step 2. Print problem size and ITS parameters
    if (print_level>PRINT_NONE) {
        printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
        printf("b: n = %d\n", b.row);
        fasp_param_solver_print(&itparam);
    }
    
    // Step 3. Solve the system with ITS as an iterative solver
    // Set the initial guess to be zero and then solve it using standard
    // iterative methods, without applying any preconditioners
    fasp_dvec_alloc(A.row, &x); fasp_dvec_set(A.row,&x,0.0);
    fasp_solver_dcsr_itsolver(&A, &b, &x, NULL, &itparam);
    
    // Step 4. Clean up memory
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);
    
    return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
