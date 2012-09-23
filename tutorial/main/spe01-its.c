/*! \file spe01-its.c
 *  \brief The fourth test example for FASP: using ITS_BSR to solve 
 *         the Jacobian equation from reservoir simulation benchmark
 *         problem SPE01.
 *
 *  \note  ITS_BSR example for FASP: C version
 *
 *  Solving the SPE01 benchmark with iterative methods
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for the fourth example.
 *
 * \author Feiteng Huang, Chensong Zhang
 * \date   05/22/2012
 *
 * Modified by Chensong Zhang on 09/22/2012
 */
int main (int argc, const char * argv[]) 
{
    input_param          inparam;  // parameters from input files
    itsolver_param       itparam;  // parameters for itsolver
    ILU_param            iluparam; // parameters for ILU
    
    printf("\n========================================");
    printf("\n||   FASP: SPE01 -- ITS BSR version   ||");
    printf("\n========================================\n\n");
    
    // Step 0. Set parameters
    // Read input parameters from a disk file
    // In this example, we read everything from a disk file:
    //          "ini/its_bsr.dat"
    // See the reference manual for details of the parameters. 
    fasp_param_init("ini/its_bsr.dat",&inparam,&itparam,NULL,&iluparam,NULL);
    
    // Set local parameters
    const int print_level = inparam.print_level;
    
    // Step 1. Get stiffness matrix and right-hand side
    // Read A and b -- P1 FE discretization for Poisson. The location
    // of the data files is given in "its.dat".
    dBSRmat A;
    dvector b, x;
    char filename1[512], *datafile1;
    char filename2[512], *datafile2;
    
    // Read the stiffness matrix from bsrmat_SPE01.dat
    strncpy(filename1,inparam.workdir,128);    
    datafile1="bsrmat_SPE01.dat"; strcat(filename1,datafile1);
    fasp_dbsr_read(filename1, &A);
    
    // Read the RHS from rhs_SPE01.dat
    strncpy(filename2,inparam.workdir,128);
    datafile2="rhs_SPE01.dat"; strcat(filename2,datafile2);
    fasp_dvec_read(filename2, &b);
    
    // Step 2. Print problem size and ITS_bsr parameters
    if (print_level>PRINT_NONE) {
        printf("A: m = %d, n = %d, nnz = %d\n", A.ROW, A.COL, A.NNZ);
        printf("b: n = %d\n", b.row);
        fasp_param_solver_print(&itparam);
        fasp_param_ilu_print(&iluparam);
    }
    
    // Step 3. Solve the system with ITS_BSR as an iterative solver
    // Set the initial guess to be zero and then solve it using standard
    // iterative methods, without applying any preconditioners
    fasp_dvec_alloc(b.row, &x);
    fasp_dvec_set(b.row,&x,0.0);
    
    fasp_solver_dbsr_krylov_ilu(&A, &b, &x, &itparam, &iluparam);
    
    // Step 4. Clean up memory
    fasp_dbsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);
    
    return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
