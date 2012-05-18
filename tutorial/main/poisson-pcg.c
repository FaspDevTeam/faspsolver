/*! \file poisson-pcg.c
 *  \brief The third test example for FASP: using PCG to solve 
 *         the discrete Poisson equation from P1 finite element
 *
 *  \note  PCG example for FASP: C version
 *
 *  Solving the Poisson equation (P1 FEM) with preconditioned CG methods
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for the third example.
 *
 * \author Feiteng Huang
 * \date   05/17/2012
 */
int main (int argc, const char * argv[]) 
{
    input_param          inparam;  // parameters from input files
    itsolver_param       itparam;  // parameters for itsolver
    AMG_param            amgparam; // parameters for AMG
    ILU_param            iluparam; // parameters for ILU

    printf("\n========================================");
    printf("\n||   FASP: PCG example -- C version   ||");
    printf("\n========================================\n\n");
    
    // Step 0. Set parameters
    // Read input and precondition parameters from a disk file
    // In this example, we read everything from a disk file:
    //          "./ini/pcg.dat"
    // See the reference manual for details of the parameters. 
    fasp_param_init("ini/pcg.dat",&inparam,&itparam,&amgparam,&iluparam);

    // Set local parameters
    const SHORT print_level = itparam->print_level;    
    const SHORT pc_type     = inparam.precond_type;
    const SHORT stop_type   = itparam.stop_type;
    const INT   maxit       = itparam.maxit;
    const REAL  tol         = itparam.tol;
    
    // Step 1. Get stiffness matrix and right-hand side
    // Read A and b -- P1 FE discretization for Poisson.
    // The location of the data files are given in "pcg.dat".
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
    
    // Step 2. Print problem size and PCG parameters
    if (print_level>PRINT_NONE) {
        printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
        printf("b: n = %d\n", b.row);
        fasp_param_solver_print(&itparam);
    }
    
    // Setp 3. Setup preconditioner
    // Preconditioiner type is determined by pc_type
    precond *pc = fasp_precond_setup(pc_type, &amgparam, &iluparam, &A);

    // Step 4. Solve the system with PCG as an iterative solver
    // Set the initial guess to be zero and then solve it using pcg solver
    // Note that we call PCG interface directly. There is another way which
    // calls the abstract iterative method interface; see possion-its.c for
    // more details. 
    fasp_dvec_alloc(A.row, &x); fasp_dvec_set(A.row, &x, 0.0);
    fasp_solver_dcsr_pcg(&A, &b, &x, pc, tol, maxit, stop_type, print_level);
    
    // Step 5. Clean up memory
    if (pc_type!=PREC_NULL) fasp_mem_free(pc->data);
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);
    
    return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
