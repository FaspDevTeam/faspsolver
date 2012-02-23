/**
 *		First example for FASP:
 *
 *      Solving the Poisson equation (P1 FEM) 
 *      with AMG methods
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 12/21/2011.
 *
 *------------------------------------------------------
 *
 */

/*! \file poisson-amg.c
 *  \brief The first test example for FASP: using AMG to solve 
 *         the discrete Poisson equation from P1 finite element.
 */

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for the first example.
 */
int main (int argc, const char * argv[]) 
{
    //! Step 0. Set parameters
    //!         We read everything from a disk file "amg.dat"
    input_param     inparam;  // parameters from input files
    itsolver_param  itparam;  // parameters for itsolver
    AMG_param       amgparam; // parameters for AMG
    ILU_param       iluparam; // parameters for ILU
    
    // Read input parameters from a disk file
    fasp_param_init("ini/amg.dat",&inparam,&itparam,&amgparam,&iluparam);
    
    // Set local parameters
    const int print_level   = inparam.print_level;
    const int problem_num   = inparam.problem_num;
    const int itsolver_type = inparam.itsolver_type;
    const int precond_type  = inparam.precond_type;
    const int output_type   = inparam.output_type;
    
    //! Step 1. Get stiffness matrix and right-hand side
    //!         Both stiffness matrix and right-hand side are written on disk files
    dCSRmat A;
    dvector b, uh;
    char filename1[512], *datafile1;
    char filename2[512], *datafile2;
	
    strncpy(filename1,inparam.workdir,128);
    strncpy(filename2,inparam.workdir,128);
    
    //! Read A and b -- P1 FE discretization for Poisson.
    datafile1="matP1.dat";
    strcat(filename1,datafile1);
    datafile2="rhsP1.dat";
    strcat(filename2,datafile2);
    fasp_dcoo_read(filename1, &A);
    fasp_dvecind_read(filename2, &b);    
    
    // Print problem size
    if (print_level>PRINT_NONE) {
        printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
        printf("b: n = %d\n", b.row);
        fasp_mem_usage();
    }
    
    //! Step 2. Solve the system with AMG
    
    // Set initial guess
    fasp_dvec_alloc(A.row, &uh); 
    fasp_dvec_set(A.row,&uh,0.0);
	
    // Print out solver parameters
    if (print_level>PRINT_NONE) fasp_param_amg_print(&amgparam);
    
    // AMG as the iterative solver
    fasp_solver_amg(&A, &b, &uh, &amgparam); 
    
    // Clean up memory
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&uh);
    
    return SUCCESS;
}
