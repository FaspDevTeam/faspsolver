/**
 *		Solver test for FASP. 
 *
 *------------------------------------------------------
 *
 *		Created  by Feiteng Huang on 06/12/2012
 *
 *------------------------------------------------------
 *
 */

/*! \file testmm.c
 *  \brief solver test for FASP with some matrix from Matrix-Market. 
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#define num_prob      12                  /**< how many problems to be used */
#define num_solvers    7                  /**< how many methods  to be used */
unsigned INT  ntest_diag[num_solvers];    /**< number of tests all together for diag preconditioner */
unsigned INT  nfail_diag[num_solvers];    /**< number of failed tests for diag preconditioner */
unsigned INT  ntest_ilu[num_solvers];     /**< number of tests all together for ilu preconditioner */
unsigned INT  nfail_ilu[num_solvers];     /**< number of failed tests for ilu preconditioner */
unsigned INT  ntest_amg[num_solvers];     /**< number of tests all together for amg preconditioner */
unsigned INT  nfail_amg[num_solvers];     /**< number of failed tests for amg preconditioner */
unsigned INT  ntest_amg_solver;           /**< number of tests all together for amg solver */
unsigned INT  nfail_amg_solver;           /**< number of failed tests for amg solver */

/**
 * \fn static void check_solu(dvector *x, dvector *sol, double tol)
 *
 * This function compares x and sol to a given tolerance tol. 
 */
static void check_solu(dvector *x, dvector *sol, double tol, int *nt, int *nf)
{
    double diff_u = fasp_dvec_maxdiff(x, sol);
    (*nt)++;
    
    if ( diff_u < tol ) {
        printf("Max diff %.4e smaller than tolerance................. [PASS]\n", diff_u);
    }
    else {
        (*nf)++;
        printf("### WARNING: Max diff %.4e BIGGER than tolerance..... [REQUIRES ATTENTION!!!]\n", diff_u);
    }	
}

/**
 * \fn int main (int argc, const char * argv[])
 * 
 * \brief This is the main function for solver test. 
 *
 * We pick a few test problems from the Matrix Market randomly. 
 * Some are symmetric and some are non-symmetric. We use these
 * test problems to check the package. 
 *
 * \author Feiteng Huang 
 * \date   06/12/2012
 */
int main (int argc, const char * argv[]) 
{
    const INT      print_level = 1;    // how much information to print out
    const REAL     tolerance   = 1e-4; // tolerance for accepting the solution
    
    /* Local Variables */
    itsolver_param itparam;      // input parameters for iterative solvers
    ILU_param      iluparam; 	 // input parameters for AMG
    AMG_param      amgparam; 	 // input parameters for AMG
    dCSRmat        A;            // coefficient matrix
    dvector        b, x, sol;    // rhs, numerical sol, exact sol 
    INT            indp;         // index for test problems
    INT            indm;         // index for test methods
    const char     solvers[num_solvers][128] = {"CG", "BiCGstab", "MinRes", "GMRES", "VGMRES", "VFGMRES", "GCG"};
    
    time_t         lt  = time(NULL);
    
    printf("\n\n");
    printf("------------------------- Test starts at -------------------------\n");
    printf("%s",asctime(localtime(&lt))); // output starting local time
    printf("------------------------------------------------------------------\n");
    
    memset(ntest_diag, 0x0, num_solvers*sizeof(int));
    memset(nfail_diag, 0x0, num_solvers*sizeof(int));
    memset(ntest_ilu,  0x0, num_solvers*sizeof(int));
    memset(nfail_ilu,  0x0, num_solvers*sizeof(int));
    memset(ntest_amg,  0x0, num_solvers*sizeof(int));
    memset(nfail_amg,  0x0, num_solvers*sizeof(int));
    ntest_amg_solver = 0;
    nfail_amg_solver = 0;
    
    /*******************************************/
    /* Step 1. Get matrix and right-hand side  */ 
    /*******************************************/
    for ( indp = 1; indp <= num_prob; indp++ ) {        
        
        printf("\n==================================================================\n");        
        printf("Test Problem Number %d ...\n", indp);	
        
        switch (indp) {
                
            case 1: //     - Problem 1. MatrixMarket Driven cavity E05R0500.
                // driven cavity, 5x5 elements, Re=500
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/e05r0500.mtx", &A);
                
                printf("MatrixMarket Driven cavity E05R0500\n");	
                printf("||  Condition Number:      4.8e+6   ||\n");
                printf("||        Unsymmetric               ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 2: //     - Problem 2. MatrixMarket Finite element analysis of cylindrical shells.
                // Cylindrical shell, uniform 30x30 quadrilateral mesh, 
                // stabilized MITC4 elements, R/t=100
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtxsym_read("../data/s2rmq4m1.mtx", &A);
                
                printf("MatrixMarket Finite element analysis of cylindrical shells\n");	
                printf("||  Condition Number:     1.15e+8   ||\n");
                printf("||   Symmetric positive definite    ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 3: //     - Problem 3. MatrixMarket Oil reservoir simulation - generated problems.
                // oil reservoir simulation for 21x21x5 full grid
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/orsreg_1.mtx", &A);
                
                printf("MatrixMarket Oil reservoir simulation - generated problems\n");	
                printf("||  Condition Number:        1e+2   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 4: //     - Problem 4. MatrixMarket Enhanced oil recovery.
                // 3D steam model of oil res. -5x5x6 -4 DOF
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/steam2.mtx", &A);
                
                printf("MatrixMarket Enhanced oil recovery\n");	
                printf("||  Condition Number:      3.5e+6   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 5: //     - Problem 5. MatrixMarket BCS Structural Engineering Matrices.
                // U.S. Army Corps of Engineers dam
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtxsym_read("../data/bcsstk16.mtx", &A);
                
                printf("MatrixMarket BCS Structural Engineering Matrices\n");	
                printf("||  Condition Number:          65   ||\n");
                printf("||   Symmetric positive definite    ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 6: //     - Problem 6. MatrixMarket Circuit physics modeling.
                // Computer random simulation of a circuit physics model
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/jpwh_991.mtx", &A);
                
                printf("MatrixMarket Circuit physics modeling\n");	
                printf("||  Condition Number:      7.3e+2   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 7: //     - Problem 7. MatrixMarket Simulation of computer systems.
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/gre__115.mtx", &A);
                
                printf("MatrixMarket Simulation of computer systems\n");	
                printf("||  Condition Number:      1.5e+2   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 8: //     - Problem 8. MatrixMarket Computer component design.
                // 32-bit adder
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/add32.mtx", &A);
                
                printf("MatrixMarket Computer component design\n");	
                printf("||  Condition Number:     2.14e+2   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 9: //     - Problem 9. MatrixMarket Fluid flow modeling.
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/lns__131.mtx", &A);
                
                printf("MatrixMarket Fluid flow modeling\n");	
                printf("||  Condition Number:    1.15e+15   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 10: //     - Problem 10. MatrixMarket Nuclear reactor models.
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/nnc1374.mtx", &A);
                
                printf("MatrixMarket Nuclear reactor models\n");	
                printf("||  Condition Number:        1e+2   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 11: //     - Problem 11. MatrixMarket Oil reservoir simulation challenge matrics.
                // Black oil simulation, shale barriers(NX = NY = NZ = 10, NC = 1)
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtxsym_read("../data/sherman1.mtx", &A);
                
                printf("MatrixMarket Oil reservoir simulation challenge matrics\n");	
                printf("||  Condition Number:      2.3e+4   ||\n");
                printf("||            Symmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
                
            case 12: //     - Problem 12. MatrixMarket Petroleum engineering.
                
                // Read A in MatrixMarket COO format. 
                fasp_dmtx_read("../data/watt__1.mtx", &A);
                
                printf("MatrixMarket Petroleum engineering\n");	
                printf("||  Condition Number:     5.38e+9   ||\n");
                printf("||          Unsymmetric             ||\n");
                printf("|| row:%5d, col:%5d, nnz:%6d ||\n", A.row, A.col, A.nnz);
                printf("==================================================================\n");        
                
                // Generate an exact solution randomly
                sol = fasp_dvec_create(A.row);
                fasp_dvec_rand(A.row, &sol);
                
                // Form the right-hand-side b = A*sol
                b = fasp_dvec_create(A.row);
                fasp_blas_dcsr_mxv(&A, sol.val, b.val);
                
                break;
        }
        
        /************************************/
        /* Step 2. Check matrix properties  */
        /************************************/
        fasp_check_symm(&A);     // check symmetry
        fasp_check_diagpos(&A);  // check sign of diagonal entries
        fasp_check_diagdom(&A);  // check diagonal dominance
        
        /*****************************/	
        /* Step 3. Solve the system  */ 
        /*****************************/
        fasp_dvec_alloc(b.row, &x);  // allocate mem for numerical solution
        
        if (1) {
            /* Using diagonal preconditioner for Krylov methods */
            printf("------------------------------------------------------------------\n");
            printf("Diagonal preconditioned krylov solver ...\n");	
            
            fasp_param_solver_init(&itparam);
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            iluparam.ILU_type     = ILUk;
            for (indm = 0; indm<7; indm++) {
                fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
                itparam.itsolver_type = indm+1;
                fasp_solver_dcsr_krylov_diag(&A, &b, &x, &itparam);
                
                check_solu(&x, &sol, tolerance, &(ntest_diag[indm]), &(nfail_diag[indm]));
                printf("\n");
            }
        }
        
        if (1) {
            /* Using ILUk as preconditioner for Krylov methods */
            printf("------------------------------------------------------------------\n");
            printf("ILUk preconditioned krylov solver ...\n");	
            
            fasp_param_solver_init(&itparam);
            fasp_param_ilu_init(&iluparam);
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            iluparam.ILU_type     = ILUk;
            for (indm = 0; indm<7; indm++) {
                fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
                itparam.itsolver_type = indm+1;
                fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &itparam, &iluparam);
                
                check_solu(&x, &sol, tolerance, &(ntest_ilu[indm]), &(nfail_ilu[indm]));
                printf("\n");
            }
        }
        
        if (1) {
            /* Using classical AMG as preconditioner for Krylov methods */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned krylov solver ...\n");	
            
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            for (indm = 0; indm<7; indm++) {
                fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
                itparam.itsolver_type = indm+1;
                fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
                
                check_solu(&x, &sol, tolerance, &(ntest_amg[indm]), &(nfail_amg[indm]));
                printf("\n");
            }
        }
        
        if (1) {
            /* Using classical AMG as a solver */
            amgparam.maxit        = 20;
            amgparam.tol          = 1e-10;
            amgparam.print_level  = print_level;
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            printf("Calling AMG solver ...\n");
            fasp_solver_amg(&A, &b, &x,&amgparam);
            
            check_solu(&x, &sol, tolerance, &ntest_amg_solver, &nfail_amg_solver);
            printf("\n");
        }
        
        /* clean up memory */
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b);
        fasp_dvec_free(&x);
        fasp_dvec_free(&sol);	
        
    } // end of for indp
    
	/* all done */
	lt = time(NULL);    
	printf("\n---------------------- All test finished at ----------------------\n");
	printf("%s",asctime(localtime(&lt))); // output ending local time
    printf("------------------------------------------------------------------\n\n");
    
	printf("==================================================================\n");
    printf("=======================Diagonal preconditioner====================\n");
    for (indm=0; indm<7; indm++) {
        printf("Solver %10s:  %d tests finished: %d failed, %d succeeded!\n", 
               solvers[indm], ntest_diag[indm], nfail_diag[indm], ntest_diag[indm]-nfail_diag[indm]);
    }
    printf("=========================ILU preconditioner======================\n");
    for (indm=0; indm<7; indm++) {
        printf("Solver %10s:  %d tests finished: %d failed, %d succeeded!\n", 
               solvers[indm], ntest_ilu[indm], nfail_ilu[indm], ntest_ilu[indm]-nfail_ilu[indm]);
    }
    printf("=========================AMG preconditioner======================\n");
    for (indm=0; indm<7; indm++) {
        printf("Solver %10s:  %d tests finished: %d failed, %d succeeded!\n", 
               solvers[indm], ntest_amg[indm], nfail_amg[indm], ntest_amg[indm]-nfail_amg[indm]);
    }
	printf("==================================================================\n");
    printf("Solver %10s:  %d tests finished: %d failed, %d succeeded!\n", "AMG", 
           ntest_amg_solver, nfail_amg_solver, ntest_amg_solver-nfail_amg_solver);
	printf("==================================================================\n");
	
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
