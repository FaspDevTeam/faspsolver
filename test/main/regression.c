/**
 *		Regression test for FASP. 
 *
 *------------------------------------------------------
 *
 *		Created  by Chensong Zhang on 03/20/2010
 *      Modified by Chensong Zhang on 09/02/2011
 *      Modified by Chenosng Zhang on 12/03/2011
 *      Modified by Chensong Zhang on 03/20/2012
 *
 *------------------------------------------------------
 *
 */

/*! \file regression.c
 *  \brief Regression Testing Function for CSR Iterative Solvers
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

unsigned INT  ntest;    /**< number of tests all together */
unsigned INT  nfail;    /**< number of failed tests */

/**
 * \fn static void check_solu(dvector *x, dvector *sol, double tol)
 *
 * This function compares x and sol to a given tolerance tol. 
 */
static void check_solu(dvector *x, dvector *sol, double tol)
{
    double diff_u = fasp_dvec_maxdiff(x, sol);
    ntest++;
    
    if ( diff_u < tol ) {
        printf("Max diff %.4e smaller than tolerance................. [PASS]\n", diff_u);
    }
    else {
        nfail++;
        printf("### WARNING: Max diff %.4e BIGGER than tolerance..... [REQUIRES ATTENTION!!!]\n", diff_u);
    }	
}

/**
 * \fn int main (int argc, const char * argv[])
 * 
 * \brief This is the main function for regression test. 
 *
 * \author Chensong Zhang
 * \date   03/20/2009
 * 
 * Modified by Chensong Zhang on 09/02/2011 
 * Modified by Chensong Zhang on 12/03/2011
 * Modified by Chensong Zhang on 03/20/2012
 */

int main (int argc, const char * argv[]) 
{
    const INT      print_level = 1;    // how much information to print out
    const INT      num_prob    = 3;    // how many problems to be used
    const REAL     tolerance   = 1e-4; // tolerance for accepting the solution
    
    /* Local Variables */
    itsolver_param itparam;      // input parameters for iterative solvers
    AMG_param      amgparam; 	 // input parameters for AMG
    dCSRmat        A;            // coefficient matrix
    dvector        b, x, sol;    // rhs, numerical sol, exact sol 
    INT            indp;         // index for test problems
    
    time_t         lt  = time(NULL);
    
    printf("\n\n");
    printf("------------------------- Test starts at -------------------------\n");
    printf("%s",asctime(localtime(&lt))); // output starting local time
    printf("------------------------------------------------------------------\n");
    
    ntest = nfail = 0;
    
    /*******************************************/
    /* Step 1. Get matrix and right-hand side  */ 
    /*******************************************/
    for ( indp = 1; indp <= num_prob; indp++ ) {        
        
        printf("\n=====================================================\n");        
        printf("Test Problem Number %d ...\n", indp);	
        
        switch (indp) {
                
            case 1: //     - Problem 1. 10X10 5pt FD for Poisson
                
                printf("10X10 5-point finite difference for Poisson");	
                printf("\n=====================================================\n");        
                
                // Read A and b from two files in IJ format. 
                fasp_dcsrvec2_read("../data/csrmat_FD.dat", "../data/rhs_FD.dat", &A, &b);
                
                // Read ref. sol. from a non-indexed vec file.
                fasp_dvecind_read("../data/sol_FD.dat", &sol);
                
                break;
                
            case 2: //     - Problem 2. P1 FE for Poisson.
                
                printf("P1 finite element for Poisson");	
                printf("\n=====================================================\n");        
                
                // Read A and b from two files in IJ format. 
                fasp_dcsrvec2_read("../data/csrmat_FE.dat", "../data/rhs_FE.dat", &A, &b);
                
                // Read ref. sol. from an indexed vec file.
                fasp_dvecind_read("../data/sol_FE.dat", &sol);
                
                break;
                
            case 3: //     - Problem 3. MatrixMarket finite element analysis NOS7.
                // Finite difference approximation to diffusion equation with varying
                // diffusivity in a 3D unit cube with Dirichlet boundary conditions.
                
                printf("MatrixMarket finite element analysis NOS7");	
                printf("\n=====================================================\n");        
                
                // Read A in MatrixMarket SYM COO format. 
                fasp_dmtxsym_read("../data/nos7.mtx", &A);
                
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
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG V-cycle (Direct interpolation) with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG (direct interp) V-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG V-cycle (Standard interpolation) with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG (standard interp) V-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            amgparam.interpolation_type = INTERP_STD;
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG V-cycle (EM interpolation) with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG (energy-min interp) V-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            amgparam.interpolation_type = INTERP_ENG_MIN;
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-12;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* FMG V-cycle (Direct interpolation) with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("FMG (direct interp) V-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.print_level = print_level;
            fasp_solver_famg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG W-cycle with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG W-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.cycle_type  = W_CYCLE;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG AMLI-cycle with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG AMLI-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.cycle_type  = AMLI_CYCLE;
            amgparam.amli_degree = 3;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG Nonlinear AMLI-cycle with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG Nonlinear AMLI-cycle as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.cycle_type  = NL_AMLI_CYCLE;
            amgparam.amli_degree = 3;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG V-cycle with SGS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG V-cycle with SGS smoother as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.smoother    = SGS;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG V-cycle with L1_DIAG smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG V-cycle with L1_DIAG smoother as iterative solver ...\n");	
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 500;
            amgparam.tol         = 1e-10;
            amgparam.smoother    = L1_DIAG;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* AMG V-cycle with SOR smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG V-cycle with SOR smoother as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 20;
            amgparam.tol         = 1e-10;
            amgparam.smoother    = SOR;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* SA AMG V-cycle with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("SA AMG V-cycle with GS smoother as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 500;
            amgparam.tol         = 1e-10;
            amgparam.AMG_type    = SA_AMG;
            amgparam.smoother    = GS;
            amgparam.print_level = print_level;
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* UA AMG V-cycle with GS smoother as a solver */			
            printf("------------------------------------------------------------------\n");
            printf("UA AMG V-cycle with GS smoother as iterative solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgparam);
            amgparam.maxit       = 500;
            amgparam.tol         = 1e-10;
            amgparam.AMG_type    = UA_AMG;
            amgparam.smoother    = GS;
            amgparam.print_level = print_level;        
            fasp_solver_amg(&A, &b, &x, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 ) {
            /* CG */
            printf("------------------------------------------------------------------\n");
            printf("CG solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            itparam.precond_type  = PREC_NULL;	
            itparam.maxit         = 5000;
            itparam.tol           = 1e-12;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov(&A, &b, &x, &itparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 ) {
            /* BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("BiCGstab solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            itparam.precond_type  = PREC_NULL;	
            itparam.itsolver_type = SOLVER_BiCGstab;
            itparam.maxit         = 5000;
            itparam.tol           = 1e-12;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov(&A, &b, &x, &itparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 ) {
            /* BiCGstab in BSR */
            dBSRmat A_bsr = fasp_format_dcsr_dbsr (&A, 1);
            
            printf("------------------------------------------------------------------\n");
            printf("BiCGstab solver in BSR format ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            itparam.precond_type  = PREC_NULL;	
            itparam.itsolver_type = SOLVER_BiCGstab;
            itparam.maxit         = 500;
            itparam.tol           = 1e-12;
            itparam.print_level   = print_level;
            fasp_solver_dbsr_krylov(&A_bsr, &b, &x, &itparam);
            fasp_dbsr_free(&A_bsr);

            check_solu(&x, &sol, tolerance);            
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using diag(A) as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("Diagonal preconditioned CG solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            itparam.precond_type  = PREC_DIAG;
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_diag(&A, &b, &x, &itparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned CG solver ...\n");	
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned BiCGstab solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.itsolver_type = SOLVER_BiCGstab;
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for MinRes */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned MinRes solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.itsolver_type = SOLVER_MinRes;
            itparam.maxit         = 500;
            itparam.tol           = 1e-9;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned GMRes solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.itsolver_type = SOLVER_GMRES;
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for vGMRes */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned vGMRes solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.itsolver_type = SOLVER_VGMRES;
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for vFGMRes */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned vFGMRes solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.itsolver_type = SOLVER_VFGMRES;
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using classical AMG as preconditioner for GCG */
            printf("------------------------------------------------------------------\n");
            printf("AMG preconditioned GCG solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_amg_init(&amgparam);
            itparam.itsolver_type = SOLVER_GCG;
            itparam.maxit         = 500;
            itparam.tol           = 1e-10;
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            
            check_solu(&x, &sol, tolerance);
        }
        
#if FASP_USE_ILU
        if ( indp==1 || indp==2 || indp==3 ) {
            /* Using ILUk as preconditioner for CG */
            ILU_param      iluparam;
            printf("------------------------------------------------------------------\n");
            printf("ILUk preconditioned CG solver ...\n");	
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itparam);
            fasp_param_ilu_init(&iluparam);
            itparam.print_level   = print_level;
            fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &itparam, &iluparam);
            
            check_solu(&x, &sol, tolerance);
        }
#endif	
        
        /* clean up memory */
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b);
        fasp_dvec_free(&x);
        fasp_dvec_free(&sol);	
        
    } // end of for indp
    
	/* all done */
	lt = time(NULL);    
	printf("---------------------- All test finished at ----------------------\n");
    printf("%d tests finished: %d failed, %d succeeded!\n", ntest, nfail, ntest-nfail);
	printf("%s",asctime(localtime(&lt))); // output ending local time
	printf("------------------------------------------------------------------\n");
	
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
