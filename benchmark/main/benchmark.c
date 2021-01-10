/*! \file  benchmark.c
 *
 *  \brief Benchmark tests for iterative solvers
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2021 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 * 
 * \brief This is the main function for regression test. 
 *
 * \author Chensong Zhang
 * \date   12/28/2020
 */
int main (int argc, const char * argv[]) 
{
    INT   num_prob      = 235;  // how many problems to be used
    INT   prob_startID  = 1;    // start solve problem index
    const INT   print_level   = 1;    // how much information to print out
    const SHORT coarse_solver = 32;   // 0:  default coarsest-level solver
                                      // 31: SuperLU
                                      // 32: UMFPack
                                      // 33: MUMPS
                                      // 34: PARDISO
                                       
    const INT   maxit         = 500;  // maximal iteration number
    const REAL  tolerance     = 1e-6; // tolerance for accepting the solution
    
    /* Local Variables */
    INT        indp;    // index for test problems
    ITS_param  itspar;  // input parameters for iterative solvers
    AMG_param  amgpar;  // input parameters for AMG
    dCSRmat    A;       // coefficient matrix
    dvector    b, x;    // rhs and numerical solution
    int        status;  // iteration number or error message

	char logname[64];

	// Handle command line options
	int i = 1;
	while (i < argc) {
		if (!strcmp(argv[i], "-num_prob")) {
			num_prob = atoi(argv[i + 1]);
		 	printf("%s %d\n", argv[i], num_prob);
		}
		if (!strcmp(argv[i], "-prob_startID")) {
			prob_startID = atoi(argv[i + 1]);
		 	printf("%s %d\n", argv[i], prob_startID);
		}
		i += 2;
	}

    sprintf(logname, "./BenchmarkResults%d-%d.log", prob_startID, num_prob);  
    FILE *fp = fopen(logname, "w"); // save results

    REAL       Timer0, Timer1, TimeEachPoissonIter, Score; // performance timers
    time_t     lt = time(NULL);

    printf("------------------------- Test starts at -------------------------\n");
    printf("%s",asctime(localtime(&lt))); // output starting local time
    printf("------------------------------------------------------------------\n\n");

    printf("--------- Obtain local machine Poisson Performance Unit ----------\n");

    // Read A and b -- 5pt FD stencil for Poisson, 1M DoF, natural ordering
    //  fasp_dcoo_read("../data/Poisson/baseline.mat", &A);
    //  fasp_dvec_read("../data/Poisson/baseline.rhs", &b);
    fasp_dcoo_read("Poisson/baseline.mat", &A);
    fasp_dvec_read("Poisson/baseline.rhs", &b);
    // TODO Currently 2D, use 3D baseline in the future.

    // Call PCG
    fasp_dvec_alloc(b.row, &x);
    fasp_dvec_set(b.row, &x, 0.0);
    fasp_param_solver_init(&itspar); // default choice is CG
    itspar.precond_type  = PREC_NULL;
    itspar.maxit         = 1000;
    itspar.tol           = 1e-12;

    // Record wall time for 2D Poisson for benchmarking
    fasp_gettime(&Timer0);
    status = fasp_solver_dcsr_krylov(&A, &b, &x, &itspar);
    fasp_gettime(&Timer1);
    if ( status > 0 ) {
        TimeEachPoissonIter = (Timer1 - Timer0) / b.row / status;
        printf("It costs %d iterations or %.4e seconds.\n", status, Timer1-Timer0);
    }
    else {
        TimeEachPoissonIter = (Timer1 - Timer0) / b.row / itspar.maxit;
        printf("It costs %d iterations or %.4e seconds.\n", itspar.maxit, Timer1-Timer0);
    }

    printf("Local stencil performance spMV Unit (lMVU) of this computer is %.4e.\n",
           TimeEachPoissonIter);
    printf("------------------------------------------------------------------\n");

    fprintf(fp, "Tests performed at %s", asctime(localtime(&lt)));
    fprintf(fp, "lMVU of this computer is %.4e\n", TimeEachPoissonIter);
    
    // zhaoli 
    FILE *fpread = fopen("./MatrixCollectionName.txt", "r"); // read MatrixCollectionName.txt
    char buf[32];
    char dir[48];
    char name[96];
    int pageId = 1;
    indp = 0;
    while (!feof(fpread))
    {
        // pagek
        fscanf(fpread, "%s\n", buf);
        // printf("buf = %s\n", buf);

        // 跳过 page*
        sprintf(dir, "page%d", pageId);  
        // printf("dir = %s, buf=?dir = %d\n", dir, strcmp(buf, dir));  
        if (!strcmp(buf, dir)) {
            // printf("dir = %s\n", dir);
            pageId++;
            continue;
        }
        indp++;
        if (indp > num_prob) break;

        sprintf(dir, "./mtx/page%d", pageId - 1); 
        sprintf(name, "%s/%s/%s.mtx", dir, buf, buf); 
        // printf("probem id = %d, name = %s\n", indp, name);
        
        if(indp < prob_startID) continue;

#if 1
        /*****************************/
        /* Step 1. Read the systems  */
        /*****************************/

        printf("\n=====================================================\n");
        printf("Test Problem Number %d ...\n", indp);   

        printf("SuiteSparse Matrix Collection %s problem", buf);
        printf("\n=====================================================\n");

        // Read A in MatrixMarket SYM COO format.
        fasp_dmtxsym_read(name, &A);

        // printf("%d, %s, （n=%d）\n", indp, buf, A.row); 
        printf("\nMatrix size =  %d\n", A.row);

		// Generate an exact solution randomly
        dvector sol = fasp_dvec_create(A.row);
        fasp_dvec_rand(A.row, &sol);

        // Form the right-hand-side b = A*sol
        b = fasp_dvec_create(A.row);
        fasp_blas_dcsr_mxv(&A, sol.val, b.val);

        fprintf(fp, "==================================\n");
        fprintf(fp, "Test Problem %d: %s \n", indp, buf);
        fprintf(fp, "==================================\n");

        /*****************************/
        /* Step 2. Solve the systems */
        /*****************************/
        fasp_dvec_alloc(b.row, &x);  // allocate mem for numerical solution

        // AMG as solver
        if ( 0 < indp && indp <= num_prob ) {
            /* AMG V-cycle (Direct interpolation) with GS smoother as a solver */           
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG (direct interp) V-cycle as iterative solver ...");
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.print_level = print_level;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.maxit = maxit;
            amgpar.tol = tolerance;

            fasp_gettime(&Timer0);
            fasp_solver_amg(&A, &b, &x, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG V-cycle scores ............... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* AMG W-cycle with GS smoother as a solver */          
            printf("------------------------------------------------------------------\n");
            printf("Classical AMG W-cycle as iterative solver ...");
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgpar);
            amgpar.print_level = print_level;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.maxit = maxit;
            amgpar.tol = tolerance;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;

            fasp_gettime(&Timer0);
            fasp_solver_amg(&A, &b, &x, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG W-cycle scores ............... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* SA AMG V-cycle with GS smoother as a solver */           
            printf("------------------------------------------------------------------\n");
            printf("SAMG V-cycle with GS smoother as iterative solver ...");
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgpar);
            amgpar.print_level = print_level;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.maxit = maxit;
            amgpar.tol = tolerance;
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.aggregation_type = 2; // VMB
            amgpar.smoother = SMOOTHER_GS;

            fasp_gettime(&Timer0);
            fasp_solver_amg(&A, &b, &x, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG V-cycle scores ............... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* SA AMG W-cycle with GS smoother as a solver */
            printf("------------------------------------------------------------------\n");
            printf("SAMG W-cycle with GS smoother as iterative solver ...");

            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgpar);
            amgpar.print_level = print_level;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.maxit = maxit;
            amgpar.tol = tolerance;
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.aggregation_type = 2; // VMB
            amgpar.smoother = SMOOTHER_GS;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;

            fasp_gettime(&Timer0);
            fasp_solver_amg(&A, &b, &x, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG W-cycle scores ............... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* UA AMG V-cycle with GS smoother as a solver */           
            printf("------------------------------------------------------------------\n");
            printf("UAMG V-cycle with GS smoother as iterative solver ...");    
            
            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgpar);
            amgpar.print_level = print_level;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.maxit = maxit;
            amgpar.tol = tolerance;
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.quality_bound = 8.0;

            fasp_gettime(&Timer0);
            fasp_solver_amg(&A, &b, &x, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG V-cycle scores ............... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* UA AMG W-cycle with GS smoother as a solver */
            printf("------------------------------------------------------------------\n");
            printf("UAMG W-cycle with GS smoother as iterative solver ...");

            fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
            fasp_param_amg_init(&amgpar);
            amgpar.print_level = print_level;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.maxit = maxit;
            amgpar.tol = tolerance;
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.cycle_type = W_CYCLE;
            amgpar.quality_bound = 8.0;
            amgpar.max_levels = 8;

            fasp_gettime(&Timer0);
            fasp_solver_amg(&A, &b, &x, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG W-cycle scores ............... %.4e (lMVU).\n", Score);
        }

        // AMG preconditioned CG
        if ( 0 < indp && indp <= num_prob ) {
            /* Using classical AMG V-cycle as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("CAMG V-cycle preconditioned CG solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG V-cycle PCG scores ........... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using classical AMG W-cycle as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("CAMG W-cycle preconditioned CG solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.cycle_type = W_CYCLE;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            amgpar.max_levels = 8;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG W-cycle PCG scores ........... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using SA AMG V-cycle as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("SAMG V-cycle preconditioned CG solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG V-cycle PCG scores ........... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using SA AMG W-cycle as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("SAMG W-cycle preconditioned CG solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG W-cycle PCG scores ........... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using UA AMG V-cycle as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("UAMG V-cycle preconditioned CG solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            amgpar.quality_bound = 8.0;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG V-cycle PCG scores ........... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using UA AMG W-cycle as preconditioner for CG */
            printf("------------------------------------------------------------------\n");
            printf("UAMG W-cycle preconditioned CG solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            amgpar.quality_bound = 8.0;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG W-cycle PCG scores ........... %.4e (lMVU).\n", Score);
        }

        // AMG preconditioned BiCGstab
        if ( 0 < indp && indp <= num_prob ) {
            /* Using classical AMG V-cycle as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("CAMG V-cycle preconditioned BiCGstab solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_BiCGstab;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG V-cycle PBiCGstab scores ..... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using classical AMG W-cycle as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("CAMG W-cycle preconditioned BiCGstab solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_BiCGstab;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG W-cycle PBiCGstab scores ..... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using SA AMG V-cycle as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("SAMG V-cycle preconditioned BiCGstab solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_BiCGstab;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG V-cycle PBiCGstab scores ..... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using SA AMG W-cycle as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("SAMG W-cycle preconditioned BiCGstab solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_BiCGstab;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG W-cycle PBiCGstab scores ..... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using UA AMG V-cycle as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("UAMG V-cycle preconditioned BiCGstab solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.quality_bound = 8.0;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_BiCGstab;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG V-cycle PBiCGstab scores ..... %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using UA AMG W-cycle as preconditioner for BiCGstab */
            printf("------------------------------------------------------------------\n");
            printf("UAMG W-cycle preconditioned BiCGstab solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.quality_bound = 8.0;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_BiCGstab;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG W-cycle PBiCGstab scores ..... %.4e (lMVU).\n", Score);
        }

        // AMG preconditioned GMRES
        if ( 0 < indp && indp <= num_prob ) {
            /* Using classical AMG V-cycle as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("CAMG V-cycle preconditioned GMRes solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_GMRES;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG V-cycle PGMRes scores ........ %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using classical AMG W-cycle as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("CAMG W-cycle preconditioned GMRes solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_GMRES;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "CAMG W-cycle PGMRes scores ........ %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using SA AMG V-cycle as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("SAMG V-cycle preconditioned GMRes solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_GMRES;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG V-cycle PGMRes scores ........ %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using SA AMG W-cycle as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("SAMG W-cycle preconditioned GMRes solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = SA_AMG;
            amgpar.strong_coupled = 0.25; // cannot be too big
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_GMRES;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "SAMG W-cycle PGMRes scores ........ %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using UA AMG V-cycle as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("UAMG V-cycle preconditioned GMRes solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.max_levels = 8;
            amgpar.quality_bound = 8.0;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_GMRES;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG V-cycle PGMRes scores ........ %.4e (lMVU).\n", Score);
        }

        if ( 0 < indp && indp <= num_prob ) {
            /* Using UA AMG W-cycle as preconditioner for GMRes */
            printf("------------------------------------------------------------------\n");
            printf("UAMG W-cycle preconditioned GMRes solver ...");
            fasp_dvec_set(b.row,&x,0.0);
            fasp_param_solver_init(&itspar);
            fasp_param_amg_init(&amgpar);
            amgpar.AMG_type = UA_AMG;
            amgpar.smoother = SMOOTHER_GS;
            amgpar.coarse_solver = coarse_solver;
            amgpar.cycle_type = W_CYCLE;
            amgpar.max_levels = 8;
            amgpar.quality_bound = 8.0;
            itspar.print_level = print_level;
            itspar.maxit = maxit;
            itspar.tol = tolerance;
            itspar.itsolver_type = SOLVER_GMRES;

            fasp_gettime(&Timer0);
            fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
            fasp_gettime(&Timer1);

            Score = (Timer1 - Timer0) / b.row / TimeEachPoissonIter;
            printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            fprintf(fp, "UAMG W-cycle PGMRes scores ........ %.4e (lMVU).\n", Score);
        }

        /* clean up memory */
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b);
        fasp_dvec_free(&x);
#endif
    } // end of for indp


    fclose(fpread);

    /* all done */
    lt = time(NULL);    
    printf("---------------------- All test finished at ----------------------\n");
    printf("%s",asctime(localtime(&lt))); // output ending local time
    printf("------------------------------------------------------------------\n");

    fclose(fp);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
