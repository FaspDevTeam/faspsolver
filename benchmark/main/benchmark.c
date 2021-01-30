/*! \file  benchmark.c
 *
 *  \brief Benchmark tests for iterative solvers
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2020--Present by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "benchmark.h"

static double ComputeLMVUFromBaseline(char *workdir, Baseline bl);

/**
 * \fn int main (int argc, const char * argv[])
 * 
 * \brief This is the main function for benchmark test.
 *
 * \author Chensong Zhang, Li Zhao
 * \date   1/10/2021
 */
int main (int argc, const char * argv[])
{
    char        ini_file[20]  = "input.dat";
    char        blp_dir[20]   = "baseline";
    char        mat_dir[20]   = "mtx";
    char        alg_dir[20]   = "algorithm";
    
    INT         startID       = 1;     // start solve problem index
    INT         endID         = 1;     // end solve problem index
    INT         MinProbSize   = 10000; // Filter small matrices
    const INT   maxit         = 500;   // maximal iteration number
    const REAL  tolerance     = 1e-6;  // tolerance for accepting the solution
    const INT   print_level   = 2;     // how much information to print out
    const SHORT coarse_solver = 34;    // 0:  default coarsest-level solver
                                       // 31: SuperLU, 32: UMFPack, 33: MUMPS
                                       // 34: PARDISO

    /* Local Variables */
    Baseline     bl;
    Problem      pb;
    Algorithm    ag;
    INT          indp;   // index for test problems
    input_param  inipar; // parameters from input files
    ITS_param    itspar; // parameters for itsolver
    AMG_param    amgpar; // parameters for AMG
    ILU_param    ilupar; // parameters for ILU
    SWZ_param    swzpar; // parameters for Schwarz method

    int          status; // iteration number or error message
    dCSRmat      A;      // coefficient matrix
    dvector      b, x;   // rhs and numerical solution
    int          solver_type, precond_type;
    char         logname[64];
    REAL         Timer0, Timer1, lMVU, Score; // performance timers
    time_t       lt = time(NULL);

    printf("------------------------- Test starts at -------------------------\n");
    printf("%s",asctime(localtime(&lt))); // output starting local time
    printf("------------------------------------------------------------------\n\n");

    // Handle command line options
    int i = 1;
    while (i < argc) {
        if (!strcmp(argv[i], "-startID")) startID = atoi(argv[i + 1]);
        if (!strcmp(argv[i], "-endID"))   endID = atoi(argv[i + 1]);
        if (!strcmp(argv[i], "-f"))       strcpy(ini_file, argv[i + 1]);
        if (!strcmp(argv[i], "-mps"))     MinProbSize = atoi(argv[i + 1]);
        i += 2;
    }
    printf("Start ID of test problem = %d\n", startID);
    printf("Ending ID of test problem = %d\n", endID);
    printf("Test problem size larger than %d\n", MinProbSize);

    //-----------------------------//
    //     Read input.dat file     //
    //-----------------------------//
    printf("Reading parameters from input file: %s\n\n", ini_file);
    status = ReadInputFile(ini_file, &bl, &pb, &ag);
    if (status < 0) return ERROR_READ_FILE;

    //--------------------------------------------//
    //     Compute lMVU from Baseline problem     //
    //--------------------------------------------//
    lMVU = ComputeLMVUFromBaseline(blp_dir, bl);

    sprintf(logname, "./BenchmarkResults-%d-%d.log", startID, endID);
    FILE *fp = fopen(logname, "w"); // save results
    fprintf(fp, "Tests performed at %s", asctime(localtime(&lt)));
    fprintf(fp, "lMVU of this computer is %.4e\n", lMVU);

    FILE *fpCheck = NULL;
    char matrix_file_name[128];
    char algorithm_file_name[128];

    // Loop through the problem to be tested
    for (indp = 1; indp < pb->num+1; indp++)
    {
        sprintf(matrix_file_name, "%s/%s/%s.mtx", mat_dir, pb->prob[indp-1], pb->prob[indp-1]);

        if (indp < startID || indp > endID) continue;
        printf("\n=====================================================\n");
        printf("Test Problem ID:      %d\n", indp);
        printf("Name of Problem:      %s\n", pb->prob[indp-1]);

        // Check if the file exists
        fpCheck = fopen(matrix_file_name, "r");
        if (!fpCheck)
        {
            printf("### WARNING: The file %s does not exist! Skip...\n", matrix_file_name);
            continue;
        }
        fclose(fpCheck);

        /*****************************/
        /* Step 1. Read the systems  */
        /*****************************/
        // Read A in MatrixMarket SYM COO format.
        fasp_dmtxsym_read(matrix_file_name, &A);

        // Filter small matrix
        if (A.row < MinProbSize) {
            printf("### WARNING: Skip matrices of size less than %d!\n", MinProbSize);
            continue;
        }

        // Generate an exact solution randomly
        dvector sol = fasp_dvec_create(A.row);
        fasp_dvec_rand(A.row, &sol);

        // Form the right-hand-side b = A*sol
        b = fasp_dvec_create(A.row);
        fasp_blas_dcsr_mxv(&A, sol.val, b.val);

        printf("Problem Size:         %d\n", A.row);
        printf("=====================================================\n");

        fprintf(fp, "======================================================================\n");
        fprintf(fp, "Test Problem %3d: %12s, matrix size is %6d\n", indp, pb->prob[indp-1], A.row);
        fprintf(fp, "======================================================================\n");

        /*****************************/
        /* Step 2. Solve the systems */
        /*****************************/
        fasp_dvec_alloc(b.row, &x);  // allocate mem for numerical solution

        // Loop the algorithm file and call the corresponding algorithm
        for (i = 0; i < ag->num; i++)
        {
            printf("\n\n-----------------------------------------------\n");
            sprintf(algorithm_file_name, "%s/%s.dat", alg_dir, ag->para[i]);
            printf("Algorithm ID: %d, param filename: %s", i+1, algorithm_file_name);

            // Set solver parameters
            fasp_param_input(algorithm_file_name, &inipar);
            fasp_param_init(&inipar, &itspar, &amgpar, &ilupar, &swzpar);

            // Set local parameters
            solver_type   = inipar.solver_type;
            precond_type  = inipar.precond_type;

            // Set initial guess
            fasp_dvec_set(A.row, &x, 0.0);

            // Preconditioned Krylov methods
            if (solver_type >= 1 && solver_type <= 20) {

                // Using no preconditioner for Krylov iterative methods
                if (precond_type == PREC_NULL) {
                    fasp_gettime(&Timer0);
                    status = fasp_solver_dcsr_krylov(&A, &b, &x, &itspar);
                    fasp_gettime(&Timer1);
                }

                // Using diag(A) as preconditioner for Krylov iterative methods
                else if (precond_type == PREC_DIAG) {
                    fasp_gettime(&Timer0);
                    status = fasp_solver_dcsr_krylov_diag(&A, &b, &x, &itspar);
                    fasp_gettime(&Timer1);
                }

                // Using AMG as preconditioner for Krylov iterative methods
                else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
                    if (print_level>PRINT_NONE) fasp_param_amg_print(&amgpar);
                    fasp_gettime(&Timer0);
                    status = fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itspar, &amgpar);
                    fasp_gettime(&Timer1);
                }

                // Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
                else if (precond_type == PREC_ILU) {
                    if (print_level>PRINT_NONE) fasp_param_ilu_print(&ilupar);
                    fasp_gettime(&Timer0);
                    status = fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &itspar, &ilupar);
                    fasp_gettime(&Timer1);
                }

                // Using Schwarz as preconditioner for Krylov iterative methods
                else if (precond_type == PREC_SCHWARZ){
                    if (print_level>PRINT_NONE) fasp_param_swz_print(&swzpar);
                    fasp_gettime(&Timer0);
                    status = fasp_solver_dcsr_krylov_swz(&A, &b, &x, &itspar, &swzpar);
                    fasp_gettime(&Timer1);
                }

                else {
                    printf("### ERROR: Unknown preconditioner type %d!!!\n", precond_type);
                    status = ERROR_SOLVER_PRECTYPE;
                }

            }

            // AMG as the iterative solver
            else if (solver_type == SOLVER_AMG) {
                if (print_level>PRINT_NONE) fasp_param_amg_print(&amgpar);
                fasp_gettime(&Timer0);
                fasp_solver_amg(&A, &b, &x, &amgpar);
                fasp_gettime(&Timer1);
            }

            // Full AMG as the iterative solver
            else if (solver_type == SOLVER_FMG) {
                if (print_level>PRINT_NONE) fasp_param_amg_print(&amgpar);
                fasp_gettime(&Timer0);
                fasp_solver_famg(&A, &b, &x, &amgpar);
                fasp_gettime(&Timer1);
            }

#if WITH_SuperLU // use SuperLU directly
            else if (solver_type == SOLVER_SUPERLU) {
                status = fasp_solver_superlu(&A, &b, &x, print_level);
            }
#endif

#if WITH_UMFPACK // use UMFPACK directly
            else if (solver_type == SOLVER_UMFPACK) {
                    
                dCSRmat A_tran;
                fasp_dcsr_trans(&A, &A_tran);
                fasp_dcsr_sort(&A_tran);
                fasp_dcsr_cp(&A_tran, &A);
                fasp_dcsr_free(&A_tran);
                    
                void *Numeric = fasp_umfpack_factorize(&A, print_level);
                status = fasp_umfpack_solve(&A, &b, &x, Numeric, print_level);
                fasp_umfpack_free_numeric(Numeric);
                        
            }
#endif

            else {
                printf("### ERROR: Unknown solver type %d!!!\n", solver_type);
                status = ERROR_SOLVER_TYPE;
            }

            if (status<0) {
                printf("\n### ERROR: Solver failed! Exit status = %d.\n\n", status);
            }

            // Calculate Score
            Score = (Timer1 - Timer0) / b.row / lMVU;
            // printf("Total cost is \033[31;43m%.4e\033[0m (lMVU).\n", Score);
            printf("Total cost is %.4e (lMVU).\n", Score);
            fprintf(fp, "%20s scores ............... %.4e (lMVU), iteration is %5d, status is %5d.\n", ag->para[i], \
                    Score, status > 0 ? status: inipar.itsolver_maxit, status);
        }

        // free
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b);
        fasp_dvec_free(&x);
        printf("\n");

    } // end of for indp

    /* all done */
    lt = time(NULL);
    printf("---------------------- All test finished at ----------------------\n");
    printf("%s",asctime(localtime(&lt))); // output ending local time
    printf("------------------------------------------------------------------\n");

    // free
    fclose(fp);
    FreeBaseline(bl);
    FreeProblem(pb);
    FreeAlgorithm(ag);

    // return
    return 1;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

static double ComputeLMVUFromBaseline(char *workdir, Baseline bl)
{
    ITS_param  itspar;  // input parameters for iterative solvers
    AMG_param  amgpar;  // input parameters for AMG
    dCSRmat    A;       // coefficient matrix
    dvector    b, x;    // rhs and numerical solution
    int        status;  // iteration number or error message
    double     Timer0, Timer1;

    int baselineID, baselineNum, callnum, i;
    int baselineTotNum = 0;
    double TotlMVU = 0, lMVU;
    char baseAfile[64];
    char basebfile[64];

    FILE *fpCheckBaseA = NULL, *fpCheckBaseb = NULL;
    double t1, t2, baselineTime = 0;

    printf("--------- Obtain local machine Baseline Performance Unit ----------\n");
    fasp_gettime(&t1);
    for (i = 0; i < bl->num; i++)
    {
        sprintf(baseAfile, "%s/%s/%s.mat", workdir, bl->prob[i], bl->prob[i]);
        sprintf(basebfile, "%s/%s/%s.rhs", workdir, bl->prob[i], bl->prob[i]);
        printf("Baseline matrix filename          = %s\n", baseAfile);
        printf("Baseline right-hand side filename = %s\n", basebfile);
        printf("Baseline test calling times       = %d\n", bl->callnums[i]);

        // Check if the file exists
        fpCheckBaseA = fopen(baseAfile, "r");
        fpCheckBaseb = fopen(basebfile, "r");
        if (!fpCheckBaseA || !fpCheckBaseb)
        {
            if (!fpCheckBaseA)
                printf("Warning: The file %s does not exist!!!\n", baseAfile);

            if (!fpCheckBaseb)
                printf("Warning: The file %s does not exist!!!\n", basebfile);

            printf("\n");
            continue;
        }
        fclose(fpCheckBaseA);
        fclose(fpCheckBaseb);

        // Read A and b -- 5pt FD stencil for Poisson, 1M DoF, natural ordering
        fasp_dcoo_read(baseAfile, &A);
        fasp_dvec_read(basebfile, &b);
        // TODO Currently 2D, use 3D baseline in the future.

        // Call PCG
        fasp_dvec_alloc(b.row, &x);
        fasp_param_solver_init(&itspar); // default choice is CG
        itspar.precond_type  = PREC_NULL;
        itspar.maxit         = 1000;
        itspar.tol           = 1e-20;

        printf("Matrix size is %d\n", b.row);
        for (callnum = 0; callnum < bl->callnums[i]; callnum++)
        {
            printf("Number of calls: %d\n", callnum + 1);
            fasp_dvec_set(b.row, &x, 0.0);
            // Record wall time for 2D Poisson for benchmarking
            fasp_gettime(&Timer0);
            status = fasp_solver_dcsr_krylov(&A, &b, &x, &itspar);
            fasp_gettime(&Timer1);
            if ( status > 0 ) {
                lMVU = (Timer1 - Timer0) / b.row / status;
                printf("\tIt costs %d iterations or %.4e seconds.\n", status, Timer1-Timer0);
            }
            else {
                lMVU = (Timer1 - Timer0) / b.row / itspar.maxit;
                printf("\tIt costs %d iterations or %.4e seconds.\n", itspar.maxit, Timer1-Timer0);
            }
            printf("\tLocal stencil spMV Unit (lMVU) of this computer is %.4e.\n", lMVU);
            TotlMVU += lMVU;
        }
        baselineTotNum += bl->callnums[i];
        // free
        fasp_dcsr_free(&A);
        fasp_dvec_free(&b);
        fasp_dvec_free(&x);
        printf("\n");
    }
    fasp_gettime(&t2);
    baselineTime = t2 - t1;

    lMVU = TotlMVU / baselineTotNum;
    printf("Baseline test called %d, CPU time: %.2fs.\n", baselineTotNum, baselineTime);
    printf("Average local stencil spMV Unit (lMVU) of this computer is %.4e.\n", lMVU);
    printf("------------------------------------------------------------------\n");

    // return
    return lMVU;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
