/*! \file testfem.c
 *  \brief The main test function for FASP FEM assembling.
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "misc.h"
#include "mesh.h"
#include "poisson_fem.h"

/* Test functions f and u for the Poisson's equation */
#include "testfct_poisson.inl"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for testing FASP FEM assembling.
 *
 * \author Chensong Zhang
 * \date   09/11/2011
 *
 * Modified by Chensong Zhang and Feiteng Huang on 03/07/2012
 * Modified by Feiteng Huang on 04/01/2012: l2 error
 * Modified by Feiteng Huang on 04/05/2012: refine & assemble
 * Modified by Chensong Zhang on 10/15/2012: add different solvers
 */
int main (int argc, const char * argv[])
{
    // Set default values
    int status = FASP_SUCCESS;
    int print_usage;

    FEM_param fempar;// parameter for testfem
    FEM_param_init(&fempar);
    print_usage = FEM_param_set(argc, argv, &fempar);

    if ( print_usage ) {
        printf("\nUsage: %s [<ofemparions>]\n", argv[0]);
        printf("  -output <val>    : mesh output switch [default: 0]\n");
        printf("  -meshin <val>    : input mesh  [default: ../data/mesh.dat]\n");
        printf("  -meshout <val>   : output mesh [default: ../data/mesh_?.dat]\n");
        printf("  -refine <val>    : refine level [default: 8]\n");
        printf("  -assemble <val>  : assemble ofemparion [default: ab]\n");
        printf("                     ab  |  assemble the mat & rhs;\n");
        printf("                      a  |  assemble the mat;\n");
        printf("                      b  |  assemble the rhs;\n");
        printf("  -quad_rhs <val>  : quad points for rhs [default: 3]\n");
        printf("  -quad_mat <val>  : quad points for mat [default: 1]\n");
        printf("  -help            : print this help message\n\n");
        exit(status);
    }

    // Local variables
    Mesh       mesh;       // Mesh info
    Mesh_aux   mesh_aux;   // Auxiliary mesh info
    dCSRmat    A;          // Stiffness matrix
    dvector    b;          // Right-hand side
    dvector    uh;         // Solution including boundary
    ivector    dof;        // Degree of freedom
    int        i;          // Loop index

    // Step 1: read initial mesh info
    mesh_init (&mesh, "../data/mesh.dat");

    // Step 1.5: read or build auxiliary mesh info
    // If there is already mesh_aux data available, you can use the following fct to init it:
    //     mesh_aux_init (&mesh, &mesh_aux, "mesh.aux");
    // otherwise, just build the auxiliary information for the mesh
    mesh_aux_build(&mesh, &mesh_aux);

    // Step 2: refine mesh
    clock_t mesh_refine_s = clock();
    for ( i=0; i < fempar.refine_lvl; ++i ) mesh_refine(&mesh, &mesh_aux);
    clock_t mesh_refine_e = clock();
    double mesh_refine_time = (double)(mesh_refine_e - mesh_refine_s)/(double)(CLOCKS_PER_SEC);
    printf("Mesh refinement costs... %8.4f seconds\n", mesh_refine_time);

    // Step 3: assemble the linear system and right-hand side
    clock_t assemble_s = clock();
    setup_poisson(&A, &b, &mesh, &mesh_aux, &fempar, &uh, &dof);
    clock_t assemble_e = clock();
    double assemble_time = (double)(assemble_e - assemble_s)/(double)(CLOCKS_PER_SEC);
    printf("Assembling costs........ %8.4f seconds\n", assemble_time);

    // Step 3.5: clean up auxiliary mesh info
    mesh_aux_free(&mesh_aux);
#if 0    
// test by chunsheng
     char matfile[80];
     sprintf(matfile, "out/fem_mat_%03d.dat",fempar.refine_lvl); 
     fasp_dcsr_write_coo(matfile, &A);
     sprintf(matfile, "out/fem_rhs_%03d.dat",fempar.refine_lvl); 
     fasp_dvec_write (matfile,&b);
#endif

    // Step 4: solve the linear system with AMG
    {
        dvector x; // just on DOF, not including boundary

        input_param     inparam;  // parameters from input files
        itsolver_param  itparam;  // parameters for itsolver
        AMG_param       amgparam; // parameters for AMG
        ILU_param       iluparam; // parameters for ILU
        Schwarz_param   swzparam; // parameters for Schwarz method

        // Read input parameters from a disk file
        fasp_param_input_init(&inparam);
        fasp_param_init(&inparam,&itparam,&amgparam,&iluparam,&swzparam);

		amgparam.ILU_type             = ILUtp;
		amgparam.ILU_levels           = 2;
#if 0
// test by chunsheng
	itparam.itsolver_type= SOLVER_BiCGstab;
	itparam.itsolver_type= SOLVER_VBiCGstab;
        itparam.precond_type =  PREC_NULL;
	//itparam.precond_type =  PREC_AMG;
	inparam.print_level = 3;
	inparam.solver_type = SOLVER_VBiCGstab;
//	inparam.precond_type = PREC_AMG;
	inparam.precond_type = PREC_NULL;
        itparam.maxit = 1000;
// test by chunsheng
#endif

	const int print_level   = inparam.print_level;
        const int solver_type   = inparam.solver_type;
        const int precond_type  = inparam.precond_type;

        // Set initial guess
        fasp_dvec_alloc(A.row, &x);
        fasp_dvec_set(A.row,&x,0.0);
          
        // Preconditioned Krylov methods
        if ( solver_type >= 1 && solver_type <= 20) {

            // Using no preconditioner for Krylov iterative methods
            if ( precond_type == PREC_NULL ) {
                status = fasp_solver_dcsr_krylov(&A, &b, &x, &itparam);
            //iter = fasp_solver_dcsr_pvbcgs(A, b, x, pc, tol, MaxIt, stop_type, prtlvl);
            }

            // Using diag(A) as preconditioner for Krylov iterative methods
            else if ( precond_type == PREC_DIAG ) {
                status = fasp_solver_dcsr_krylov_diag(&A, &b, &x, &itparam);
            }

            // Using AMG as preconditioner for Krylov iterative methods
            else if ( precond_type == PREC_AMG || precond_type == PREC_FMG ) {
                if ( print_level > PRINT_NONE ) fasp_param_amg_print(&amgparam);
                status = fasp_solver_dcsr_krylov_amg(&A, &b, &x, &itparam, &amgparam);
            }

            // Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
            else if ( precond_type == PREC_ILU ) {
                if ( print_level > PRINT_NONE ) fasp_param_ilu_print(&iluparam);
                status = fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &itparam, &iluparam);
            }

            // Using Schwarz as preconditioner for Krylov iterative methods
            else if ( precond_type == PREC_SCHWARZ ){
                if ( print_level > PRINT_NONE ) fasp_param_schwarz_print(&swzparam);
                status = fasp_solver_dcsr_krylov_schwarz(&A, &b, &x, &itparam, &swzparam);
            }

            else {
                printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);
                exit(ERROR_SOLVER_PRECTYPE);
            }

        }

        // AMG as the iterative solver
        else if ( solver_type == SOLVER_AMG ) {
            if ( print_level > PRINT_NONE ) fasp_param_amg_print(&amgparam);
            fasp_solver_amg(&A, &b, &x, &amgparam);
        }

        // Full AMG as the iterative solver
        else if ( solver_type == SOLVER_FMG ) {
            if ( print_level > PRINT_NONE ) fasp_param_amg_print(&amgparam);
            fasp_solver_famg(&A, &b, &x, &amgparam);
        }

        else {
            printf("### ERROR: Wrong solver type %d!!!\n", solver_type);
            status = ERROR_SOLVER_TYPE;
        }

        for ( i = 0; i < dof.row; ++i) uh.val[dof.val[i]] = x.val[i];

        fasp_dvec_free(&x);
    }

    // Step 5: estimate L2 error
    double l2error = get_l2_error_poisson(&(mesh.node), &(mesh.elem), &uh, fempar.num_qp_rhs);

    printf("\n==============================================================\n");
    printf("L2 error of FEM is %g\n", l2error);
    printf("==============================================================\n\n");

    // clean up memory
    mesh_free(&mesh);
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&uh);
    fasp_ivec_free(&dof);

    return FASP_SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
