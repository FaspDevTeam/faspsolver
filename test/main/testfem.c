/**
 *    Testing the FASP FEM assembling routines. 
 *
 *------------------------------------------------------
 *
 */

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
 * Modified by Feiteng Huang on 04/01/2012, for l2 error
 * Modified by Feiteng Huang on 04/05/2012, for refine & assemble
 */
int main (int argc, const char * argv[]) 
{
    // Set default values
    int status = SUCCESS;
	int print_usage;
	
    FEM_param fempar;// parameter for testfem
	FEM_param_init(&fempar);
	print_usage = FEM_param_set(argc, argv, &fempar);
    
	if (print_usage) {
        printf("\nUsage: %s [<ofemparions>]\n", argv[0]);
        printf("  -meshin <val>    : input mesh file  [default: ./data/mesh.dat]\n");
        printf("  -meshout <val>   : output mesh file [default: ./data/mesh_?.dat]\n");
        printf("  -assemble <val>  : assemble ofemparion [default: ab]\n");	
        printf("                     ab  |  assemble the mat & rhs;\n");
        printf("                      a  |  assemble the mat;\n");
        printf("                      b  |  assemble the rhs;\n");
        printf("  -refine <val>    : refine level [default: 8]\n");
        printf("  -output <val>    : mesh output ofemparion [default: 0]\n");
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
	mesh_init (&mesh, "./data/mesh.dat");

	// Step 1.5: read or build auxiliary mesh info
    // If there is already mesh_aux data available, you can use the following fct to init it:    
    //	  mesh_aux_init (&mesh, &mesh_aux, "mesh.aux");
    // otherwise, just build the auxiliary information for the mesh
	mesh_aux_build(&mesh, &mesh_aux);
    
	// Step 2: refine mesh
	clock_t mesh_refine_s = clock();
	for (i=0; i<fempar.refine_lvl; ++i) mesh_refine(&mesh, &mesh_aux);
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
	
	// Step 4: solve the linear system with AMG
    {
        AMG_param amgparam; // parameters for AMG
        dvector x; // just on DOF, not including boundary
        
        // Print problem size
        printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
        printf("b: n = %d\n", b.row);
	
        fasp_param_amg_init(&amgparam);    // set AMG param with default values
        amgparam.print_level = PRINT_SOME; // print some AMG message
        amgparam.maxit       = 25;         // max number of iterations
        
        fasp_dvec_alloc(A.row, &x); fasp_dvec_set(A.row, &x, 0.0); // initial guess
        
        fasp_solver_amg(&A, &b, &x, &amgparam); // solve Ax = b with AMG

        for ( i=0; i<dof.row; ++i) uh.val[dof.val[i]] = x.val[i]; // save x to uh

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
    
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
