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
	param_test pt;// parameter for testfem
	param_init(&pt);
	print_usage = param_set(argc, argv, &pt);
    
	if (print_usage)
	{
		printf("\nUsage: %s [<options>]\n", argv[0]);
		printf("  -meshin <val>    : input mesh file [default: ./data/testmesh.dat]\n");
		printf("  -meshout <val>   : output mesh file [default: ./data/mesh_?.dat]\n");
		printf("  -assemble <val>  : assemble option [default: ab]\n");	
        printf("                     ab  |  assemble the mat & rhs;\n");
		printf("                      a  |  assemble the mat;\n");
		printf("                      b  |  assemble the rhs;\n");
		printf("  -refine <val>    : refine level [default: 0]\n");
		printf("  -output <val>    : mesh output option [default: 0]\n");
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
    
	// Step 1: reading mesh info
	mesh_init (&mesh, "./data/testmesh.dat");

    // If there is already mesh_aux data available, you can use the following fct to init it:    
    //	  mesh_aux_init (&mesh, &mesh_aux, "dd.dat");
    // otherwise, just build the auxiliary information for the mesh
	mesh_aux_build(&mesh, &mesh_aux);
    
	// Step 2: refinement
	clock_t mesh_refine_s = clock();
	for (i=0; i<pt.refine_lvl; ++i) mesh_refine(&mesh, &mesh_aux);
	clock_t mesh_refine_e = clock();
    double mesh_refine_time = (double)(mesh_refine_e - mesh_refine_s)/(double)(CLOCKS_PER_SEC);
    printf("Mesh refinement costs... %8.4f seconds\n", mesh_refine_time);
    
	// Step 3: assembling the linear system and right-hand side
    clock_t assemble_s = clock();
    setup_poisson(&A, &b, &mesh, &mesh_aux, &pt, &uh, &dof);
    clock_t assemble_e = clock();
    double assemble_time = (double)(assemble_e - assemble_s)/(double)(CLOCKS_PER_SEC);
    printf("Assembling costs........ %8.4f seconds\n", assemble_time);
    
	// Clean auxiliary mesh info
	mesh_aux_free(&mesh_aux);
	
    // Print problem size
    printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
    printf("b: n = %d\n", b.row);
	
	// Step 4: solve the linear system with AMG
    dvector x; // just DOF, not including boundary
    {
        AMG_param amgparam; // parameters for AMG
        
        fasp_param_amg_init(&amgparam); // set AMG param with default values
        amgparam.print_level = PRINT_SOME; // print some AMG message
        amgparam.maxit = 100; // max iter number = 100 
        
        fasp_dvec_alloc(A.row, &x); 
        fasp_dvec_set(A.row, &x, 0.0);
        
        fasp_solver_amg(&A, &b, &x, &amgparam);
    }
	
	// Step 5: estimate L2 error
	for ( i=0; i<dof.row; ++i) uh.val[dof.val[i]] = x.val[i];
    
	double l2error = get_l2_error_poisson(&(mesh.node), &(mesh.elem), &uh, pt.num_qp_rhs);
	
    printf("\n==============================================================\n");
	printf("L2 error of FEM is %g\n", l2error);
    printf("==============================================================\n\n");
    
	// clean 
	mesh_free(&mesh);
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
    fasp_ivec_free(&dof);
    fasp_dvec_free(&x);
    
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
