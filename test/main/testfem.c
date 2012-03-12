/**
 *    Testing the FASP FEM assembling routines. 
 *
 *------------------------------------------------------
 *
 */

/*! \file testfem.c
 *  \brief The main test function for FASP FEM assembling.
 */

#include "fasp.h"
#include "fasp_functs.h"
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
 * \note   Modified by Chensong Zhang and Feiteng Huang on 03/07/2012
 */
int main (int argc, const char * argv[]) 
{
    // Set default values
    const char *assemble_option = "a&b";
    int status      = SUCCESS;
	int arg_index   = 1;
	int print_usage = 0;
	int oo          = 0;
	int refine_lvl  = 8;
	int input_flag  = 1;
	int num_qp_rhs  = 3; // enough for P1 smooth right-hand-side
	int num_qp_mat  = 1;
    
    // Set default input/output mesh files
	const char *meshIn  = "./data/testmesh.dat";
	const char *meshOut = "./data/mesh_";
    
    // Local variables
	char filename[128];
	char *matFile = "./out/mat";
	char *rhsFile = "./out/rhs";
	dCSRmat A;
	dvector b;
    
    // Read input from arguments
	while (arg_index < argc)
	{
		if (argc%2 == 0)
		{
			print_usage = 1;
			break;
		}
        
		if ( strcmp(argv[arg_index], "-help") == 0 )
		{
			print_usage = 1;
			break;
		}
        
		if ( strcmp(argv[arg_index], "-meshin") == 0 )
		{
			arg_index ++;
			meshIn = argv[arg_index++];
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if ( strcmp(argv[arg_index], "-meshout") == 0 )
		{
			arg_index ++;
			meshOut = argv[arg_index++];
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if ( strcmp(argv[arg_index], "-assemble") == 0 )
		{
			arg_index ++;
			assemble_option = argv[arg_index++];
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if ( strcmp(argv[arg_index], "-refine") == 0 )
		{
			arg_index ++;
			refine_lvl = atoi(argv[arg_index++]);
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if ( strcmp(argv[arg_index], "-output") == 0 )
		{
			arg_index ++;
			oo = atoi(argv[arg_index++]);
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if ( strcmp(argv[arg_index], "-quad_rhs") == 0 )
		{
			arg_index ++;
			num_qp_rhs = atoi(argv[arg_index++]);
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if ( strcmp(argv[arg_index], "-quad_mat") == 0 )
		{
			arg_index ++;
			num_qp_mat = atoi(argv[arg_index++]);
			input_flag = 0;
		}
		if (arg_index >= argc) break;
        
		if (input_flag)
		{
			print_usage = 1;
			break;
		}
	}
    
	if (print_usage)
	{
		printf("\nUsage: %s [<options>]\n", argv[0]);
		printf("  -meshin <val>    : input mesh file [default: ./data/testmesh.dat]\n");
		printf("  -meshout <val>   : output mesh file [default: ./data/mesh_?.dat]\n");
		printf("  -assemble <val>  : assemble option [default: a&b]\n");	
        printf("                     a&b |  assemble the mat & rhs;\n");
		printf("                      a  |  assemble the mat;\n");
		printf("                      b  |  assemble the rhs;\n");
		printf("  -refine <val>    : refine level [default: 0]\n");
		printf("  -output <val>    : mesh output option [default: 0]\n");
		printf("  -quad_rhs <val>  : quad points for rhs [default: 16]\n");
		printf("  -quad_mat <val>  : quad points for mat [default: 1]\n");
		printf("  -help            : print this help message\n\n");
		exit(status);
	}
    
	// Assemble A and b -- P1 FE discretization for Poisson.
    setup_poisson(&A, &b, refine_lvl+1, meshIn, meshOut, oo, assemble_option, 
                  num_qp_rhs, num_qp_mat);
    
    // Write A and b
	sprintf(filename, "%s_coo_%d.dat", matFile, refine_lvl+1);
	fasp_dcsr_write (filename, &A);
    
	sprintf(filename, "%s_%d.dat", rhsFile, refine_lvl+1);
	fasp_dvec_write (filename, &b);
    
    // Print problem size
    printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
    printf("b: n = %d\n", b.row);
    
    // Solve A x = b with AMG
    {
        AMG_param       amgparam; // parameters for AMG
        fasp_param_amg_init(&amgparam); // set AMG param with default values
        amgparam.print_level = PRINT_SOME; // print some AMG message
        amgparam.maxit = 20; // max iteration number = 20 
        
        dvector x;
        fasp_dvec_alloc(A.row, &x); 
        fasp_dvec_set(A.row,&x,0.0);
        
        fasp_solver_amg(&A, &b, &x, &amgparam); 
        fasp_dvec_free(&x);
    }
    
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
    
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
