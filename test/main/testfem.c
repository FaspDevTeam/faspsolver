/**
 *		Test FASP FEM assembling routines. 
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 09/12/2011.
 *
 *------------------------------------------------------
 *
 */

/*! \file testfem.c
 *  \brief The main test function for FASP FEM assembling.
 */

#include "fasp.h"
#include "fasp_functs.h"

// Declare P1 assembling function
extern int assemble(dCSRmat *ptr_A, dvector *ptr_b, int levelNum);

/**
 * \fn int main (int argc, const char * argv[])
 *
 * This is the main function for test.
 */
int main (int argc, const char * argv[]) 
{
	dCSRmat A;
	dvector b;

	printf("Test P1 FEM assmbling routine ...\n");
    
	// Assemble A and b -- P1 FE discretization for Poisson.
    assemble(&A,&b,9);
    
    // Print problem size
    printf("A: m = %d, n = %d, nnz = %d\n", A.row, A.col, A.nnz);
    printf("b: n = %d\n", b.row);
    
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
    
	return SUCCESS;
}
