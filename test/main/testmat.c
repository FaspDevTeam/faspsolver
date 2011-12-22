/**
 *		Test matrix properties 
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/31/2009.
 *      Modified by Chensong Zhang on 09/02/2011.
 *
 *------------------------------------------------------
 *
 */

/*! \file testmat.c
 *  \brief Test matrix properties
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * This is the main function for test purpose. It contains five steps:
 */
int main(int argc, const char * argv[]) 
{		
	char *inputfile="ini/input.dat";
	input_param Input;
	fasp_param_input(inputfile,&Input);
	
	char filename1[512], *datafile1;	
	strncpy(filename1,Input.workdir,128);
	
	/** Step 1. Assemble matrix and right-hand side */ 
	dCSRmat A;
	datafile1="matP1.dat";
	strcat(filename1,datafile1);
	fasp_dcoo_read(filename1, &A);
	
	/** Step 2. Check matrix properties */
	
	/* check sparse pattern */
	char *bmpfile="out/matrix.bmp";	/* Output the matrix as BMP file */
	fasp_dcsr_plot(&A, bmpfile, 200);
	
	/* get diagonal part */
	fasp_check_diagpos(&A);
	
	/* check symmetry */
	fasp_check_symm(&A);
	
	/* check diagnoal dominance */
	fasp_check_diagdom(&A);
	
	/* check max and min eigenvalues */
	
	
	/** Step 3. Select a solver */
	
	/** Step 4. Solve the system */ 
	
	/** Step 5.	Post processing */	
	char *matfile="out/matrix.out";	/* Output the matrix in coordinate format */	
	fasp_dcsr_write(matfile, &A);
	fasp_dcsr_free(&A);
	
	return 0;
}
