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
 * \brief This is the main function for testing matrix propertities.
 *
 * \author Chensong Zhang
 * \date   03/31/2009
 * 
 * Modified by Chensong Zhang on 03/19/2012 
 */
int main(int argc, const char * argv[]) 
{		
	char *inputfile="ini/input.dat";
	input_param Input;
	fasp_param_input(inputfile,&Input);
	
	char filename1[512], *datafile1;	
	strncpy(filename1,Input.workdir,128);
	
	// Read matrix for testing
	dCSRmat A;
	datafile1="matP1.dat";
	strcat(filename1,datafile1);
	fasp_dcoo_read(filename1, &A);
		
	// Check sparse pattern
	char *bmpfile="out/matrix.bmp";	/* Output the matrix as BMP file */
	fasp_dcsr_plot(&A, bmpfile, 200);
	
	// Check symmetry
	fasp_check_symm(&A);
	
	// Check diagonal positivity
	fasp_check_diagpos(&A);
	
	// Check diagnoal dominance
	fasp_check_diagdom(&A);
	
	// Output matrix in COO format	
	char *matfile="out/matrix.out";	/* Output the matrix in coordinate format */	
	fasp_dcsr_write(matfile, &A);
	
    // Clean up memory
    fasp_dcsr_free(&A);
	
	return SUCCESS;
}
