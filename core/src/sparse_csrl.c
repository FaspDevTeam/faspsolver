/*! \file sparse_csrl.c
 *  \brief Functions for sparse matrices in CSRL format
 *
 *  \note For details of CSRL format, refer to 
 *        Optimizaing sparse matrix vector product computations using unroll and jam
 *        by John Mellor-Crummey and John Garvin, Tech Report Rice Univ, Aug 2002.
 *
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dCSRLmat * fasp_dcsrl_create (INT num_rows, INT num_cols, INT num_nonzeros) 
 *
 * \brief Create a dCSRLmat object
 *
 * \param num_rows      Number of rows
 * \param num_cols      Number of cols
 * \param num_nonzeros  Number of nonzero entries
 *
 * \author Zhou Zhiyang
 * \date   2011/01/07
 */
dCSRLmat * fasp_dcsrl_create (INT num_rows, 
                              INT num_cols, 
                              INT num_nonzeros)
{
	dCSRLmat *A = NULL;
	A = (dCSRLmat *)fasp_mem_calloc(1, sizeof(dCSRLmat));
	A -> num_rows     = num_rows;
	A -> num_cols     = num_cols;
	A -> num_nonzeros = num_nonzeros;
	A -> nzdifnum     = NULL;
	A -> rowindex     = NULL;
	A -> rowstart     = NULL;
	A -> ja           = NULL;
	A -> data         = NULL;
	
	return A;
}

/**
 * \fn void fasp_dcsrl_free ( dCSRLmat *A )
 *
 * \brief Destroy a dCSRLmat object
 *
 * \param A   Pointer to the dCSRLmat type matrix
 *
 * \author Zhou Zhiyang
 * \date   2011/01/07
 */
void fasp_dcsrl_free (dCSRLmat *A)
{
	if (A)
	{  
		if (A -> nzdifnum) free(A -> nzdifnum);
		if (A -> rowindex) free(A -> rowindex);
		if (A -> rowstart) free(A -> rowstart);
		if (A -> ja)       free(A -> ja);
		if (A -> data)     free(A -> data);
		free(A);
	}
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
