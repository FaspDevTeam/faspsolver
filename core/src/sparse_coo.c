/*! \file sparse_coo.c
 *  \brief Functions for COO sparse matrices. 
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dCOOmat fasp_dcoo_create (INT m, INT n, INT nnz)
 *
 * \brief Create IJ sparse matrix data memory space
 *
 * \param m    number of rows
 * \param n    number of columns
 * \param nnz  number of nonzeros
 * \return A   the new dCOOmat matrix
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
dCOOmat fasp_dcoo_create (INT m, 
                          INT n, 
                          INT nnz)
{			
	dCOOmat A;
	
	A.I   = (INT *)fasp_mem_calloc(nnz, sizeof(INT)); 	
	A.J   = (INT *)fasp_mem_calloc(nnz, sizeof(INT)); 	
	A.val = (REAL *)fasp_mem_calloc(nnz, sizeof(REAL)); 
	
#if CHMEM_MODE	
	total_alloc_mem += nnz*(sizeof(REAL)+sizeof(INT)*2);
#endif
	
	A.row=m; A.col=n; A.nnz=nnz;
	
	return A;
}

/**
 * \fn void fasp_dcoo_free (dCOOmat *A)
 *
 * \brief Free IJ sparse matrix data memeory space
 *
 * \param A   pointer to the dCOOmat matrix
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_dcoo_free (dCOOmat *A)
{			
	if (A==NULL) return;
	
	fasp_mem_free(A->I);   A->I   = NULL;
	fasp_mem_free(A->J);   A->J   = NULL;
	fasp_mem_free(A->val); A->val = NULL;
}

/**
 * \fn void fasp_dcoo_shift (dCOOmat *A, INT offset)
 *
 * \brief Reindex a REAL matrix in IJ format to make the index starting from 0 or 1.
 *
 * \param A       pointer to IJ matrix
 * \param  offset  size of offset (1 or -1)
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
void fasp_dcoo_shift (dCOOmat *A, INT offset)
{
	const INT nnz=A->nnz;
	INT i, *ai=A->I, *aj=A->J;
	
	if (offset == 0) offset = ISTART;
	
	for (i=0;i<nnz;++i) {	
		ai[i]+=offset; aj[i]+=offset;
	}
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
