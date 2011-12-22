/*! \file sparse_str.c
 *  \brief Simple operations for STR sparse matrices.
 *
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_dstr_init (dSTRmat *A)
 * \brief Initialize sparse matrix on structured grid
 *
 * \param *A pointer to the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010 
 */
void fasp_dstr_init (dSTRmat *A)
{		
	A->nx=0;
	A->ny=0;
	A->nz=0;
	A->nxy=0;
	A->ngrid=0;
	A->nband=0;
	A->nc=0;	
	A->offsets=NULL;
	A->diag=NULL;
	A->offdiag=NULL;
}

/**
 * \fn dSTRmat fasp_dstr_create (int nx, int ny, int nz, int nc, int nband, int *offsets)
 * \brief Create STR sparse matrix data memory space
 *
 * \param nx      integer, number of grids in x direction
 * \param ny      integer, number of grids in y direction
 * \param nz      integer, number of grids in z direction
 * \param nc      integer, number of components 
 * \param nband   integer, number of off-diagonal bands 
 * \return A      the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010 
 */
dSTRmat fasp_dstr_create (int nx, 
													int ny, 
													int nz, 
													int nc, 
													int nband, 
													int *offsets)
{	
	dSTRmat A;
	
	unsigned int i;
	
	A.nx=nx; A.ny=ny; A.nz=nz;
	A.nc=nc;
	A.nxy=A.nx*A.ny;
	A.ngrid=A.nxy*A.nz;
	A.nband=nband;
	
	A.offsets=(int*)fasp_mem_calloc(nband, sizeof(int));
	
	for (i=0;i<nband;++i) A.offsets[i]=offsets[i];
	
	A.diag=(double*)fasp_mem_calloc(A.ngrid*A.nc*A.nc, sizeof(double));
	
	A.offdiag=(double**)fasp_mem_calloc(nband, sizeof(double*));
	
	for(i=0;i<A.nband;++i) {
		A.offdiag[i]=(double*)fasp_mem_calloc((A.ngrid-ABS(A.offsets[i]))*A.nc*A.nc, sizeof(double));
	}
	
	return(A);
}


/**
 * \fn void fasp_dstr_alloc(int nx, int ny, int nz, int nxy, int ngrid, int nband, int nc, int *offsets, dSTRmat *A)
 * \brief Allocate STR sparse matrix memory space
 *
 * \param nx     integer, number of grids in x direction
 * \param ny     integer, number of grids in y direction
 * \param nz     integer, number of grids in z direction
 * \param nxy    integer, number of grids in x-y plane
 * \param ngrid  integer, number of grids
 * \param nband  integer, number of off-diagonal bands 
 * \param nc     integer, number of components 
 * \param *A     pointer to the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010  
 */
void fasp_dstr_alloc(int nx, 
										 int ny, 
										 int nz, 
										 int nxy, 
										 int ngrid, 
										 int nband, 
										 int nc,
										 int *offsets, 
										 dSTRmat *A)
{	
	int i;
	
	A->nx=nx;
	A->ny=ny;
	A->nz=nz;
	A->nxy=nxy;
	A->ngrid=ngrid;
	A->nband=nband;
	A->nc=nc;
	
	A->offsets=(int*)fasp_mem_calloc(nband, sizeof(int));
	
	for (i=0;i<nband;++i) A->offsets[i]=offsets[i];
	
	A->diag=(double*)fasp_mem_calloc(ngrid*nc*nc, sizeof(double));
	
	A->offdiag = (double **)fasp_mem_calloc(A->nband, sizeof(double*));
	
	for (i=0;i<nband;++i) {
		A->offdiag[i]=(double*)fasp_mem_calloc((ngrid-ABS(offsets[i]))*nc*nc, sizeof(double));
	}
}

/**
 * \fn void fasp_dstr_free (dSTRmat *A)
 * \brief Free STR sparse matrix data memeory space
 *
 * \param *A   pointer to the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010 
 */
void fasp_dstr_free (dSTRmat *A)
{		
	unsigned int i;
	
	fasp_mem_free(A->offsets);
	fasp_mem_free(A->diag);
	for (i=0;i<A->nband;++i) fasp_mem_free(A->offdiag[i]);
	
	A->nx=0;
	A->ny=0;
	A->nz=0;
	A->nxy=0;
	A->ngrid=0;
	A->nband=0;
	A->nc=0;	
}

/**
 * \fn void fasp_dstr_cp (dSTRmat *A, dSTRmat *A1)
 * \brief copy a dSTRmat to a new one A1=A
 * \param *A pointer to the dSTRmat matrix
 * \param *A1 pointer to the dSTRmat matrix
 * 
 * \author Zhiyang Zhou
 * \date 04/21/2010  
 */
void fasp_dstr_cp (dSTRmat *A, 
									 dSTRmat *A1)
{		
	int i;
	A1->nx=A->nx;
	A1->ny=A->ny;
	A1->nz=A->nz;
	A1->nxy=A->nxy;
	A1->ngrid=A->ngrid;
	A1->nc=A->nc;
	A1->nband=A->nband;
	
	int nc2 = (A->nc)*(A->nc);
	
	memcpy(A1->offsets,A->offsets,(A->nband)*sizeof(int));
	memcpy(A1->diag,A->diag,(A->ngrid*nc2)*sizeof(double));
	for (i=0;i<A->nband;++i)
		memcpy(A1->offdiag[i],A->offdiag[i],((A->ngrid - ABS(A->offsets[i]))*nc2)*sizeof(double)); 
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
