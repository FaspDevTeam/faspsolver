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
 *
 * \brief Initialize sparse matrix on structured grid
 *
 * \param A pointer to the dSTRmat matrix
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
 * \fn dSTRmat fasp_dstr_create (INT nx, INT ny, INT nz, INT nc, INT nband, INT *offsets)
 *
 * \brief Create STR sparse matrix data memory space
 *
 * \param nx      integer, number of grids in x direction
 * \param ny      integer, number of grids in y direction
 * \param nz      integer, number of grids in z direction
 * \param nc      integer, number of components 
 * \param nband   integer, number of off-diagonal bands 
 * \param offsets integer, shift from diagonal
 *
 * \return        the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010 
 */
dSTRmat fasp_dstr_create (INT nx, 
                          INT ny, 
                          INT nz, 
                          INT nc, 
                          INT nband, 
                          INT *offsets)
{	
	dSTRmat A;
	
	unsigned INT i;
	
	A.nx=nx; A.ny=ny; A.nz=nz;
	A.nc=nc;
	A.nxy=A.nx*A.ny;
	A.ngrid=A.nxy*A.nz;
	A.nband=nband;
	
	A.offsets=(int*)fasp_mem_calloc(nband, sizeof(INT));
	
	for (i=0;i<nband;++i) A.offsets[i]=offsets[i];
	
	A.diag=(REAL*)fasp_mem_calloc(A.ngrid*A.nc*A.nc, sizeof(REAL));
	
	A.offdiag=(REAL**)fasp_mem_calloc(nband, sizeof(REAL*));
	
	for(i=0;i<A.nband;++i) {
		A.offdiag[i]=(REAL*)fasp_mem_calloc((A.ngrid-ABS(A.offsets[i]))*A.nc*A.nc, sizeof(REAL));
	}
	
	return(A);
}


/**
 * \fn void fasp_dstr_alloc (INT nx, INT ny, INT nz, INT nxy, INT ngrid, INT nband, 
 *                           INT nc, INT *offsets, dSTRmat *A)
 *
 * \brief Allocate STR sparse matrix memory space
 *
 * \param nx     integer, number of grids in x direction
 * \param ny     integer, number of grids in y direction
 * \param nz     integer, number of grids in z direction
 * \param nxy    integer, number of grids in x-y plane
 * \param ngrid  integer, number of grids
 * \param nband  integer, number of off-diagonal bands 
 * \param nc     integer, number of components 
 * \param offsets integer, shift from diagonal
 * \param A       pointer to the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010  
 */
void fasp_dstr_alloc(INT nx, 
                     INT ny, 
                     INT nz, 
                     INT nxy, 
                     INT ngrid, 
                     INT nband, 
                     INT nc,
                     INT *offsets, 
                     dSTRmat *A)
{	
	INT i;
	
	A->nx=nx;
	A->ny=ny;
	A->nz=nz;
	A->nxy=nxy;
	A->ngrid=ngrid;
	A->nband=nband;
	A->nc=nc;
	
	A->offsets=(int*)fasp_mem_calloc(nband, sizeof(INT));
	
	for (i=0;i<nband;++i) A->offsets[i]=offsets[i];
	
	A->diag=(REAL*)fasp_mem_calloc(ngrid*nc*nc, sizeof(REAL));
	
	A->offdiag = (REAL **)fasp_mem_calloc(A->nband, sizeof(REAL*));
	
	for (i=0;i<nband;++i) {
		A->offdiag[i]=(REAL*)fasp_mem_calloc((ngrid-ABS(offsets[i]))*nc*nc, sizeof(REAL));
	}
}

/**
 * \fn void fasp_dstr_free (dSTRmat *A)
 *
 * \brief Free STR sparse matrix data memeory space
 *
 * \param A   pointer to the dSTRmat matrix
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 05/17/2010 
 */
void fasp_dstr_free (dSTRmat *A)
{		
	unsigned INT i;
	
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
 *
 * \brief copy a dSTRmat to a new one A1=A
 *
 * \param A pointer to the dSTRmat matrix
 * \param A1 pointer to the dSTRmat matrix
 * 
 * \author Zhiyang Zhou
 * \date 04/21/2010  
 */
void fasp_dstr_cp (dSTRmat *A, 
                   dSTRmat *A1)
{		
    const INT nc2 = (A->nc)*(A->nc);

    INT i;
	A1->nx=A->nx;
	A1->ny=A->ny;
	A1->nz=A->nz;
	A1->nxy=A->nxy;
	A1->ngrid=A->ngrid;
	A1->nc=A->nc;
	A1->nband=A->nband;
		
	memcpy(A1->offsets,A->offsets,(A->nband)*sizeof(INT));
	memcpy(A1->diag,A->diag,(A->ngrid*nc2)*sizeof(REAL));
	for (i=0;i<A->nband;++i)
		memcpy(A1->offdiag[i],A->offdiag[i],((A->ngrid - ABS(A->offsets[i]))*nc2)*sizeof(REAL)); 
	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
