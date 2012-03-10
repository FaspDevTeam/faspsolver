/*! \file blas_str.c
 *  \brief BLAS operations for STR sparse matrices.
 *
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

static inline void smat_amxv_nc3(REAL alpha, REAL *a, REAL *b, REAL *c);
static inline void smat_amxv_nc5(REAL alpha, REAL *a, REAL *b, REAL *c);
static inline void smat_amxv(REAL alpha, REAL *a, REAL *b,INT n, REAL *c);
static inline void blkcontr_str(INT start_data, INT start_vecx, INT start_vecy, INT nc, 
                                REAL *data, REAL *x, REAL *y);
static inline void spaaxpy_str_2D_scalar(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_2D_nc3(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_2D_nc5(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_2D_block(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_3D_scalar(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_3D_nc3(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_3D_nc5(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_3D_block(REAL alpha, dSTRmat *A, REAL *x, REAL *y);
static inline void spaaxpy_str_general(REAL alpha, dSTRmat *A, REAL *x, REAL *y);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_dstr_aAxpy (REAL alpha, dSTRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \author Zhiyang Zhou, Xiaozhe Hu, Shiquan Zhang
 * \date 2010/10/15
 */
void fasp_blas_dstr_aAxpy (REAL alpha, 
                           dSTRmat *A, 
                           REAL *x, 
                           REAL *y)
{
	
	switch (A->nband)
	{
		case 4:
			
			switch (A->nc)
		{
			case 1:
				spaaxpy_str_2D_scalar(alpha, A, x, y);
				break;
				
			case 3:
				spaaxpy_str_2D_nc3(alpha, A, x, y);
				break;
				
			case 5:
				spaaxpy_str_2D_nc5(alpha, A, x, y);
				break;
				
			default:
				spaaxpy_str_2D_block(alpha, A, x, y);
				break;
		}
			
			break;
			
		case 6:
			
			switch (A->nc)
		{
			case 1:
				spaaxpy_str_3D_scalar(alpha, A, x, y);
				break;
				
			case 3:
				spaaxpy_str_3D_nc3(alpha, A, x, y);
				break;
				
			case 5:
				spaaxpy_str_3D_nc5(alpha, A, x, y);
				break;
				
			default:
				spaaxpy_str_3D_block(alpha, A, x, y);
				break;
		}	
			break;
			
		default:
			spaaxpy_str_general(alpha, A, x, y);
			break;
	}
	
}

/*!
 * \fn INT fasp_dstr_diagscale (dSTRmat *A, dSTRmat *B)
 *
 * \brief B=D^{-1}A
 *
 * \param A pointer to a 'dSTRmat' type matrix A
 * \param B pointer to a 'dSTRmat' type matrix B
 * 
 * \author Shiquan Zhang
 * \date 2010/10/15
 */
INT fasp_dstr_diagscale (dSTRmat *A, 
                         dSTRmat *B)
{
	INT ngrid=A->ngrid, nc=A->nc,nband=A->nband;
	INT nc2=nc*nc, size=ngrid*nc2;
	INT i,j,ic2,nb,nb1;
	REAL *diag=(REAL *)fasp_mem_calloc(size,sizeof(REAL));
	
	fasp_array_cp(size,A->diag,diag);
	
	fasp_dstr_alloc(A->nx, A->ny, A->nz,A->nxy,ngrid, nband,nc,A->offsets, B);
	
	//compute diagnal elements of B
	for (i=0;i<ngrid;++i)
	{  ic2=i*nc2;
        for (j=0;j<nc2;++j)
		{ if(j/nc == j%nc) B->diag[ic2+j]=1; else B->diag[ic2+j]=0;}
	}
    
	for (i=0;i<ngrid;++i) fasp_blas_smat_inv(&(diag[i*nc2]),nc);
	
	for (i=0;i<nband;++i) 
	{
		nb=A->offsets[i];
		nb1=abs(nb);
		if (nb<0)
			for (j=0;j<ngrid-nb1;++j) fasp_blas_smat_mul(&(diag[(j+nb1)*nc2]),&(A->offdiag[i][j*nc2]),&(B->offdiag[i][j*nc2]),nc);
		else
			for (j=0;j<ngrid-nb1;++j) fasp_blas_smat_mul(&(diag[j*nc2]),&(A->offdiag[i][j*nc2]),&(B->offdiag[i][j*nc2]),nc);
		
    }
	
	fasp_mem_free(diag);
	
	return (0);	
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn void smat_amxv_nc3 (REAL alpha, REAL *a, REAL *b, REAL *c)
 *
 * \brief Matrix-vector multiplication c = alpha*a*b + c where a is a 3*3 full matrix
 *
 * \param alpha real number
 * \param a pointer to the REAL vector which stands a 3*3 matrix
 * \param b pointer to the REAL vector with length 3
 * \param c pointer to the REAL vector with length 3
 */
static inline void smat_amxv_nc3 (REAL alpha, 
                                  REAL *a, 
                                  REAL *b, 
                                  REAL *c)
{ 
	c[0] += alpha*(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
	c[1] += alpha*(a[3]*b[0] + a[4]*b[1] + a[5]*b[2]);
	c[2] += alpha*(a[6]*b[0] + a[7]*b[1] + a[8]*b[2]);
}

/**
 * \fn void smat_amxv_nc5(REAL alpha, REAL *a, REAL *b, REAL *c)
 *
 * \brief  Matrix-vector multiplication c = alpha*a*b + c where a is a 5*5 full matrix
 *
 * \param alpha real number
 * \param a pointer to the REAL vector which stands a 5*5 matrix
 * \param b pointer to the REAL vector with length 5
 * \param c pointer to the REAL vector with length 5
 */
static inline void smat_amxv_nc5 (REAL alpha, 
                                  REAL *a, 
                                  REAL *b, 
                                  REAL *c)
{ 
	c[0] += alpha*(a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3] * b[3] + a[4] * b[4]);
	c[1] += alpha*(a[5]*b[0] + a[6]*b[1] + a[7]*b[2] + a[8] * b[3] + a[9] * b[4]);
	c[2] += alpha*(a[10]*b[0] + a[11]*b[1] + a[12]*b[2] + a[13] * b[3] + a[14] * b[4]);
	c[3] += alpha*(a[15]*b[0] + a[16]*b[1] + a[17]*b[2] + a[18] * b[3] + a[19] * b[4]);
	c[4] += alpha*(a[20]*b[0] + a[21]*b[1] + a[22]*b[2] + a[23] * b[3] + a[24] * b[4]);
	
}

/**
 * \fn void smat_amxv (REAL alpha, REAL *a, REAL *b, INT n, REAL *c)
 *
 * \brief  Matrix-vector multiplication c = alpha*a*b + c where a is a n*n full matrix
 *
 * \param alpha real number
 * \param a pointer to the REAL vector which stands a n*n matrix
 * \param b pointer to the REAL vector with length n
 * \param c pointer to the REAL vector with length n
 * \param n the dimension of the matrix
 */
static inline void smat_amxv (REAL alpha, 
                              REAL *a, 
                              REAL *b,
                              INT n, 
                              REAL *c)
{ 
	INT i,j;
	INT in;	
	
	for (i=0;i<n;++i)
	{
		in = i*n;
		for (j=0;j<n;++j)
			c[i] += alpha*a[in+j]*b[j];
	}  //! end for (i=0;i<n;++i)
	
	return;
}

/**
 * \fn static void blkcontr_str (INT start_data, INT start_vecx, INT start_vecy, INT nc, 
 *                               REAL *data, REAL *x, REAL *y)
 *
 * \brief contribute the block computation 'P*z' to 'y', where 'P' is a nc*nc  
 *        full matrix stored in 'data' from the address 'start_data', and 'z' 
 *        is a nc*1 vector stored in 'x' from the address 'start_vect'. 
 *
 * \param start_data start_data starting position
 * \param start_vecx start_data starting position
 * \param start_vecy start_data starting position
 * \param nc the dimension of the submatrix
 * \param data pointer to matrix data
 * \param x pointer to real array
 * \param y pointer to real array
 * \date 04/24/2010
 */
static inline void blkcontr_str (INT start_data, 
                                 INT start_vecx, 
                                 INT start_vecy, 
                                 INT nc, 
                                 REAL *data, 
                                 REAL *x, 
                                 REAL *y)
{
	INT i,j,k,m;
	for (i = 0; i < nc; i ++) 
	{
		k = start_data + i*nc;
		m = start_vecy + i;
		for (j = 0; j < nc; j ++) 
		{
			y[m] += data[k+j]*x[start_vecx+j];
		}
	}
} 

/**
 * \fn void spaaxpy_str_2D_scalar (REAL alpha, dSTRmat *A, REAL *x, REAL *y)  
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 5 banded 2D structured matrix (nc = 1)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsets of the five bands have to be (-1, +1, -nx, +nx) for nx != 1 and (-1,+1,-ny,+ny) 
 *       for nx = 1, but the order can be arbitrary. 
 */
static inline void spaaxpy_str_2D_scalar (REAL alpha, 
                                          dSTRmat *A, 
                                          REAL *x, 
                                          REAL *y)
{
	unsigned INT i;
	INT idx1,idx2;
	INT end1, end2;
	INT nline;
	
	//! information of A
	INT nx = A->nx;
	INT ngrid = A->ngrid;  //! number of grids
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, *offdiag3=NULL;
	
	if (nx == 1)
	{
		nline = A->ny;
	}
	else
	{
		nline = nx;
	}
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nline)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nline)
		{
			offdiag3 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 2D scalar case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nline;
	
	y[0] += alpha*(diag[0]*x[0] + offdiag1[0]*x[1] + offdiag3[0]*x[nline]);
	
	for (i=1; i<nline; ++i){
		idx1 = i-1;
		y[i] += alpha*(offdiag0[idx1]*x[idx1] + diag[i]*x[i] + offdiag1[i]*x[i+1] + 
                       offdiag3[i]*x[i+nline]);
	}
	
	for (i=nline; i<end2; ++i){
		idx1 = i-1; 
		idx2 = i-nline; 
		y[i] += alpha*(offdiag2[idx2]*x[idx2] + offdiag0[idx1]*x[idx1] + 
                       diag[i]*x[i] + offdiag1[i]*x[i+1] + offdiag3[i]*x[i+nline]);
	}
	
	for (i=end2; i<end1; ++i){
		idx1 = i-1; 
		idx2 = i-nline; 
		y[i] += alpha*(offdiag2[idx2]*x[idx2] + offdiag0[idx1]*x[idx1] + 
                       diag[i]*x[i] + offdiag1[i]*x[i+1]);
	}
	
	idx1 = end1-1; 
	idx2 = end1-nline; 
	y[end1] += alpha*(offdiag2[idx2]*x[idx2] + 
                      offdiag0[idx1]*x[idx1] + diag[end1]*x[end1]);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_2D_nc3(REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 5 banded 2D structured matrix (nc = 3)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsets of the five bands have to be (-1, +1, -nx, +nx) for nx != 1 and (-1,+1,-ny,+ny) 
 *       for nx = 1, but the order can be arbitrary. 
 */
static inline void spaaxpy_str_2D_nc3 (REAL alpha, 
                                       dSTRmat *A, 
                                       REAL *x, 
                                       REAL *y)
{
	INT i;
	INT idx,idx1,idx2;
	INT matidx, matidx1, matidx2;
	INT end1, end2;
	INT nline, nlinenc;
	
	//! information of A
	INT nx = A->nx;
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, *offdiag3=NULL;
	
	if (nx == 1)
	{
		nline = A->ny;
	}
	else
	{
		nline = nx;
	}
	nlinenc = nline*nc;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nline)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nline)
		{
			offdiag3 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 2D of nc=3 case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nline;
	
	smat_amxv_nc3(alpha, diag, x, y);
	smat_amxv_nc3(alpha, offdiag1, x+nc, y);
	smat_amxv_nc3(alpha, offdiag3, x+nlinenc, y);
	
	for (i=1; i<nline; ++i){
		idx = i*nc;
		matidx = idx*nc;
		idx1 = idx - nc;
		matidx1 = idx1*nc;
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc3(alpha, offdiag3+matidx, x+idx+nlinenc, y+idx);
	}
	
	
	for (i=nx; i<end2; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nlinenc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc3(alpha, offdiag3+matidx, x+idx+nlinenc, y+idx);
	}
	
	for (i=end2; i<end1; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nlinenc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
	}
	
	i=end1;
	idx = i*nc;
	idx1 = idx-nc;
	idx2 = idx-nlinenc;
	matidx = idx*nc;
	matidx1 = idx1*nc;
	matidx2 = idx2*nc;
	smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
	smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
	smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_2D_nc5 (REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 5 banded 2D structured matrix (nc = 5)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsets of the five bands have to be (-1, +1, -nx, +nx) for nx != 1 and (-1,+1,-ny,+ny) 
 *       for nx = 1, but the order can be arbitrary. 
 */
static inline void spaaxpy_str_2D_nc5(REAL alpha, dSTRmat *A, REAL *x, REAL *y)
{
	INT i;
	INT idx,idx1,idx2;
	INT matidx, matidx1, matidx2;
	INT end1, end2;
	INT nline, nlinenc;
	
	//! information of A
	INT nx = A->nx;
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, *offdiag3=NULL;
	
	if (nx == 1)
	{
		nline = A->ny;
	}
	else
	{
		nline = nx;
	}
	nlinenc = nline*nc;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nline)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nline)
		{
			offdiag3 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 2D of nc=5 case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nline;
	
	smat_amxv_nc5(alpha, diag, x, y);
	smat_amxv_nc5(alpha, offdiag1, x+nc, y);
	smat_amxv_nc5(alpha, offdiag3, x+nlinenc, y);
	
	for (i=1; i<nline; ++i){
		idx = i*nc;
		matidx = idx*nc;
		idx1 = idx - nc;
		matidx1 = idx1*nc;
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc5(alpha, offdiag3+matidx, x+idx+nlinenc, y+idx);
	}
	
	
	for (i=nx; i<end2; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nlinenc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc5(alpha, offdiag3+matidx, x+idx+nlinenc, y+idx);
	}
	
	for (i=end2; i<end1; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nlinenc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
	}
	
	i=end1;
	idx = i*nc;
	idx1 = idx-nc;
	idx2 = idx-nlinenc;
	matidx = idx*nc;
	matidx1 = idx1*nc;
	matidx2 = idx2*nc;
	smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
	smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
	smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_2D_block (REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 5 banded 2D structured matrix (nc != 1)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsets of the five bands have to be (-1, +1, -nx, +nx) for nx != 1 and (-1,+1,-ny,+ny)
 *       for nx = 1, but the order can be arbitrary. 
 */
static inline void spaaxpy_str_2D_block (REAL alpha, 
                                         dSTRmat *A, 
                                         REAL *x, 
                                         REAL *y)
{
	INT i;
	INT idx,idx1,idx2;
	INT matidx, matidx1, matidx2;
	INT end1, end2;
	INT nline, nlinenc;
	
	//! information of A
	INT nx = A->nx;
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, *offdiag3=NULL;
	
	if (nx == 1)
	{
		nline = A->ny;
	}
	else
	{
		nline = nx;
	}
	nlinenc = nline*nc;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nline)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nline)
		{
			offdiag3 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 2D block case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nline;
	
	smat_amxv(alpha, diag, x, nc, y);
	smat_amxv(alpha, offdiag1, x+nc, nc, y);
	smat_amxv(alpha, offdiag3, x+nlinenc, nc, y);
	
	for (i=1; i<nline; ++i){
		idx = i*nc;
		matidx = idx*nc;
		idx1 = idx - nc;
		matidx1 = idx1*nc;
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
		smat_amxv(alpha, offdiag3+matidx, x+idx+nlinenc, nc, y+idx);
	}
	
	
	for (i=nx; i<end2; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nlinenc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
		smat_amxv(alpha, offdiag3+matidx, x+idx+nlinenc, nc, y+idx);
	}
	
	for (i=end2; i<end1; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nlinenc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
	}
	
	i=end1;
	idx = i*nc;
	idx1 = idx-nc;
	idx2 = idx-nlinenc;
	matidx = idx*nc;
	matidx1 = idx1*nc;
	matidx2 = idx2*nc;
	smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
	smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
	smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_3D_scalar(REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 7 banded 3D structured matrix (nc = 1)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsetsoffsets of the five bands have to be -1, +1, -nx, +nx, -nxy and +nxy, but the order 
 *       can be arbitrary. 
 */
static inline void spaaxpy_str_3D_scalar (REAL alpha, 
                                          dSTRmat *A, 
                                          REAL *x, 
                                          REAL *y)
{
	INT i;
	INT idx1,idx2,idx3;
	INT end1, end2, end3;
	//! information of A
	INT nx = A->nx;
	INT nxy = A->nxy;
	INT ngrid = A->ngrid;  //! number of grids
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, 
	*offdiag3=NULL, *offdiag4=NULL, *offdiag5=NULL;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nx)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nx)
		{
			offdiag3 = A->offdiag[i];
		}
		else if (A->offsets[i] == -nxy)
		{
			offdiag4 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nxy)
		{
			offdiag5 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 3D scalar case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nx;
	end3 = ngrid-nxy;
	
	y[0] += alpha*(diag[0]*x[0] + offdiag1[0]*x[1] + offdiag3[0]*x[nx] + offdiag5[0]*x[nxy]);
	
	for (i=1; i<nx; ++i){
		idx1 = i-1;
		y[i] += alpha*(offdiag0[idx1]*x[idx1] + diag[i]*x[i] + offdiag1[i]*x[i+1] + 
                       offdiag3[i]*x[i+nx] + offdiag5[i]*x[i+nxy]);
	}
	
	for (i=nx; i<nxy; ++i){
		idx1 = i-1; 
		idx2 = i-nx; 
		y[i] += alpha*(offdiag2[idx2]*x[idx2] + offdiag0[idx1]*x[idx1] + diag[i]*x[i] + 
                       offdiag1[i]*x[i+1] + offdiag3[i]*x[i+nx] + offdiag5[i]*x[i+nxy]);
	}
	
	for (i=nxy; i<end3; ++i){
		idx1 = i-1; 
		idx2 = i-nx; 
		idx3 = i-nxy;
		y[i] += alpha*(offdiag4[idx3]*x[idx3] + offdiag2[idx2]*x[idx2] + offdiag0[idx1]*x[idx1] + 
                       diag[i]*x[i] + offdiag1[i]*x[i+1] + offdiag3[i]*x[i+nx] + offdiag5[i]*x[i+nxy]);
	}
	
	for (i=end3; i<end2; ++i){
		idx1 = i-1; 
		idx2 = i-nx; 
		idx3 = i-nxy;
		y[i] += alpha*(offdiag4[idx3]*x[idx3] + offdiag2[idx2]*x[idx2] + offdiag0[idx1]*x[idx1] + 
                       diag[i]*x[i] + offdiag1[i]*x[i+1] + offdiag3[i]*x[i+nx]);
	}
	
	for (i=end2; i<end1; ++i){
		idx1 = i-1; 
		idx2 = i-nx; 
		idx3 = i-nxy;
		y[i] += alpha*(offdiag4[idx3]*x[idx3] + offdiag2[idx2]*x[idx2] + offdiag0[idx1]*x[idx1] + 
                       diag[i]*x[i] + offdiag1[i]*x[i+1]);
	}
	
	idx1 = end1-1; 
	idx2 = end1-nx; 
	idx3 = end1-nxy;
	y[end1] += alpha*(offdiag4[idx3]*x[idx3] + offdiag2[idx2]*x[idx2] + 
                      offdiag0[idx1]*x[idx1] + diag[end1]*x[end1]);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_3D_nc3(REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 7 banded 3D structured matrix (nc = 3)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsetsoffsets of the five bands have to be -1, +1, -nx, +nx, -nxyand +nxy, but the order 
 *       can be arbitrary. 
 */
static inline void spaaxpy_str_3D_nc3 (REAL alpha, 
                                       dSTRmat *A, 
                                       REAL *x, 
                                       REAL *y)
{
	INT i;
	INT idx,idx1,idx2,idx3;
	INT matidx, matidx1, matidx2, matidx3;
	INT end1, end2, end3;
	//! information of A
	INT nx = A->nx;
	INT nxy = A->nxy;
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;
	INT nxnc = nx*nc;
	INT nxync = nxy*nc;
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, 
	*offdiag3=NULL, *offdiag4=NULL, *offdiag5=NULL;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nx)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nx)
		{
			offdiag3 = A->offdiag[i];
		}
		else if (A->offsets[i] == -nxy)
		{
			offdiag4 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nxy)
		{
			offdiag5 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 3D of nc=3 case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nx;
	end3 = ngrid-nxy;
	
	smat_amxv_nc3(alpha, diag, x, y);
	smat_amxv_nc3(alpha, offdiag1, x+nc, y);
	smat_amxv_nc3(alpha, offdiag3, x+nxnc, y);
	smat_amxv_nc3(alpha, offdiag5, x+nxync, y);
	
	for (i=1; i<nx; ++i){
		idx = i*nc;
		matidx = idx*nc;
		idx1 = idx - nc;
		matidx1 = idx1*nc;
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc3(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
		smat_amxv_nc3(alpha, offdiag5+matidx, x+idx+nxync, y+idx);
	}
	
	for (i=nx; i<nxy; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc3(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
		smat_amxv_nc3(alpha, offdiag5+matidx, x+idx+nxync, y+idx);
		
	}
	
	for (i=nxy; i<end3; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv_nc3(alpha, offdiag4+matidx3, x+idx3, y+idx);
		smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc3(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
		smat_amxv_nc3(alpha, offdiag5+matidx, x+idx+nxync, y+idx);
	}
	
	for (i=end3; i<end2; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv_nc3(alpha, offdiag4+matidx3, x+idx3, y+idx);
		smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc3(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
	}
	
	for (i=end2; i<end1; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv_nc3(alpha, offdiag4+matidx3, x+idx3, y+idx);
		smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc3(alpha, offdiag1+matidx, x+idx+nc, y+idx);
	}
	
	i=end1;
	idx = i*nc;
	idx1 = idx-nc;
	idx2 = idx-nxnc;
	idx3 = idx-nxync;
	matidx = idx*nc;
	matidx1 = idx1*nc;
	matidx2 = idx2*nc;
	matidx3 = idx3*nc;
	smat_amxv_nc3(alpha, offdiag4+matidx3, x+idx3, y+idx);
	smat_amxv_nc3(alpha, offdiag2+matidx2, x+idx2, y+idx);
	smat_amxv_nc3(alpha, offdiag0+matidx1, x+idx1, y+idx);
	smat_amxv_nc3(alpha, diag+matidx, x+idx, y+idx);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_3D_nc5(REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 7 banded 3D structured matrix (nc = 5)
 * 
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsetsoffsets of the five bands have to be -1, +1, -nx, +nx, -nxyand +nxy, but the order 
 *       can be arbitrary. 
 */
static inline void spaaxpy_str_3D_nc5 (REAL alpha, 
                                       dSTRmat *A, 
                                       REAL *x, 
                                       REAL *y)
{
	INT i;
	INT idx,idx1,idx2,idx3;
	INT matidx, matidx1, matidx2, matidx3;
	INT end1, end2, end3;
	//! information of A
	INT nx = A->nx;
	INT nxy = A->nxy;
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;
	INT nxnc = nx*nc;
	INT nxync = nxy*nc;
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, 
	*offdiag3=NULL, *offdiag4=NULL, *offdiag5=NULL;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nx)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nx)
		{
			offdiag3 = A->offdiag[i];
		}
		else if (A->offsets[i] == -nxy)
		{
			offdiag4 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nxy)
		{
			offdiag5 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 3D of nc=5 case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nx;
	end3 = ngrid-nxy;
	
	smat_amxv_nc5(alpha, diag, x, y);
	smat_amxv_nc5(alpha, offdiag1, x+nc, y);
	smat_amxv_nc5(alpha, offdiag3, x+nxnc, y);
	smat_amxv_nc5(alpha, offdiag5, x+nxync, y);
	
	for (i=1; i<nx; ++i){
		idx = i*nc;
		matidx = idx*nc;
		idx1 = idx - nc;
		matidx1 = idx1*nc;
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc5(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
		smat_amxv_nc5(alpha, offdiag5+matidx, x+idx+nxync, y+idx);
	}
	
	for (i=nx; i<nxy; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc5(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
		smat_amxv_nc5(alpha, offdiag5+matidx, x+idx+nxync, y+idx);
		
	}
	
	for (i=nxy; i<end3; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv_nc5(alpha, offdiag4+matidx3, x+idx3, y+idx);
		smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc5(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
		smat_amxv_nc5(alpha, offdiag5+matidx, x+idx+nxync, y+idx);
	}
	
	for (i=end3; i<end2; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv_nc5(alpha, offdiag4+matidx3, x+idx3, y+idx);
		smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
		smat_amxv_nc5(alpha, offdiag3+matidx, x+idx+nxnc, y+idx);
	}
	
	for (i=end2; i<end1; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv_nc5(alpha, offdiag4+matidx3, x+idx3, y+idx);
		smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
		smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
		smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
		smat_amxv_nc5(alpha, offdiag1+matidx, x+idx+nc, y+idx);
	}
	
	i=end1;
	idx = i*nc;
	idx1 = idx-nc;
	idx2 = idx-nxnc;
	idx3 = idx-nxync;
	matidx = idx*nc;
	matidx1 = idx1*nc;
	matidx2 = idx2*nc;
	matidx3 = idx3*nc;
	smat_amxv_nc5(alpha, offdiag4+matidx3, x+idx3, y+idx);
	smat_amxv_nc5(alpha, offdiag2+matidx2, x+idx2, y+idx);
	smat_amxv_nc5(alpha, offdiag0+matidx1, x+idx1, y+idx);
	smat_amxv_nc5(alpha, diag+matidx, x+idx, y+idx);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_3D_block (REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where A is a 7 banded 3D structured matrix (nc != 1)
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 * \note the offsetsoffsets of the five bands have to be -1, +1, -nx, +nx, -nxyand +nxy, but the order 
 *       can be arbitrary. 
 */
static inline void spaaxpy_str_3D_block (REAL alpha, 
                                         dSTRmat *A, 
                                         REAL *x, 
                                         REAL *y)
{
	INT i;
	INT idx,idx1,idx2,idx3;
	INT matidx, matidx1, matidx2, matidx3;
	INT end1, end2, end3;
	//! information of A
	INT nx = A->nx;
	INT nxy = A->nxy;
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;
	INT nxnc = nx*nc;
	INT nxync = nxy*nc;
	INT nband = A->nband;
	
	REAL *diag = A->diag;
	REAL *offdiag0=NULL, *offdiag1=NULL, *offdiag2=NULL, 
	*offdiag3=NULL, *offdiag4=NULL, *offdiag5=NULL;
	
	for (i=0; i<nband; ++i)
	{
		if (A->offsets[i] == -1)
		{
			offdiag0 = A->offdiag[i];
		}
		else if (A->offsets[i] == 1)
		{
			offdiag1 = A->offdiag[i];
		}	
		else if (A->offsets[i] == -nx)
		{
			offdiag2 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nx)
		{
			offdiag3 = A->offdiag[i];
		}
		else if (A->offsets[i] == -nxy)
		{
			offdiag4 = A->offdiag[i];
		}		
		else if (A->offsets[i] == nxy)
		{
			offdiag5 = A->offdiag[i];
		}		
		else
		{
			printf("### WARNING: offsets for 3D block case is illegal!\n");
			spaaxpy_str_general(alpha, A, x, y); 
			return;
		}
	}
	
	end1 = ngrid-1;
	end2 = ngrid-nx;
	end3 = ngrid-nxy;
	
	smat_amxv(alpha, diag, x, nc, y);
	smat_amxv(alpha, offdiag1, x+nc, nc, y);
	smat_amxv(alpha, offdiag3, x+nxnc, nc, y);
	smat_amxv(alpha, offdiag5, x+nxync, nc, y);
	
	for (i=1; i<nx; ++i){
		idx = i*nc;
		matidx = idx*nc;
		idx1 = idx - nc;
		matidx1 = idx1*nc;
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
		smat_amxv(alpha, offdiag3+matidx, x+idx+nxnc, nc, y+idx);
		smat_amxv(alpha, offdiag5+matidx, x+idx+nxync, nc, y+idx);
	}
	
	for (i=nx; i<nxy; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
		smat_amxv(alpha, offdiag3+matidx, x+idx+nxnc, nc, y+idx);
		smat_amxv(alpha, offdiag5+matidx, x+idx+nxync, nc, y+idx);
		
	}
	
	for (i=nxy; i<end3; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv(alpha, offdiag4+matidx3, x+idx3, nc, y+idx);
		smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
		smat_amxv(alpha, offdiag3+matidx, x+idx+nxnc, nc, y+idx);
		smat_amxv(alpha, offdiag5+matidx, x+idx+nxync, nc, y+idx);
	}
	
	for (i=end3; i<end2; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv(alpha, offdiag4+matidx3, x+idx3, nc, y+idx);
		smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
		smat_amxv(alpha, offdiag3+matidx, x+idx+nxnc, nc, y+idx);
	}
	
	for (i=end2; i<end1; ++i){
		idx = i*nc;
		idx1 = idx-nc;
		idx2 = idx-nxnc;
		idx3 = idx-nxync;
		matidx = idx*nc;
		matidx1 = idx1*nc;
		matidx2 = idx2*nc;
		matidx3 = idx3*nc;
		smat_amxv(alpha, offdiag4+matidx3, x+idx3, nc, y+idx);
		smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
		smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
		smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
		smat_amxv(alpha, offdiag1+matidx, x+idx+nc, nc, y+idx);
	}
	
	i=end1;
	idx = i*nc;
	idx1 = idx-nc;
	idx2 = idx-nxnc;
	idx3 = idx-nxync;
	matidx = idx*nc;
	matidx1 = idx1*nc;
	matidx2 = idx2*nc;
	matidx3 = idx3*nc;
	smat_amxv(alpha, offdiag4+matidx3, x+idx3, nc, y+idx);
	smat_amxv(alpha, offdiag2+matidx2, x+idx2, nc, y+idx);
	smat_amxv(alpha, offdiag0+matidx1, x+idx1, nc, y+idx);
	smat_amxv(alpha, diag+matidx, x+idx, nc, y+idx);
	
	return;
	
}

/**
 * \fn void spaaxpy_str_general(REAL alpha, dSTRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y for general cases
 *
 * \param alpha real number
 * \param A pointer to dSTRmat matrix
 * \param x pointer to real array
 * \param y pointer to real array
 *
 */
static inline void spaaxpy_str_general (REAL alpha, 
                                        dSTRmat *A, 
                                        REAL *x, 
                                        REAL *y) 
{
	//! information of A
	INT ngrid = A->ngrid;  //! number of grids
	INT nc = A->nc;        //! size of each block (number of components)
	INT nband = A->nband ; //! number of off-diag band
	INT *offsets = A->offsets; //! offsets of the off-diagals
	REAL  *diag = A->diag;       //! Diagonal entries
	REAL **offdiag = A->offdiag; //! Off-diagonal entries
	
	//! local variables
	INT k;
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT size = nc*ngrid; 
	INT nc2  = nc*nc;
	INT ncw  = 0;
	INT start_data = 0;
	INT start_vecx = 0;
	INT start_vecy = 0;
	INT start_vect = 0;
	REAL beta = 0.0;
	
	if (alpha == 0)
	{
		return; //! nothing should be done
	}
	
	beta = 1.0/alpha;
	
	//! y: = beta*y
	for (k = 0; k < size; ++k)
	{
		y[k] *= beta;
	}
	
	//! y: = y + A*x   
	if (nc > 1)
	{	
		//! Deal with the diagonal band
		for (block = 0; block < ngrid; ++block)
		{
			start_data = nc2*block;
			start_vect = nc*block;
			blkcontr_str(start_data,start_vect,start_vect,nc,diag,x,y);
		}
		
		//! Deal with the off-diagonal bands
		for (band = 0; band < nband; band ++)
		{  
			width = offsets[band];
			ncw   = nc*width;
			if (width < 0)
			{
				for (block = 0; block < ngrid+width; ++block)
				{
					start_data = nc2*block;
					start_vecx = nc*block;
					start_vecy = start_vecx - ncw; 
					blkcontr_str(start_data,start_vecx,start_vecy,nc,offdiag[band],x,y);
				}
			}
			else
			{
				for (block = 0; block < ngrid-width; ++block)
				{
					start_data = nc2*block;
					start_vecy = nc*block;
					start_vecx = start_vecy + ncw;
					blkcontr_str(start_data,start_vecx,start_vecy,nc,offdiag[band],x,y);
				}
			}
		}
	}
	else if (nc == 1)
	{
		//! Deal with the diagonal band
		for (point = 0; point < ngrid; point ++)
		{
			y[point] += diag[point]*x[point];
		}
		
		//! Deal with the off-diagonal bands
		for (band = 0; band < nband; band ++)
		{  
			width = offsets[band];
			if (width < 0)
			{
				for (point = 0; point < ngrid+width; point ++)
				{
					y[point-width] += offdiag[band][point]*x[point];
				}
			}
			else
			{
				for (point = 0; point < ngrid-width; point ++)
				{  
					y[point] += offdiag[band][point]*x[point+width];
				}
			}
		}
	}
	else
	{
		printf("### WARNING: nc is illegal!\n");
		return;
	}
	
	//! y: = alpha*y
	for (k = 0; k < size; ++k)
	{
		y[k] *= alpha;
	}        	
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
