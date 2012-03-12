/*! \file blas_vec.c
 *  \brief BLAS operations for vectors
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_dvec_axpy (const REAL a, dvector *x, dvector *y)
 *
 * \brief y = a*x + y
 *
 * \param a   REAL factor a
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
void fasp_blas_dvec_axpy (const REAL a, 
                          dvector *x, 
                          dvector *y)
{
	unsigned INT i, m=x->row;
	REAL *xpt=x->val, *ypt=y->val;
	
	if ((y->row-m)!=0) {
		printf("### ERROR: Two vectors have different length!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	for (i=0; i<m; ++i) ypt[i] += a*xpt[i];
}

/**
 * \fn void fasp_blas_dvec_axpyz (const REAL a, dvector *x, dvector *y, dvector *z) 
 *
 * \brief z = a*x + y, z is a third vector (z is cleared)
 *
 * \param a   REAL factor a
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 * \param z   Pointer to dvector z
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
void fasp_blas_dvec_axpyz (const REAL a, 
                           dvector *x, 
                           dvector *y, 
                           dvector *z)
{
	unsigned INT i;
	const INT m=x->row;
	REAL *xpt=x->val, *ypt=y->val, *zpt=z->val;
	
	if ((y->row-m)!=0) {
		printf("### ERROR: Two vectors have different length!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	z->row = m;
	memcpy(zpt,ypt,m*sizeof(REAL));
	for (i=0; i<m; ++i) zpt[i] += a*xpt[i];
}

/**
 * \fn REAL fasp_blas_dvec_dotprod (dvector *x, dvector *y) 
 *
 * \brief Inner product of two vectors (x,y)
 *
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \return Inner product
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
REAL fasp_blas_dvec_dotprod (dvector *x, 
                             dvector *y)
{
	unsigned INT i;
	const INT length=x->row;
	REAL *xpt=x->val, *ypt=y->val;	
	
    register REAL value=0;
	for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
	return value;
}


/**
 * \fn REAL fasp_dvec_relerr (dvector *x, dvector *y) 
 *
 * \brief Relative error of two dvector x and y
 *
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \return relative error ||x-y||/||x||
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
REAL fasp_dvec_relerr (dvector *x, 
                       dvector *y)
{
	unsigned INT i;
	const INT length=x->row;
	REAL diff=0, temp=0;
	REAL *xpt=x->val, *ypt=y->val;
	
	if (length!=y->row) {
		printf("### ERROR: Two vectors have different length!\n");
		exit(ERROR_DUMMY_VAR);	
	}
	
	for (i=0;i<length;++i) {
		temp += xpt[i]*xpt[i];
		diff += pow(xpt[i]-ypt[i],2);
	}
	
	return sqrt(diff/temp);
}

/**
 * \fn REAL fasp_blas_dvec_norm1 (dvector *x) 
 *
 * \brief L1 norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L1 norm of x
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
REAL fasp_blas_dvec_norm1 (dvector *x)
{
	unsigned INT i;
	const INT length=x->row;
	REAL *xpt=x->val;	
	
    register REAL onenorm=0;
	for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
	return onenorm;
}

/**
 * \fn REAL fasp_blas_dvec_norm2 (dvector *x) 
 *
 * \brief L2 norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L2 norm of x
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
REAL fasp_blas_dvec_norm2 (dvector *x)
{
	unsigned INT i;
	const INT length=x->row;
	REAL *xpt=x->val;	
	
    register REAL twonorm=0;
	for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
	return sqrt(twonorm);
}

/**
 * \fn REAL fasp_blas_dvec_norminf (dvector *x) 
 *
 * \brief Linf norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L_inf norm of x
 *
 * \author Chensong Zhang
 * \data   07/01/209
 */
REAL fasp_blas_dvec_norminf (dvector *x)
{
	unsigned INT i;
	const INT length=x->row;
	REAL *xpt=x->val;
    
	register REAL infnorm=0;
	for (i=0;i<length;++i) infnorm=MAX(infnorm,ABS(xpt[i]));
	return infnorm;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
