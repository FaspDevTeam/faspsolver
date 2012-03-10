/*! \file blas_array.c
 *  \brief BLAS operations for arrays.
 * 
 *  Some simple array operations ...
 *
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_array_ax (const INT n, const REAL a, REAL *x)
 *
 * \brief x = a*x
 *
 * \param n    number of variables
 * \param a    a real number
 * \param x   pointer to the original vector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_array_ax (const INT n, 
                         const REAL a, 
                         REAL *x)
{
	unsigned INT i;
	
	if (a==1.0) {
	}
	else {
		for (i=0; i<n; ++i) x[i] = a*x[i];
	}
}

/**
 * \fn void fasp_blas_array_axpy (const INT n, const REAL a, REAL *x, REAL *y)
 *
 * \brief y = a*x + y
 *
 * \param n   number of variables
 * \param a   a real number
 * \param x  pointer to the original vector
 * \param y  pointer to the destination vector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_array_axpy (const INT n, 
                           const REAL a, 
                           REAL *x, 
                           REAL *y)
{
	unsigned INT i;
	
	if (a==1.0) {
		for (i=0; i<n; ++i) y[i] += x[i];
	}
	else if (a==-1.0) {
		for (i=0; i<n; ++i) y[i] -= x[i];
	}
	else {
		for (i=0; i<n; ++i) y[i] += a*x[i];
	}
}

/**
 * \fn void fasp_blas_array_axpyz(const INT n, const REAL a, REAL *x,
 *                                REAL *y, REAL *z)
 *
 * \brief z = a*x + y, z is the third vector
 *
 * \param n   number of variables
 * \param a   a real number
 * \param x  pointer to the original vector 1
 * \param y  pointer to the original vector 2
 * \param z  pointer to the destination vector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_array_axpyz (const INT n, 
                            const REAL a, 
                            REAL *x, 
                            REAL *y, 
                            REAL *z)
{
	unsigned INT i;
	for (i=0; i<n; ++i) z[i] = a*x[i]+y[i];
}

/**
 * \fn void fasp_blas_array_axpby (const INT n, const REAL a, REAL *x, 
 *                                 const REAL b, REAL *y)
 *
 * \brief y = a*x + b*y
 *
 * \param n   number of variables
 * \param a   real number
 * \param b   real number
 * \param x  pointer to the origianl vector
 * \param y  pointer to the destination vector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_array_axpby (const INT n, 
                            const REAL a, 
                            REAL *x, 
                            const REAL b, 
                            REAL *y)
{
	unsigned INT i;
	for (i=0; i<n; ++i) y[i] = a*x[i]+b*y[i];
}

/**
 * \fn REAL fasp_blas_array_dotprod (const INT n, REAL *x, REAL *y) 
 *
 * \brief Inner product of two arraies (x,y)
 *
 * \param n   number of variables
 * \param x  pointer to vector 1
 * \param y  pointer to vector 2
 *
 * \return    inner product
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
REAL fasp_blas_array_dotprod (const INT n, 
                              REAL *x, 
                              REAL *y)
{
	unsigned INT i; 
	REAL value = 0.0;
	
	for (i=n;i--;) value+=x[i]*y[i];    // modified by Xiaozhe Hu, 03/11/2011
	
	//const INT one=1;
	//value=ddot_(&n,x,&one,y,&one); /* requires common.h and cblas.h */
	
	return value;
}

/**
 * \fn REAL fasp_blas_array_norm1 (const INT n, REAL *x)
 *
 * \brief L1 norm of array x
 *
 * \param n   number of variables
 * \param x  pointer to the original vector
 *
 * \return    L1 norm of x
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
REAL fasp_blas_array_norm1 (const INT n, 
                            REAL *x)
{
	unsigned INT i;
	REAL onenorm=0;
	for (i=0;i<n;++i) onenorm+=ABS(x[i]);
	return onenorm;
}

/**
 * \fn REAL fasp_blas_array_norm2 (const INT n, REAL *x) 
 *
 * \brief L2 norm of array x
 *
 * \param n   number of variables
 * \param x  pointer to the original vector
 *
 * \return    L2 norm of x
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
REAL fasp_blas_array_norm2 (const INT n, 
                            REAL *x)
{
	unsigned INT i;
	REAL twonorm=0;
	for (i=n;i--;) twonorm+=x[i]*x[i];
	return sqrt(twonorm);
}

/**
 * \fn REAL fasp_blas_array_norminf (const INT n, REAL *x) 
 *
 * \brief Linf norm of array x
 *
 * \param n   number of variables
 * \param x  pointer to the original vector
 *
 * \return    Linf norm of x
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
REAL fasp_blas_array_norminf (const INT n, 
                              REAL *x)
{
	unsigned INT i;
	REAL infnorm=0;
	for (i=0;i<n;++i) infnorm=MAX(infnorm,ABS(x[i]));
	return infnorm;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
