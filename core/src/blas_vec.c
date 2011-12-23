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
 * \fn void fasp_blas_dvec_axpy (const double a, dvector *x, dvector *y) 
 * \brief y = a*x + y
 * \param a real number
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_dvec_axpy (const double a, 
                          dvector *x, 
                          dvector *y)
{
	unsigned int i, m=x->row;
	double *xpt=x->val, *ypt=y->val;
	
	if ((y->row-m)!=0) {
		printf("Error: two vectors have different length!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	for (i=0; i<m; ++i) ypt[i] += a*xpt[i];
}

/**
 * \fn void fasp_blas_dvec_axpyz (const double a, dvector *x, dvector *y, dvector *z) 
 * \brief z = a*x + y, z is a third vector (z is cleared)
 * \param a real number
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param *z pointer to dvector
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
void fasp_blas_dvec_axpyz (const double a, 
                           dvector *x, 
                           dvector *y, 
                           dvector *z)
{
	unsigned int i;
	const int m=x->row;
	double *xpt=x->val, *ypt=y->val, *zpt=z->val;
	
	if ((y->row-m)!=0) {
		printf("Error: two vectors have different length!\n");
		exit(ERROR_DATA_STRUCTURE);
	}
	
	z->row = m;
	memcpy(zpt,ypt,m*sizeof(double));
	for (i=0; i<m; ++i) zpt[i] += a*xpt[i];
}

/**
 * \fn double fasp_blas_dvec_dotprod (dvector *x, dvector *y) 
 * \brief Inner product of two vectors (x,y)
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return Inner product
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
double fasp_blas_dvec_dotprod (dvector *x, 
                               dvector *y)
{
	unsigned int i;
	const int length=x->row;
	double *xpt=x->val, *ypt=y->val;	
	double value=0;
	for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
	return value;
}


/**
 * \fn double fasp_dvec_relerr (dvector *x, dvector *y) 
 * \brief Relative error of two dvector x and y
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \return relative error ||x-y||/||x||
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
double fasp_dvec_relerr (dvector *x, 
                         dvector *y)
{
	unsigned int i;
	const int length=x->row;
	double diff=0, temp=0;
	double *xpt=x->val, *ypt=y->val;
	
	if (length!=y->row) {
		printf("Error: The lengths of vectors do not match! \n");
		exit(ERROR_DUMMY_VAR);	
	}
	
	for (i=0;i<length;++i) {
		temp += xpt[i]*xpt[i];
		diff += pow(xpt[i]-ypt[i],2);
	}
	
	return sqrt(diff/temp);
}

/**
 * \fn double fasp_blas_dvec_norm1 (dvector *x) 
 * \brief L1 norm of dvector x
 * \param *x pointer to dvector
 * \return L1 norm of x
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
double fasp_blas_dvec_norm1 (dvector *x)
{
	unsigned int i;
	const int length=x->row;
	double *xpt=x->val;	
	double onenorm=0;
	for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
	return onenorm;
}

/**
 * \fn double fasp_blas_dvec_norm2 (dvector *x) 
 * \brief L2 norm of dvector x
 * \param *x pointer to dvector
 * \return L2 norm of x
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
double fasp_blas_dvec_norm2 (dvector *x)
{
	unsigned int i;
	const int length=x->row;
	double *xpt=x->val;	
	double twonorm=0;
	for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
	return sqrt(twonorm);
}

/**
 * \fn double fasp_blas_dvec_norminf (dvector *x) 
 * \brief Linf norm of dvector x
 * \param *x pointer to dvector
 * \return Linf norm of x
 *
 * \author Chensong Zhang
 * \date 07/01/209
 */
double fasp_blas_dvec_norminf (dvector *x)
{
	unsigned int i;
	const int length=x->row;
	double *xpt=x->val;	
	double infnorm=0;
	for (i=0;i<length;++i) infnorm=MAX(infnorm,ABS(xpt[i]));
	return infnorm;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
