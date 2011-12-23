/*! \file array.c
 *  \brief Array operations
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
 * \fn void fasp_array_init (REAL *x)
 *
 * \brief Initialize an array
 *
 * \param *x pointer to the vector
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_array_init (REAL *x) 
{
	x=NULL;
}

/**
 * \fn void fasp_array_set (const INT n, REAL *x, const REAL val)
 *
 * \brief Set initial value for an array to be x=val
 *
 * \param n number of variables
 * \param *x pointer to the vector
 * \param val initial value for the REAL array
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_array_set (const INT n, 
                     REAL *x, 
                     const REAL val)
{
	unsigned INT i;
	for (i=0; i<n; ++i) x[i]=val;
}

/**
 * \fn void fasp_array_cp (const INT n, REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x
 *
 * \param n number of variables
 * \param *x pointer to the original vector
 * \param *y pointer to the destination vector
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_array_cp (const INT n, 
                    REAL *x, 
                    REAL *y)
{
	memcpy(y,x,n*sizeof(REAL));
}

/**
 * \fn void fasp_array_cp_nc3 (REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x, the length is 3
 *
 * \param *x   pointer to the original vector
 * \param *y   pointer to the destination vector
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date 05/01/2010
 *
 * \note Special unrolled routine designed for a specific application
 */
void fasp_array_cp_nc3 (REAL *x, 
                        REAL *y) 
{
	y[0] = x[0];
	y[1] = x[1];
	y[2] = x[2];
}

/**
 * \fn void fasp_array_cp_nc5(REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x, the length is 5
 *
 * \param *x   pointer to the original vector
 * \param *y   pointer to the destination vector
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date 05/01/2010
 *
 * \note Special unrolled routine designed for a specific application
 */
void fasp_array_cp_nc5 (REAL *x, 
                        REAL *y) 
{
	y[0] = x[0];
	y[1] = x[1];
	y[2] = x[2];
	y[3] = x[3];
	y[4] = x[4];
}

/**
 * \fn void fasp_array_cp_nc7(REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x, the length is 7
 *
 * \param *x   pointer to the original vector
 * \param *y   pointer to the destination vector
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date 05/01/2010
 *
 * \note Special unrolled routine designed for a specific application
 */
void fasp_array_cp_nc7 (REAL *x, 
                        REAL *y) 
{
	y[0] = x[0];
	y[1] = x[1];
	y[2] = x[2];
	y[3] = x[3];
	y[4] = x[4];
	y[5] = x[5];
	y[6] = x[6];
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

