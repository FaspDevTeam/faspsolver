/*! \file blas_bsr.c
 *  \brief BLAS operations for sparse matrices in BSR format.
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/*!
 * \fn void fasp_blas_dbsr_aAxpby ( const REAL alpha, dBSRmat *A, 
 *                                  REAL *x, const REAL beta, REAL *y )
 *
 * \brief Compute y := alpha*A*x + beta*y
 *
 * \param alpha   a real number
 * \param *A      pointer to the matrix
 * \param *x      pointer to the vector x
 * \param beta    a real number
 * \param *y      pointer to the vector y
 *
 * \author Zhou Zhiyang
 * \date 2010/10/25  
 *
 * \note: works for general nb (Xiaozhe)
 */
void fasp_blas_dbsr_aAxpby (const REAL alpha, 
                            dBSRmat *A, 
                            REAL *x, 
                            const REAL beta, 
                            REAL *y )
{
	/* members of A */
	INT     ROW = A->ROW;
	INT     nb  = A->nb;
	INT    *IA  = A->IA;
	INT    *JA  = A->JA;
	REAL   *val = A->val;
	
	/* local variables */
	INT     size = ROW*nb;
	INT     jump = nb*nb;
	INT     i,j,k,ibegin,iend;
	REAL    temp;
	REAL   *pA  = NULL;
	REAL   *px0 = NULL;
	REAL   *py0 = NULL;
	REAL   *py  = NULL;
	
	//----------------------------------------------
	//   Treat (alpha == 0.0) computation 
	//----------------------------------------------
	
	if (alpha == 0.0)
	{
		for (i = size; i--; ) y[i] *= beta;
		return;
	}
	
	//-------------------------------------------------
	//   y = (beta/alpha)*y
	//-------------------------------------------------
	
	temp = beta / alpha;
	if (temp != 1.0)
	{
		if (temp == 0.0)
		{
			memset(y, 0X0, size*sizeof(REAL));
		}
		else
		{
			for (i = size; i--; ) y[i] *= temp;  // modified by Xiaozhe, 03/11/2011
		}
	}
	
	//-----------------------------------------------------------------
	//   y += A*x (Core Computation)
	//   each non-zero block elements are stored in row-major order
	//-----------------------------------------------------------------
	
	for (i = 0; i < ROW; ++i)
	{
		py0 = &y[i*nb];
		ibegin = IA[i]; iend = IA[i+1];
		for (k = ibegin; k < iend; ++k)
		{
			j = JA[k];
			pA = val+k*jump; // &val[k*jump];
			px0 = x+j*nb; // &x[j*nb];
			py = py0;
			fasp_blas_smat_ypAx( pA, px0, py, nb );
		}  
	}
	
	//------------------------------------------
	//   y = alpha*y
	//------------------------------------------
	
	if (alpha != 1.0)
	{
		for (i = size; i--; ++i)
		{
			y[i] *= alpha;
		}
	}   
}


/*!
 * \fn void fasp_blas_dbsr_aAxpy ( const REAL alpha, dBSRmat *A, REAL *x, REAL *y )
 *
 * \brief Compute y := alpha*A*x + y
 *
 * \param alpha a real number
 * \param *A pointer to the matrix
 * \param *x pointer to the vector x
 * \param *y pointer to the vector y
 *
 * \author Zhou Zhiyang
 * \date 2010/10/25  
 *
 * \note: works for general nb (Xiaozhe)
 */
void fasp_blas_dbsr_aAxpy (const REAL alpha, 
                           dBSRmat *A, 
                           REAL *x, 
                           REAL *y)
{
	/* members of A */
	INT     ROW = A->ROW;
	INT     nb  = A->nb;
	INT    *IA  = A->IA;
	INT    *JA  = A->JA;
	REAL   *val = A->val;
	
	/* local variables */
	INT     size = ROW*nb;
	INT     jump = nb*nb;
	INT     i,j,k, ibegin, iend;
	REAL    temp = 0.0;
	REAL    *pA   = NULL;
	REAL    *px0  = NULL;
	REAL    *py0  = NULL;
	REAL    *py   = NULL;
	
	//----------------------------------------------
	//   Treat (alpha == 0.0) computation 
	//----------------------------------------------
	
	if (alpha == 0.0)
	{
		return; // Nothting to compute
	}
	
	//-------------------------------------------------
	//   y = (1.0/alpha)*y
	//-------------------------------------------------
	
	if (alpha != 1.0)
	{
		temp = 1.0 / alpha;
		for (i = size; i-- ; ) y[i] *= temp; // modified by Xiaozhe, 03/11/2011
	}
	
	//-----------------------------------------------------------------
	//   y += A*x (Core Computation)
	//   each non-zero block elements are stored in row-major order
	//-----------------------------------------------------------------
	
	switch (nb)
	{
		case 3: 
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*3];
				ibegin = IA[i]; iend = IA[i+1]; // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*9; // &val[k*jump];
					px0 = x+j*3; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx_nc3( pA, px0, py );
				}  
			}
		}
			break;
			
		case 5:
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*5];
				ibegin = IA[i]; iend = IA[i+1]; // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*25; // &val[k*jump];
					px0 = x+j*5; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx_nc5( pA, px0, py );
				}  
			}
		}
			break;
			
		case 7:
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*7];
				ibegin = IA[i]; iend = IA[i+1]; // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*49; // &val[k*jump];
					px0 = x+j*7; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx_nc7( pA, px0, py );
				}  
			}
		}
			break;
			
		default: 
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*nb];
				ibegin = IA[i]; iend = IA[i+1]; // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*jump; // &val[k*jump];
					px0 = x+j*nb; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx( pA, px0, py, nb );
				}  
			}
		}
			break;
			
	}
	
	//------------------------------------------
	//   y = alpha*y
	//------------------------------------------
	
	if (alpha != 1.0)
	{
		for (i = size; i--; ) // modified by Xiaozhe, 03/11/2011
		{
			y[i] *= alpha;
		}
	} 
	
	return;
}

/*!
 * \fn void fasp_blas_dbsr_mxv ( dBSRmat *A, REAL *x, REAL *y )
 *
 * \brief Compute y := A*x
 *
 * \param *A pointer to the matrix
 * \param *x pointer to the vector x
 * \param *y pointer to the vector y
 *
 * \author Zhou Zhiyang
 * \date 2010/10/25 
 *
 * \note: works for general nb (Xiaozhe) 
 */
void fasp_blas_dbsr_mxv (dBSRmat *A, 
                         REAL *x, 
                         REAL *y)
{
	/* members of A */
	INT     ROW = A->ROW;
	INT     nb  = A->nb;
	INT    *IA  = A->IA;
	INT    *JA  = A->JA;
	REAL   *val = A->val;
	
	/* local variables */
	INT     size = ROW*nb;
	INT     jump = nb*nb;
	INT     i,j,k, ibegin, iend;
	REAL   *pA  = NULL;
	REAL   *px0 = NULL;
	REAL   *py0 = NULL;
	REAL   *py  = NULL;
	
	//-----------------------------------------------------------------
	//  zero out 'y' 
	//-----------------------------------------------------------------
	
	memset(y, 0X0, size*sizeof(REAL));
	
	//-----------------------------------------------------------------
	//   y = A*x (Core Computation)
	//   each non-zero block elements are stored in row-major order
	//-----------------------------------------------------------------
	
	switch (nb)
	{
		case 3: 
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*3];
				ibegin = IA[i]; iend = IA[i+1];  // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*9; // &val[k*jump];
					px0 = x+j*3; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx_nc3( pA, px0, py );
				}  
			}
		}
			break;
			
		case 5:
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*5];
				ibegin = IA[i]; iend = IA[i+1];  // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*25; // &val[k*jump];
					px0 = x+j*5; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx_nc5( pA, px0, py );
				}  
			}
		}
			break;
			
		case 7:
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*7];
				ibegin = IA[i]; iend = IA[i+1];  // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*49; // &val[k*jump];
					px0 = x+j*7; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx_nc7( pA, px0, py );
				}  
			}
		}
			break;
			
		default: 
		{
			for (i = 0; i < ROW; ++i)
			{
				py0 = &y[i*nb];
				ibegin = IA[i]; iend = IA[i+1]; // modified by Xiaozhe, 03/11/2011
				for (k = ibegin; k < iend; ++k)
				{
					j = JA[k];
					pA = val+k*jump; // &val[k*jump];
					px0 = x+j*nb; // &x[j*nb];
					py = py0;
					fasp_blas_smat_ypAx( pA, px0, py, nb );
				}  
			}
		}
			break;
			
	}
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
