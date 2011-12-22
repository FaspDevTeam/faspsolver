/*! \file blas_csrl.c
 *  \brief BLAS operations for sparse matrices in CSRL format.
 *
 *  \note For details of CSRL format, refer to 
 *      "Optimizaing sparse matrix vector product computations using unroll and jam"
 *      by John Mellor-Crummey and John Garvin, Tech Report Rice Univ, Aug 2002.
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_blas_dcsrl_mxv ( dCSRLmat *A, REAL *x, REAL *y )
 *
 * \brief Compute y = A*x for a sparse matrix in CSRL format
 *
 * \param *A pointer to the dCSRLmat type matrix
 * \param *x pointer to the REAL array of vector x
 * \param *y pointer to the REAL array of vector y
 *
 * \author Zhou Zhiyang
 * \date 2011/01/07
 */
INT fasp_blas_dcsrl_mxv (dCSRLmat *A, 
                         REAL *x, 
                         REAL *y)
{
	INT     dif      = A -> dif;
	INT    *nzdifnum = A -> nzdifnum;
	INT    *rowindex = A -> rowindex;
	INT    *rowstart = A -> rowstart;
	INT    *ja       = A -> ja;
	REAL   *a        = A -> data;
	
	INT i;
	INT row, col=0;
	INT len, rowlen;
	INT firstrow, lastrow;
	
	REAL val0, val1;
    
	for (len = 0; len < dif; len ++)
	{
		firstrow = rowstart[len];
		lastrow  = rowstart[len+1] - 1;
		rowlen   = nzdifnum[len];
		
		if (lastrow > firstrow )
		{
			//----------------------------------------------------------
			// Fully-unrolled code for special case (i.g.,rowlen = 5) 
			// Note: you can also set other special case
			//----------------------------------------------------------               
			if (rowlen == 5)
			{
				for (row = firstrow; row < lastrow; row += 2)
				{
					val0 = a[col]*x[ja[col]];
					val1 = a[col+5]*x[ja[col+5]];
					col ++;
					
					val0 += a[col]*x[ja[col]];
					val1 += a[col+5]*x[ja[col+5]];
					col ++;
					
					val0 += a[col]*x[ja[col]];
					val1 += a[col+5]*x[ja[col+5]];
					col ++;
					
					val0 += a[col]*x[ja[col]];
					val1 += a[col+5]*x[ja[col+5]];
					col ++;
					
					val0 += a[col]*x[ja[col]];
					val1 += a[col+5]*x[ja[col+5]];
					col ++;
					
					y[rowindex[row]] = val0;
					y[rowindex[row+1]] = val1;
					
					col += 5;
				}
			}
			else        
			{
				//------------------------------------------------------------------
				// Unroll-and-jammed code for handling two rows at a time 
				//------------------------------------------------------------------
				
				for (row = firstrow; row < lastrow; row += 2)
				{ 
					val0 = 0.0;
					val1 = 0.0;
					for (i = 0; i < rowlen; i ++)
					{
						val0 += a[col]*x[ja[col]];
						val1 += a[col+rowlen]*x[ja[col+rowlen]];
						col ++;
					}
					y[rowindex[row]] = val0;
					y[rowindex[row+1]] = val1;
					col += rowlen;               
				}  
			}
			firstrow = row;
		}
		
		//-----------------------------------------------------------
		// Handle leftover rows that can't be handled in bundles 
		// in the unroll-and -jammed loop 
		//-----------------------------------------------------------
		
		for (row = firstrow; row <= lastrow; row ++)
		{
			val0 = 0.0;
			for (i = 0; i < rowlen; i ++)
			{
				val0 += a[col]*x[ja[col]];
				col ++;
			}
			y[rowindex[row]] = val0;
		}
		
	}
	
	return 0;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
