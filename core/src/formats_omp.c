/*! \file formats_omp.c
 *  \brief Matrix format conversion routines
 *
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

#if FASP_USE_OPENMP

/*!
 * \fn dCSRmat dBSR2dCSRMatrix_omp( dBSRmat *B, int nthreads, int openmp_holds )
 * \brief Transfer a 'dBSRmat' type matrix into a dCSRmat.
 * \param B pointer to the 'dBSRmat' type matrix
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
dCSRmat fasp_format_dbsr_dcsr_omp (dBSRmat *B, 
                                   int nthreads, 
                                   int openmp_holds )
{
	dCSRmat A;
    
	/* members of B */
	int     ROW = B->ROW;
	int     COL = B->COL;
	int     NNZ = B->NNZ;    
	int     nb  = B->nb;
	int    *IA  = B->IA;
	int    *JA  = B->JA;
	int     storage_manner = B->storage_manner;
	double *val = B->val;
	
	int jump = nb*nb;
	int rowA = ROW*nb;
	int colA = COL*nb;
	int nzA  = NNZ*jump;
	
	int     *ia = NULL;
	int     *ja = NULL;
	double  *a  = NULL;
	
	int i,j,k;
	int mr,mc;
	int rowstart0,rowstart,colstart0,colstart;
	int colblock,nzperrow; 
	double *vp = NULL;
	double *ap = NULL;
	int    *jap = NULL;
	int stride_i,mybegin,myend,myid;
    
	//--------------------------------------------------------
	// Create a CSR Matrix 
	//--------------------------------------------------------
	A = fasp_dcsr_create(rowA, colA, nzA);
	ia = A.IA;
	ja = A.JA;
	a = A.val;
	
	//--------------------------------------------------------------------------
	// Compute the number of nonzeros per row, and after this loop,
	// ia[i],i=1:rowA, will be the number of nonzeros of the (i-1)-th row.
	//--------------------------------------------------------------------------
	
	if (ROW > openmp_holds) {
		stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i, rowstart, colblock, nzperrow, j) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			mybegin = myid*stride_i;
			if(myid < nthreads-1) myend = mybegin+stride_i;
			else myend = ROW;
			for (i=mybegin; i<myend; ++i)
			{
				rowstart = i*nb + 1;
				colblock = IA[i+1] - IA[i];
				nzperrow = colblock*nb;
				for (j = 0; j < nb; ++j)
				{
					ia[rowstart+j] = nzperrow;
				}
			}
		}
	}
	else {
		for (i = 0; i < ROW; ++i)
		{
			rowstart = i*nb + 1;
			colblock = IA[i+1] - IA[i];
			nzperrow = colblock*nb;
			for (j = 0; j < nb; ++j)
			{
				ia[rowstart+j] = nzperrow;
			}
		}
	}
	
	//-----------------------------------------------------
	// Generate the real 'ia' for CSR of A
	//-----------------------------------------------------
	
	ia[0] = 0;
	for (i = 1; i <= rowA; ++i)
	{
		ia[i] += ia[i-1];
	}
	
	//-----------------------------------------------------
	// Generate 'ja' and 'a' for CSR of A
	//-----------------------------------------------------
	
	switch (storage_manner)
	{
		case 0: // each non-zero block elements are stored in row-major order
		{
			if (ROW > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i, k, j, rowstart, colstart, vp, mr, ap, jap, mc) ////num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1) myend = mybegin+stride_i;
					else myend = ROW;
					for (i=mybegin; i<myend; ++i)
					{
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							j = JA[k];
							rowstart = i*nb;
							colstart = j*nb;
							vp = &val[k*jump];
							for (mr = 0; mr < nb; mr ++)
							{
								ap  = &a[ia[rowstart]];
								jap = &ja[ia[rowstart]];
								for (mc = 0; mc < nb; mc ++)
								{
									*ap = *vp;
									*jap = colstart + mc;
									vp ++; ap ++; jap ++;
								}
								ia[rowstart] += nb;
								rowstart ++;
							}
						}
					}
				}
			}
			else {
				for (i = 0; i < ROW; ++i)
				{
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						j = JA[k];
						rowstart = i*nb;
						colstart = j*nb;
						vp = &val[k*jump];
						for (mr = 0; mr < nb; mr ++)
						{
							ap  = &a[ia[rowstart]];
							jap = &ja[ia[rowstart]];
							for (mc = 0; mc < nb; mc ++)
							{
								*ap = *vp;
								*jap = colstart + mc;
								vp ++; ap ++; jap ++;
							}
							ia[rowstart] += nb;
							rowstart ++;
						}
					}
				}
			}
		}
			break;
            
		case 1: // each non-zero block elements are stored in column-major order
		{
			for (i = 0; i < ROW; ++i)
			{
				for (k = IA[i]; k < IA[i+1]; ++k)
				{
					j = JA[k];
					rowstart0 = i*nb;
					colstart0 = j*nb;
					vp = &val[k*jump];
					for (mc = 0; mc < nb; mc ++)
					{
						rowstart = rowstart0;
						colstart = colstart0 + mc;
						for (mr = 0; mr < nb; mr ++)
						{
							a[ia[rowstart]] = *vp; 
							ja[ia[rowstart]] = colstart; 
							vp ++; ia[rowstart]++; rowstart++;
						}
					}
				}
			}
		}
			break;
	}
	
	//-----------------------------------------------------
	// Map back the real 'ia' for CSR of A
	//-----------------------------------------------------
	
	for (i = rowA; i > 0; i --)
	{
		ia[i] = ia[i-1];
	}
	ia[0] = 0; 
    
	return (A);   
}
#endif // OMP

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
