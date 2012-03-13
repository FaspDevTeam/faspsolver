/*! \file sparse_bsr_omp.c
 *  \brief Simple operations for BSR format
 *
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

#if FASP_USE_OPENMP
/**
 * \fn dBSRmat fasp_dbsr_diaginv3_omp(dBSRmat *A, double *diaginv, int nthreads, int openmp_holds)
 * \brief Compute B := D^{-1}*A, where 'D' is the block diagonal part of A.
 * \param A pointer to the dBSRmat matrix
 * \param diaginv pointer to the inverses of all the diagonal blocks
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date Nov/28/2011
 */
dBSRmat fasp_dbsr_diaginv3_omp (dBSRmat *A, 
                                double *diaginv, 
                                int nthreads, 
                                int openmp_holds)
{
	dBSRmat B;
    
	// members of A 
	int     ROW = A->ROW;
	int     ROW_plus_one = ROW+1;
	int     COL = A->COL;
	int     NNZ = A->NNZ;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	int    *IAb  = NULL;
	int    *JAb  = NULL;
	double *valb = NULL;
	
	int nb2  = nb*nb;
	int i,j,k,m;
	int myid;
	int mybegin;
	int stride_i;
	int myend;
	
	// Create a dBSRmat 'B'
	B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
	
	IAb  = B.IA;
	JAb  = B.JA;
	valb = B.val;
	
	fasp_iarray_cp_omp(ROW_plus_one, IA, IAb, nthreads,openmp_holds);
	fasp_iarray_cp_omp(NNZ, JA, JAb, nthreads,openmp_holds);
	
	switch (nb)
	{
		case 2:
			// main loop 
			if (ROW > openmp_holds) {
				stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1)  myend = mybegin+stride_i;
					else myend = ROW;
					for (i=mybegin; i < myend; ++i)
					{
						// get the diagonal sub-blocks
                        
                        k = IA[i];
                        m = k*4;
                        memcpy(diaginv+i*4, val+m, 4*sizeof(double));
                        fasp_smat_identity_nc2(valb+m);
                        
						// compute the inverses of the diagonal sub-blocks 
						fasp_blas_smat_inv_nc2(diaginv+i*4);
						// compute D^{-1}*A
						for (k = IA[i]+1; k < IA[i+1]; ++k)
						{
							m = k*4;
							j = JA[k];
                            fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
						}
					}// end of main loop
				}
			}
			else {
                // main loop 
                for (i = 0; i < ROW; ++i)
                {
                    // get the diagonal sub-blocks
                    k = IA[i];
                    m = k*4;
                    memcpy(diaginv+i*4, val+m, 4*sizeof(double));
                    fasp_smat_identity_nc2(valb+m);
					
                    // compute the inverses of the diagonal sub-blocks 
                    fasp_blas_smat_inv_nc2(diaginv+i*4);
                    // compute D^{-1}*A
                    for (k = IA[i]+1; k < IA[i+1]; ++k)
                    {
                        m = k*4;
                        j = JA[k];
                        fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
                    }
                }// end of main loop
			}
			
			break;
		case 3:
			// main loop 
			if (ROW > openmp_holds) {
				stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1)  myend = mybegin+stride_i;
					else myend = ROW;
					for (i=mybegin; i < myend; ++i)
					{
						// get the diagonal sub-blocks
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							if (JA[k] == i)
							{
								m = k*9;
								memcpy(diaginv+i*9, val+m, 9*sizeof(double));
								fasp_smat_identity_nc3(valb+m);
							}
						}
						// compute the inverses of the diagonal sub-blocks 
						fasp_blas_smat_inv_nc3(diaginv+i*9);
						// compute D^{-1}*A
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							m = k*9;
							j = JA[k];
							if (j != i) fasp_blas_smat_mul_nc3(diaginv+i*9, val+m, valb+m);
						}
					}// end of main loop
				}
			}
			else {
				for (i = 0; i < ROW; ++i)
				{
					// get the diagonal sub-blocks
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						if (JA[k] == i)
						{
							m = k*9;
							memcpy(diaginv+i*9, val+m, 9*sizeof(double));
							fasp_smat_identity_nc3(valb+m);
						}
					}
					
					// compute the inverses of the diagonal sub-blocks 
					fasp_blas_smat_inv_nc3(diaginv+i*9);
					
					// compute D^{-1}*A
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						m = k*9;
						j = JA[k];
						if (j != i) fasp_blas_smat_mul_nc3(diaginv+i*9, val+m, valb+m);
					}
				}// end of main loop
			}
			
			break;
			
		case 5: 
			// main loop 
			if (ROW > openmp_holds) {
				stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1)  myend = mybegin+stride_i;
					else myend = ROW;
					for (i=mybegin; i < myend; ++i)
					{
						// get the diagonal sub-blocks
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							if (JA[k] == i)
							{
								m = k*25;
								memcpy(diaginv+i*25, val+m, 25*sizeof(double));
								fasp_smat_identity_nc5(valb+m);
							}
						}
						
						// compute the inverses of the diagonal sub-blocks 
						fasp_blas_smat_inv_nc5(diaginv+i*25);
						
						// compute D^{-1}*A
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							m = k*25;
							j = JA[k];
							if (j != i) fasp_blas_smat_mul_nc5(diaginv+i*25, val+m, valb+m);
						}
					}// end of main loop
				}
			}
			else {
				for (i = 0; i < ROW; ++i)
				{
					// get the diagonal sub-blocks
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						if (JA[k] == i)
						{
							m = k*25;
							memcpy(diaginv+i*25, val+m, 25*sizeof(double));
							fasp_smat_identity_nc5(valb+m);
						}
					}
					
					// compute the inverses of the diagonal sub-blocks 
					fasp_blas_smat_inv_nc5(diaginv+i*25);
					
					// compute D^{-1}*A
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						m = k*25;
						j = JA[k];
						if (j != i) fasp_blas_smat_mul_nc5(diaginv+i*25, val+m, valb+m);
					}
				}// end of main loop
			}
			
			break;
			
		case 7:
			// main loop
			if (ROW > openmp_holds) {
				stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1)  myend = mybegin+stride_i;
					else myend = ROW;
					for (i=mybegin; i < myend; ++i)
					{
						// get the diagonal sub-blocks
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							if (JA[k] == i)
							{
								m = k*49;
								memcpy(diaginv+i*49, val+m, 49*sizeof(double));
								fasp_smat_identity_nc7(valb+m);
							}
						}
						
						// compute the inverses of the diagonal sub-blocks 
						fasp_blas_smat_inv_nc7(diaginv+i*49);
						
						// compute D^{-1}*A
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							m = k*49;
							j = JA[k];
							if (j != i) fasp_blas_smat_mul_nc7(diaginv+i*49, val+m, valb+m);
						}
					}// end of main loop
				}
			}
			else {
				for (i = 0; i < ROW; ++i)
				{
					// get the diagonal sub-blocks
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						if (JA[k] == i)
						{
							m = k*49;
							memcpy(diaginv+i*49, val+m, 49*sizeof(double));
							fasp_smat_identity_nc7(valb+m);
						}
					}
					
					// compute the inverses of the diagonal sub-blocks 
					fasp_blas_smat_inv_nc7(diaginv+i*49);
					
					// compute D^{-1}*A
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						m = k*49;
						j = JA[k];
						if (j != i) fasp_blas_smat_mul_nc7(diaginv+i*49, val+m, valb+m);
					}
				}// end of main loop
			}
			
			break;
			
		default:
			// main loop
			if (ROW > openmp_holds) {
				stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					mybegin = myid*stride_i;
					if(myid < nthreads-1)  myend = mybegin+stride_i;
					else myend = ROW;
					for (i=mybegin; i < myend; ++i)
					{
						// get the diagonal sub-blocks
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							if (JA[k] == i)
							{
								m = k*nb2;
								memcpy(diaginv+i*nb2, val+m, nb2*sizeof(double));
								fasp_smat_identity(valb+m, nb, nb2);
							}
						}
						
						// compute the inverses of the diagonal sub-blocks 
						fasp_blas_smat_inv(diaginv+i*nb2, nb);
						
						// compute D^{-1}*A
						for (k = IA[i]; k < IA[i+1]; ++k)
						{
							m = k*nb2;
							j = JA[k];
							if (j != i) fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
						}
					}// end of main loop
				}
			}
			else {
				for (i = 0; i < ROW; ++i)
				{
					// get the diagonal sub-blocks
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						if (JA[k] == i)
						{
							m = k*nb2;
							memcpy(diaginv+i*nb2, val+m, nb2*sizeof(double));
							fasp_smat_identity(valb+m, nb, nb2);
						}
					}
					
					// compute the inverses of the diagonal sub-blocks 
					fasp_blas_smat_inv(diaginv+i*nb2, nb);
					
					// compute D^{-1}*A
					for (k = IA[i]; k < IA[i+1]; ++k)
					{
						m = k*nb2;
						j = JA[k];
						if (j != i) fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
					}
				}// end of main loop
			}
			
			break;
	}
    
	return (B);
}
#endif // OMP

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
