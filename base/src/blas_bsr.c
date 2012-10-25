/*! \file blas_bsr.c
 *  \brief BLAS operations for sparse matrices in BSR format.
 */

#include <math.h>
#include <omp.h>

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
 * \param alpha  REAL factor alpha
 * \param A      Pointer to the dBSR matrix
 * \param x      Pointer to the array x
 * \param beta   REAL factor beta
 * \param y      Pointer to the array y
 *
 * \author Zhiyang Zhou
 * \date   10/25/2010  
 *
 * Modified by Chunsheng Feng, Zheng Li
 * \date   06/29/2012
 *
 * \note Works for general nb (Xiaozhe)
 */
void fasp_blas_dbsr_aAxpby (const REAL alpha, 
                                  dBSRmat *A, 
                                  REAL *x, 
                            const REAL beta, 
                                  REAL *y )
{
    /* members of A */
    INT  ROW  = A->ROW;
    INT  nb   = A->nb;
    INT  *IA  = A->IA;
    INT  *JA  = A->JA;
    REAL *val = A->val;

    /* local variables */
    INT     size = ROW*nb;
    INT     jump = nb*nb;
    INT     i,j,k,iend;
    REAL    temp;
    REAL   *pA  = NULL;
    REAL   *px0 = NULL;
    REAL   *py0 = NULL;
    REAL   *py  = NULL;

    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP 
    if ( ROW > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif

    //----------------------------------------------
    //   Treat (alpha == 0.0) computation 
    //----------------------------------------------

    if (alpha == 0.0) {
        fasp_blas_array_ax(size, beta, y);
        return;
    }

    //-------------------------------------------------
    //   y = (beta/alpha)*y
    //-------------------------------------------------

    temp = beta / alpha;
    if (temp != 1.0) {
        if (temp == 0.0) {
            memset(y, 0X0, size*sizeof(REAL));
        }
        else {
	    //for (i = size; i--; ) y[i] *= temp;  // modified by Xiaozhe, 03/11/2011
            fasp_blas_array_ax(size, temp, y);
        }
    }

    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------

    switch (nb)
    {
        case 2: 
        {
            if (use_openmp) {
                INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend)
#endif
                for (myid =0; myid < nthreads; myid++) {
                     FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                     for (i=mybegin; i < myend; ++i) {
                          py0 = &y[i*2];
                          iend = IA[i+1];
                          for (k = IA[i]; k < iend; ++k) {
                               j = JA[k];
                               pA = val+k*4; // &val[k*jump];
                               px0 = x+j*2; // &x[j*nb];
                               py = py0;
                               fasp_blas_smat_ypAx_nc2( pA, px0, py );
                          }
                     }
                }
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    py0 = &y[i*2];
                    iend = IA[i+1];
                    for (k = IA[i]; k < iend; ++k) {
                        j = JA[k];
                        pA = val+k*4; // &val[k*jump];
                        px0 = x+j*2; // &x[j*nb];
                        py = py0;
                        fasp_blas_smat_ypAx_nc2( pA, px0, py );
                    }
                }
            }
        }
        break;
        case 3: 
        {
            if (use_openmp) {
                INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend ) 
#endif
                for (myid =0; myid < nthreads; myid++) {
                     FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                     for (i=mybegin; i < myend; ++i) {
                          py0 = &y[i*3];
                          iend = IA[i+1];
                          for (k = IA[i]; k < iend; ++k) {
                               j = JA[k];
                               pA = val+k*9; // &val[k*jump];
                               px0 = x+j*3; // &x[j*nb];
                               py = py0;
                               fasp_blas_smat_ypAx_nc3( pA, px0, py );
                          }
                     }
                }
            }
            else {
                 for (i = 0; i < ROW; ++i) {
                      py0 = &y[i*3];
                      iend = IA[i+1];
                      for (k = IA[i]; k < iend; ++k) {
                           j = JA[k];
                           pA = val+k*9; // &val[k*jump];
                           px0 = x+j*3; // &x[j*nb];
                           py = py0;
                           fasp_blas_smat_ypAx_nc3( pA, px0, py );
                      }
                 }
            }
        }
        break;

        case 5:
        {
            if (use_openmp) {
                INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend )
#endif
                for (myid =0; myid < nthreads; myid++) {
                     FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                     for (i=mybegin; i < myend; ++i) {
                          py0 = &y[i*5];
                          iend = IA[i+1];
                          for (k = IA[i]; k < iend; ++k) {
                               j = JA[k];
                               pA = val+k*25; // &val[k*jump];
                               px0 = x+j*5; // &x[j*nb];
                               py = py0;
                               fasp_blas_smat_ypAx_nc5( pA, px0, py );
                          }
                     }
                }
            }
            else {
                 for (i = 0; i < ROW; ++i) {
                      py0 = &y[i*5];
                      iend = IA[i+1];
                      for (k = IA[i]; k < iend; ++k) {
                           j = JA[k];
                           pA = val+k*25; // &val[k*jump];
                           px0 = x+j*5; // &x[j*nb];
                           py = py0;
                           fasp_blas_smat_ypAx_nc5( pA, px0, py );
                      }
                 }
            }
        }
        break;
        case 7:
        {
            if (use_openmp) {
                INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend )
#endif
                for (myid =0; myid < nthreads; myid++) {
                    FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                    for (i=mybegin; i < myend; ++i) {
                        py0 = &y[i*7];
                        iend = IA[i+1];
                        for (k = IA[i]; k < iend; ++k) {
                            j = JA[k];
                            pA = val+k*49; // &val[k*jump];
                            px0 = x+j*7; // &x[j*nb]；
                            py = py0;
                            fasp_blas_smat_ypAx_nc7( pA, px0, py );
                        }
                    }
                } 
            }
            else {
                 for (i = 0; i < ROW; ++i) {
                      py0 = &y[i*7];
                      iend = IA[i+1];
                      for (k = IA[i]; k < iend; ++k) {
                           j = JA[k];
                           pA = val+k*49; // &val[k*jump];
                           px0 = x+j*7; // &x[j*nb];
                           py = py0;
                           fasp_blas_smat_ypAx_nc7( pA, px0, py );
                      }
                 }
            }
        }
        break;
        
        default: 
        {
             if (use_openmp) {
                 INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend)
#endif
                 for (myid =0; myid < nthreads; myid++) {
                      FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                      for (i=mybegin; i < myend; ++i) {
                           py0 = &y[i*nb];
                           iend = IA[i+1];
                           for (k = IA[i]; k < iend; ++k) {
                                j = JA[k];
                                pA = val+k*jump; // &val[k*jump];
                                px0 = x+j*nb; // &x[j*nb];
                                py = py0;
                                fasp_blas_smat_ypAx( pA, px0, py, nb );
                           }
                      }
                 }
             }
             else {
                  for (i = 0; i < ROW; ++i) {
                      py0 = &y[i*nb];
                      iend = IA[i+1];
                      for (k = IA[i]; k < iend; ++k) {
                           j = JA[k];
                           pA = val+k*jump; // &val[k*jump];
                           px0 = x+j*nb; // &x[j*nb];
                           py = py0;
                           fasp_blas_smat_ypAx( pA, px0, py, nb );
                      }
                  }
             }
        }
        break;
    }

    //------------------------------------------
    //   y = alpha*y
    //------------------------------------------

    if (alpha != 1.0) {
        fasp_blas_array_ax(size, alpha, y);
    }   
}


/*!
 * \fn void fasp_blas_dbsr_aAxpy ( const REAL alpha, dBSRmat *A, REAL *x, REAL *y )
 *
 * \brief Compute y := alpha*A*x + y
 *
 * \param alpha  REAL factor alpha
 * \param A      Pointer to the dBSR matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 * \author Zhiyang Zhou
 * \date   10/25/2010  
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue  
 * \date   05/23/2012    
 *
 * \note Works for general nb (Xiaozhe)
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
	REAL *val = A->val;

	/* local variables */
	INT     size = ROW*nb;
	INT     jump = nb*nb;
	INT     i,j,k, iend;
	REAL  temp = 0.0;
	REAL *pA   = NULL;
	REAL *px0  = NULL;
	REAL *py0  = NULL;
	REAL *py   = NULL;

	INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP 
	if ( ROW > OPENMP_HOLDS ) {
            use_openmp = TRUE;
            nthreads = FASP_GET_NUM_THREADS();
	}
#endif

	//----------------------------------------------
	//   Treat (alpha == 0.0) computation 
	//----------------------------------------------

	if (alpha == 0.0){
		return; // Nothting to compute
	}

	//-------------------------------------------------
	//   y = (1.0/alpha)*y
	//-------------------------------------------------

	if (alpha != 1.0){
		temp = 1.0 / alpha;
		fasp_blas_array_ax(size, temp, y);
	}

	//-----------------------------------------------------------------
	//   y += A*x (Core Computation)
	//   each non-zero block elements are stored in row-major order
	//-----------------------------------------------------------------

	switch (nb)
	{
		case 2: 
			{
				if (use_openmp) {
					INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend)
#endif
					for (myid =0; myid < nthreads; myid++) {
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i) {
							py0 = &y[i*2];
							iend = IA[i+1];
							for (k = IA[i]; k < iend; ++k) {
								j = JA[k];
								pA = val+k*4; // &val[k*jump];
								px0 = x+j*2; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc2( pA, px0, py );
							}
						}
					}
				}
				else {
					for (i = 0; i < ROW; ++i) {
						py0 = &y[i*2];
						iend = IA[i+1];
						for (k = IA[i]; k < iend; ++k) {
							j = JA[k];
							pA = val+k*4; // &val[k*jump];
							px0 = x+j*2; // &x[j*nb];
							py = py0;
							fasp_blas_smat_ypAx_nc2( pA, px0, py );
						}
					}
				}
			}
			break;

		case 3: 
			{
				if (use_openmp) {
					INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend ) 
#endif
					for (myid =0; myid < nthreads; myid++) {
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i) {
							py0 = &y[i*3];
							iend = IA[i+1];
							for (k = IA[i]; k < iend; ++k) {
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );
							}    
						}
					}
				}
				else {
					for (i = 0; i < ROW; ++i){
						py0 = &y[i*3];
						iend = IA[i+1];
						for (k = IA[i]; k < iend; ++k) {
							j = JA[k];
							pA = val+k*9; // &val[k*jump];
							px0 = x+j*3; // &x[j*nb];
							py = py0;
							fasp_blas_smat_ypAx_nc3( pA, px0, py );
						}    
					}
				}
			}
			break;

		case 5:
			{
				if (use_openmp) {
					INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for  private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend )
#endif
					for (myid =0; myid < nthreads; myid++) {
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i) {
							py0 = &y[i*5];
							iend = IA[i+1];
							for (k = IA[i]; k < iend; ++k) {
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );
							}  
						}
					}
				}
				else {
					for (i = 0; i < ROW; ++i){
						py0 = &y[i*5];
						iend = IA[i+1];
						for (k = IA[i]; k < iend; ++k) {
							j = JA[k];
							pA = val+k*25; // &val[k*jump];
							px0 = x+j*5; // &x[j*nb];
							py = py0;
							fasp_blas_smat_ypAx_nc5( pA, px0, py );
						}
					}
				}
			}
			break;
		case 7:
			{
				if (use_openmp) {
					INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend )
#endif
					for (myid =0; myid < nthreads; myid++) {
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i) {
							py0 = &y[i*7];
							iend = IA[i+1];
							for (k = IA[i]; k < iend; ++k) { 
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py ); 
							}  
						}
					}
				}
				else {
					for (i = 0; i < ROW; ++i) {
						py0 = &y[i*7];
						iend = IA[i+1];
						for (k = IA[i]; k < iend; ++k) {
							j = JA[k];
							pA = val+k*49; // &val[k*jump];
							px0 = x+j*7; // &x[j*nb];
							py = py0;
							fasp_blas_smat_ypAx_nc7( pA, px0, py );
						}

					}
				}
			}
			break;

		default: 
			{
				if (use_openmp) {
					INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, py0, k, j, pA, px0, py,iend)
#endif
					for (myid =0; myid < nthreads; myid++) {
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i) {
							py0 = &y[i*nb];
							iend = IA[i+1];
							for (k = IA[i]; k < iend; ++k) {
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );
							}  

						}
					}
				}
				else {
					for (i = 0; i < ROW; ++i) {
						py0 = &y[i*nb];
						iend = IA[i+1];
						for (k = IA[i]; k < iend; ++k) {
							j = JA[k];
							pA = val+k*jump; // &val[k*jump];
							px0 = x+j*nb; // &x[j*nb];
							py = py0;
							fasp_blas_smat_ypAx( pA, px0, py, nb );
						}

					}
				}
			}
			break;
	}

	//------------------------------------------
	//   y = alpha*y
	//------------------------------------------

	if (alpha != 1.0){
		fasp_blas_array_ax(size, alpha, y);
	} 
	return;
}

/*!
 * \fn void fasp_blas_dbsr_mxv ( dBSRmat *A, REAL *x, REAL *y )
 *
 * \brief Compute y := A*x
 *
 * \param A      Pointer to the dBSR matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 * \author Zhiyang Zhou
 * \date   10/25/2010 
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue 
 * \date   05/23/2012    
 *
 * \note Works for general nb (Xiaozhe) 
 */
void fasp_blas_dbsr_mxv(dBSRmat *A, 
		                REAL *x, 
		                REAL *y) 
{
	/* members of A */
	INT     ROW = A->ROW;
	INT     nb  = A->nb;
	INT    *IA  = A->IA;
	INT    *JA  = A->JA;
	REAL *val = A->val;

	/* local variables */
	INT     size = ROW*nb;
	INT     jump = nb*nb;
	INT     i,j,k, num_nnz_row;
	REAL *pA  = NULL;
	REAL *px0 = NULL;
	REAL *py0 = NULL;
	REAL *py  = NULL;

	INT use_openmp = FALSE;

#ifdef _OPENMP 
        INT myid, mybegin, myend, nthreads;
	if ( ROW > OPENMP_HOLDS ) {
            use_openmp = TRUE;
            nthreads = FASP_GET_NUM_THREADS();
	}
#endif

	//-----------------------------------------------------------------
	//  zero out 'y' 
	//-----------------------------------------------------------------
	fasp_array_set(size, y, 0.0);

	//-----------------------------------------------------------------
	//   y = A*x (Core Computation)
	//   each non-zero block elements are stored in row-major order
	//-----------------------------------------------------------------

	switch (nb)
	{
		case 3: 
			{
				if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py)
					{
						myid = omp_get_thread_num();
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i)
						{
							py0 = &y[i*3];
							num_nnz_row = IA[i+1] - IA[i];
							switch(num_nnz_row)
							{
								case 3:
									k = IA[i];
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									break;
								case 4:
									k = IA[i];
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									break;
								case 5:
									k = IA[i];
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									break;
								case 6:
									k = IA[i];
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									break;
								case 7:
									k = IA[i];
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*9;
									px0 = x+j*3;
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );

									break;
								default:
									for (k = IA[i]; k < IA[i+1]; ++k)
									{
										j = JA[k];
										pA = val+k*9;
										px0 = x+j*3;
										py = py0;
										fasp_blas_smat_ypAx_nc3( pA, px0, py );
									}
									break;
							}
						}
					}
#endif
				}
				else {
					for (i = 0; i < ROW; ++i)
					{
						py0 = &y[i*3];
						num_nnz_row = IA[i+1] - IA[i];
						switch(num_nnz_row)
						{
							case 3:
								k = IA[i];
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								break;
							case 4:
								k = IA[i];
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								break;
							case 5:
								k = IA[i];
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								break;
							case 6:
								k = IA[i];
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								break;
							case 7:
								k = IA[i];
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*9; // &val[k*jump];
								px0 = x+j*3; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc3( pA, px0, py );

								break;
							default:
								for (k = IA[i]; k < IA[i+1]; ++k)
								{
									j = JA[k];
									pA = val+k*9; // &val[k*jump];
									px0 = x+j*3; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc3( pA, px0, py );
								}
								break;
						}
					}
				}
			}
			break;

		case 5:
			{
				if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py)
					{
						myid = omp_get_thread_num();
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i)
						{
							py0 = &y[i*5];
							num_nnz_row = IA[i+1] - IA[i];
							switch(num_nnz_row)
							{
								case 3:
									k = IA[i];
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									break;
								case 4:
									k = IA[i];
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									break;
								case 5:
									k = IA[i];
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									break;
								case 6:
									k = IA[i];
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									break;
								case 7:
									k = IA[i];
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );

									break;
								default:
									for (k = IA[i]; k < IA[i+1]; ++k)
									{
										j = JA[k];
										pA = val+k*25; // &val[k*jump];
										px0 = x+j*5; // &x[j*nb];
										py = py0;
										fasp_blas_smat_ypAx_nc5( pA, px0, py );
									}
									break;
							}
						}
					}
#endif
				}
				else {
					for (i = 0; i < ROW; ++i)
					{
						py0 = &y[i*5];
						num_nnz_row = IA[i+1] - IA[i];
						switch(num_nnz_row)
						{
							case 3:
								k = IA[i];
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								break;
							case 4:
								k = IA[i];
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								break;
							case 5:
								k = IA[i];
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								break;
							case 6:
								k = IA[i];
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								break;
							case 7:
								k = IA[i];
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*25; // &val[k*jump];
								px0 = x+j*5; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc5( pA, px0, py );

								break;
							default:
								for (k = IA[i]; k < IA[i+1]; ++k)
								{
									j = JA[k];
									pA = val+k*25; // &val[k*jump];
									px0 = x+j*5; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc5( pA, px0, py );
								}
								break;
						}
					}
				}
			}
			break;

		case 7:
			{
				if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py)
					{
						myid = omp_get_thread_num();
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i)
						{
							py0 = &y[i*7];
							num_nnz_row = IA[i+1] - IA[i];
							switch(num_nnz_row)
							{
								case 3:
									k = IA[i];
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									break;
								case 4:
									k = IA[i];
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									break;
								case 5:
									k = IA[i];
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									break;
								case 6:
									k = IA[i];
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									break;
								case 7:
									k = IA[i];
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									k ++;
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );

									break;
								default:
									for (k = IA[i]; k < IA[i+1]; ++k)
									{
										j = JA[k];
										pA = val+k*49; // &val[k*jump];
										px0 = x+j*7; // &x[j*nb];
										py = py0;
										fasp_blas_smat_ypAx_nc7( pA, px0, py );
									}
									break;
							}
						}
					}
#endif
				}
				else {
					for (i = 0; i < ROW; ++i)
					{
						py0 = &y[i*7];
						num_nnz_row = IA[i+1] - IA[i];
						switch(num_nnz_row)
						{
							case 3:
								k = IA[i];
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								break;
							case 4:
								k = IA[i];
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								break;
							case 5:
								k = IA[i];
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								break;
							case 6:
								k = IA[i];
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								break;
							case 7:
								k = IA[i];
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								k ++;
								j = JA[k];
								pA = val+k*49; // &val[k*jump];
								px0 = x+j*7; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx_nc7( pA, px0, py );

								break;
							default:
								for (k = IA[i]; k < IA[i+1]; ++k)
								{
									j = JA[k];
									pA = val+k*49; // &val[k*jump];
									px0 = x+j*7; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx_nc7( pA, px0, py );
								}
								break;
						}
					}
				}
			}
			break;

		default: 
			{
				if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel private(myid, mybegin, myend, i, py0, num_nnz_row, k, j, pA, px0, py)
					{
						myid = omp_get_thread_num();
						FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
						for (i=mybegin; i < myend; ++i)
						{
							py0 = &y[i*nb];
							num_nnz_row = IA[i+1] - IA[i];
							switch(num_nnz_row)
							{
								case 3:
									k = IA[i];
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									break;
								case 4:
									k = IA[i];
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									break;
								case 5:
									k = IA[i];
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									break;
								case 6:
									k = IA[i];
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									break;
								case 7:
									k = IA[i];
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									k ++;
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );

									break;
								default:
									for (k = IA[i]; k < IA[i+1]; ++k)
									{
										j = JA[k];
										pA = val+k*jump; // &val[k*jump];
										px0 = x+j*nb; // &x[j*nb];
										py = py0;
										fasp_blas_smat_ypAx( pA, px0, py, nb );
									}
									break;
							}
						}  
					}
#endif
				}
				else {
					for (i = 0; i < ROW; ++i)
					{
						py0 = &y[i*nb];
						num_nnz_row = IA[i+1] - IA[i];
						switch(num_nnz_row)
						{
							case 3:
								k = IA[i];
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								break;
							case 4:
								k = IA[i];
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								break;
							case 5:
								k = IA[i];
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								break;
							case 6:
								k = IA[i];
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								break;
							case 7:
								k = IA[i];
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								k ++;
								j = JA[k];
								pA = val+k*jump; // &val[k*jump];
								px0 = x+j*nb; // &x[j*nb];
								py = py0;
								fasp_blas_smat_ypAx( pA, px0, py, nb );

								break;
							default:
								for (k = IA[i]; k < IA[i+1]; ++k)
								{
									j = JA[k];
									pA = val+k*jump; // &val[k*jump];
									px0 = x+j*nb; // &x[j*nb];
									py = py0;
									fasp_blas_smat_ypAx( pA, px0, py, nb );
								}
								break;
						}
					}
				}
			}
			break;
	}
}


/**
 * \fn void fasp_blas_dbsr_rap1 (dBSRmat *R, dBSRmat *A, dBSRmat *P, dBSRmat *B)
 *
 * \brief dBSRmat sparse matrix multiplication B=R*A*P
 *
 * \param R   Pointer to the dBSRmat matrix
 * \param A   Pointer to the dBSRmat matrix
 * \param P   Pointer to the dBSRmat matrix
 * \param B   Pointer to dBSRmat matrix equal to R*A*P (output)
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Xiaozhe Hu
 * \date   08/08/2011
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *            Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void fasp_blas_dbsr_rap1 (dBSRmat *R,
                          dBSRmat *A, 
                          dBSRmat *P, 
                          dBSRmat *B)
{
    const INT row=R->ROW, col=P->COL,nb=A->nb, nb2=A->nb*A->nb;
    INT nB=A->NNZ;
    INT i,i1,j,jj,k,length;    
    INT begin_row,end_row,begin_rowA,end_rowA,begin_rowR,end_rowR;
    INT istart,iistart,count;

    REAL *rj=R->val, *aj=A->val, *pj=P->val, *acj;
    INT *ir=R->IA, *ia=A->IA, *ip=P->IA, *iac;
    INT *jr=R->JA, *ja=A->JA, *jp=P->JA, *jac;

    INT *index=(INT *)fasp_mem_calloc(A->COL,sizeof(INT));

    REAL *smat_tmp=(REAL *)fasp_mem_calloc(nb2,sizeof(REAL));

    INT *iindex=(INT *)fasp_mem_calloc(col,sizeof(INT));    

    for (i=0; i<A->COL; ++i) index[i] = -2;

    memcpy(iindex,index,col*sizeof(INT));

    jac=(INT*)fasp_mem_calloc(nB,sizeof(INT));    

    iac=(INT*)fasp_mem_calloc(row+1,sizeof(INT));    

    REAL *temp=(REAL*)fasp_mem_calloc(A->COL*nb2,sizeof(REAL));

    iac[0] = 0;

    // First loop: form sparsity partern of R*A*P
    for (i=0; i < row; ++i) {
        // reset istart and length at the begining of each loop
        istart = -1; length = 0; i1 = i+1;

       // go across the rows in R
       begin_rowR=ir[i]; end_rowR=ir[i1];
       for (jj=begin_rowR; jj<end_rowR; ++jj) {
           j = jr[N2C(jj)];
           // for each column in A
           begin_rowA=ia[j]; end_rowA=ia[j+1];    
           for (k=begin_rowA; k<end_rowA; ++k) {
               if (index[N2C(ja[N2C(k)])] == -2) {
                   index[N2C(ja[N2C(k)])] = istart;
                   istart = ja[N2C(k)];
                   ++length;
               }
           }
       }    

       // book-keeping [reseting length and setting iistart]
       count = length; iistart = -1; length = 0;

       // use each column that would have resulted from R*A
       for (j=0; j < count; ++j) {
           jj = istart;
           istart = index[istart];
           index[N2C(jj)] = -2;

           // go across the row of P
           begin_row=ip[jj]; end_row=ip[jj+1];
           for (k=begin_row; k<end_row; ++k) {
               // pull out the appropriate columns of P
               if (iindex[N2C(jp[N2C(k)])] == -2){
                   iindex[N2C(jp[N2C(k)])] = iistart;
                   iistart = jp[N2C(k)];
                   ++length;
               }
           } // end for k
       } // end for j

       // set B->IA
       iac[i1]=iac[i]+length;

       if (iac[i1]>nB) {
           nB=nB*2;
           jac=(INT*)fasp_mem_realloc(jac, nB*sizeof(INT));
       }

       // put the correct columns of p into the column list of the products
       begin_row=iac[i]; end_row=iac[i1];
       for (j=begin_row; j<end_row; ++j) {
           // put the value in B->JA
           jac[N2C(j)] = iistart;
           // set istart to the next value
           iistart = iindex[N2C(iistart)];
           // set the iindex spot to 0
           iindex[N2C(jac[j])] = -2;
       } // end j
    } // end i: First loop

    jac=(INT*)fasp_mem_realloc(jac,(iac[row])*sizeof(INT));

    acj=(REAL*)fasp_mem_calloc(iac[row]*nb2,sizeof(REAL));

    INT *BTindex=(INT*)fasp_mem_calloc(col,sizeof(INT));

    // Second loop: compute entries of R*A*P
    for (i=0; i<row; ++i) {
        i1 = i+1;
        // each col of B
        begin_row=iac[i]; end_row=iac[i1];
        for (j=begin_row; j<end_row; ++j) {
            BTindex[N2C(jac[N2C(j)])]=j;
        }
        // reset istart and length at the begining of each loop
        istart = -1; length = 0;

        // go across the rows in R
        begin_rowR=ir[i]; end_rowR=ir[i1];    
        for ( jj=begin_rowR; jj<end_rowR; ++jj ) {
            j = jr[N2C(jj)];
            // for each column in A
            begin_rowA=ia[j]; end_rowA=ia[j+1];
            for (k=begin_rowA; k<end_rowA; ++k) {
                if (index[N2C(ja[N2C(k)])] == -2) {
                    index[N2C(ja[N2C(k)])] = istart;
                    istart = ja[N2C(k)];
                    ++length;
                }
                fasp_blas_smat_mul(&rj[N2C(jj)*nb2],&aj[N2C(k)*nb2],smat_tmp,nb);
                //fasp_array_xpy(nb2,&temp[N2C(ja[N2C(k)])*nb2], smat_tmp );
                fasp_blas_array_axpy (nb2, 1.0, smat_tmp, &temp[N2C(ja[N2C(k)])*nb2]);

                //temp[N2C(ja[N2C(k)])]+=rj[N2C(jj)]*aj[N2C(k)];
                // change to   X = X+Y*Z
            }
        } 

        // book-keeping [reseting length and setting iistart]
        // use each column that would have resulted from R*A
        for (j=0; j<length; ++j) {
             jj = N2C(istart);
             istart = index[N2C(istart)];
             index[N2C(jj)] = -2;

             // go across the row of P
             begin_row=ip[jj]; end_row=ip[jj+1];    
             for (k=begin_row; k<end_row; ++k) {
                 // pull out the appropriate columns of P
                 //acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
                 fasp_blas_smat_mul(&temp[jj*nb2],&pj[k*nb2],smat_tmp,nb);
                 //fasp_array_xpy(nb2,&acj[BTindex[N2C(jp[N2C(k)])]*nb2], smat_tmp );
                 fasp_blas_array_axpy (nb2, 1.0, smat_tmp, &acj[BTindex[N2C(jp[N2C(k)])]*nb2]);

                 // change to   X = X+Y*Z
             }
             //temp[jj]=0.0; // change to   X[nb,nb] = 0;
             fasp_array_set(nb2,&temp[jj*nb2],0.0);
         }
     } // end for i: Second loop
     // setup coarse matrix B
     B->ROW=row; B->COL=col;
     B->IA=iac; B->JA=jac; B->val=acj;
     B->NNZ=B->IA[B->ROW]-B->IA[0];

     B->nb=A->nb;
     B->storage_manner = A->storage_manner;
     fasp_mem_free(temp);
     fasp_mem_free(index);
     fasp_mem_free(iindex);
     fasp_mem_free(BTindex);
     fasp_mem_free(smat_tmp);
}

/**
 * \fn void fasp_blas_dbsr_rap (dBSRmat *R, dBSRmat *A, dBSRmat *P, dBSRmat *B)
 *
 * \brief dBSRmat sparse matrix multiplication B=R*A*P
 *
 * \param R   Pointer to the dBSRmat matrix
 * \param A   Pointer to the dBSRmat matrix
 * \param P   Pointer to the dBSRmat matrix
 * \param B   Pointer to dBSRmat matrix equal to R*A*P (output)
 *
 * \author Chunsheng Feng, Zheng Li
 * \date   10/24/2012
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *            Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void fasp_blas_dbsr_rap (dBSRmat *R,
                         dBSRmat *A, 
                         dBSRmat *P, 
                         dBSRmat *B)
{
    const INT row=R->ROW, col=P->COL,nb=A->nb, nb2=A->nb*A->nb;
    INT nB=A->NNZ;

    REAL *rj=R->val, *aj=A->val, *pj=P->val, *acj;
    INT *ir=R->IA, *ia=A->IA, *ip=P->IA, *iac;
    INT *jr=R->JA, *ja=A->JA, *jp=P->JA, *jac;

    INT *Ps_marker = NULL;
    INT *As_marker = NULL;

#ifdef _OPENMP
    INT *P_marker = NULL;
    INT *A_marker = NULL;
    REAL *smat_tmp = NULL;
#endif

    INT i, i1, i2, i3, jj1, jj2, jj3;
    INT counter, jj_row_begining;

    INT nthreads = 1;

#ifdef _OPENMP
    INT myid, mybegin, myend, Ctemp;
    nthreads = FASP_GET_NUM_THREADS();
#endif

    INT n_coarse = row;
    INT n_fine   = A->ROW;
    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;

    Ps_marker = (INT *)fasp_mem_calloc(total_calloc, sizeof(INT));
    As_marker = Ps_marker + coarse_mul_nthreads;

    /*------------------------------------------------------*
     *  First Pass: Determine size of B and set up B_i  *
     *------------------------------------------------------*/
    iac = (INT *)fasp_mem_calloc(n_coarse+1, sizeof(INT));

    fasp_iarray_set(minus_one_length, Ps_marker, -1);

    REAL *tmp=(REAL *)fasp_mem_calloc(2*nthreads*nb2, sizeof(REAL));

#ifdef _OPENMP
    INT * RAP_temp = As_marker + fine_mul_nthreads;
    INT * part_end = RAP_temp + coarse_add_nthreads;

    if (n_coarse > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, Ctemp, P_marker, A_marker, counter, i, jj_row_begining, jj1, i1, jj2, i2, jj3, i3)
        for (myid = 0; myid < nthreads; myid++) {
            FASP_GET_START_END(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            counter  = 0;
            for (i = mybegin; i < myend; ++i) {
                P_marker[i] = counter;
                jj_row_begining = counter;
                counter ++;
                for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
                    i1 = jr[jj1];
                    for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                        i2 = ja[jj2];
                        if (A_marker[i2] != i) {
                            A_marker[i2] = i;
                            for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                                i3 = jp[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3] = counter;
                                    counter ++;
                                }
                            }
                        }
                    }
                }
                RAP_temp[i+myid] = jj_row_begining;
            }
            RAP_temp[myend+myid] = counter;
            part_end[myid] = myend + myid + 1;
        }
        fasp_iarray_cp(part_end[0], RAP_temp, iac);
        counter = part_end[0];
        Ctemp = 0;
        for (i1 = 1; i1 < nthreads; i1 ++) {
            Ctemp += RAP_temp[part_end[i1-1]-1];
            for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++) {
                iac[counter] = RAP_temp[jj1] + Ctemp;
                counter ++;
            }
        }
    }
    else {
#endif
        counter = 0;
        for (i = 0; i < row; ++ i) {
            Ps_marker[i] = counter;
            jj_row_begining = counter;
            counter ++;

            for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
                i1 = jr[jj1];
                for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                    i2 = ja[jj2];
                    if (As_marker[i2] != i) {
                        As_marker[i2] = i;
                        for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                            i3 = jp[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3] = counter;
                                counter ++;
                            }
                        }
                    }
                }
            }
            iac[i] = jj_row_begining;
        }
#ifdef _OPENMP
    }
#endif

    iac[row] = counter;

    jac=(INT*)fasp_mem_calloc(iac[row], sizeof(INT));

    acj=(REAL*)fasp_mem_calloc(iac[row]*nb2, sizeof(REAL));

    fasp_iarray_set(minus_one_length, Ps_marker, -1);

    /*------------------------------------------------------*
     *  Second Pass: compute entries of B=R*A*P             *
     *------------------------------------------------------*/
#ifdef _OPENMP
    if (n_coarse > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, Ctemp, P_marker, A_marker, counter, i, jj_row_begining, jj1, i1, jj2, i2, jj3, i3, smat_tmp)
        for (myid = 0; myid < nthreads; myid++) {
            FASP_GET_START_END(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            smat_tmp = tmp + myid*2*nb2; 
            counter = iac[mybegin];
            for (i = mybegin; i < myend; ++i) {
                P_marker[i] = counter;
                jj_row_begining = counter;
                jac[counter] = i;
                fasp_array_set(nb2, &acj[counter*nb2], 0x0);
                counter ++;

                for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
                    i1 = jr[jj1];
                    for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                        fasp_blas_smat_mul(&rj[jj1*nb2],&aj[jj2*nb2], smat_tmp, nb);
                        i2 = ja[jj2];
                        if (A_marker[i2] != i) {
                            A_marker[i2] = i;
                            for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                                i3 = jp[jj3];
                                fasp_blas_smat_mul(smat_tmp, &pj[jj3*nb2], smat_tmp+nb2, nb);
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3] = counter;
                                    fasp_array_cp(nb2, smat_tmp+nb2, &acj[counter*nb2]);
                                    jac[counter] = i3;
                                    counter ++;
                                }
                                else {
                                    fasp_blas_array_axpy (nb2, 1.0, smat_tmp+nb2, &acj[P_marker[i3]*nb2]);
                                }
                            }
                        }
                        else {
                            for (jj3 = ip[i2]; jj3 < ip[i2+1]; jj3 ++) {
                                i3 = jp[jj3];
                                fasp_blas_smat_mul(smat_tmp, &pj[jj3*nb2], smat_tmp+nb2, nb);
                                fasp_blas_array_axpy (nb2, 1.0, smat_tmp+nb2, &acj[P_marker[i3]*nb2]);
                            }
                        }
                    }
                }
            }
        }
    }
    else {
#endif
        counter = 0;
        for (i = 0; i < row; ++i) {
            Ps_marker[i] = counter;
            jj_row_begining = counter;
            jac[counter] = i;
            fasp_array_set(nb2, &acj[counter*nb2], 0x0);
            counter ++;

            for (jj1 = ir[i]; jj1 < ir[i+1]; ++jj1) {
                i1 = jr[jj1];
                for (jj2 = ia[i1]; jj2 < ia[i1+1]; ++jj2) {
                    fasp_blas_smat_mul(&rj[jj1*nb2],&aj[jj2*nb2], tmp, nb);
                    i2 = ja[jj2];
                    if (As_marker[i2] != i) {
                        As_marker[i2] = i;
                        for (jj3 = ip[i2]; jj3 < ip[i2+1]; ++jj3) {
                            i3 = jp[jj3];
                            fasp_blas_smat_mul(tmp, &pj[jj3*nb2], tmp+nb2, nb);
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3] = counter;
                                fasp_array_cp(nb2, tmp+nb2, &acj[counter*nb2]);
                                jac[counter] = i3;
                                counter ++;
                            }
                            else {
                                fasp_blas_array_axpy (nb2, 1.0, tmp+nb2, &acj[Ps_marker[i3]*nb2]);
                            }
                        }
                    }
                    else {
                        for (jj3 = ip[i2]; jj3 < ip[i2+1]; jj3 ++) {
                            i3 = jp[jj3];
                            fasp_blas_smat_mul(tmp, &pj[jj3*nb2], tmp+nb2, nb);
                            fasp_blas_array_axpy (nb2, 1.0, tmp+nb2, &acj[Ps_marker[i3]*nb2]);
                        }
                    }
                }
            }  
        }
#ifdef _OPENMP
    }
#endif
    // setup coarse matrix B
    B->ROW=row; B->COL=col;
    B->IA=iac; B->JA=jac; B->val=acj;
    B->NNZ=B->IA[B->ROW]-B->IA[0];
    B->nb=A->nb;
    B->storage_manner = A->storage_manner;
    fasp_mem_free(Ps_marker);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
