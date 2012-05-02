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
 * \param alpha  REAL factor alpha
 * \param A      Pointer to the dBSR matrix
 * \param x      Pointer to the array x
 * \param beta   REAL factor beta
 * \param y      Pointer to the array y
 *
 * \author Zhou Zhiyang
 * \date   2010/10/25  
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
    if (temp != 1.0) {
        if (temp == 0.0) {
            memset(y, 0X0, size*sizeof(REAL));
        }
        else {
            for (i = size; i--; ) y[i] *= temp;  // modified by Xiaozhe, 03/11/2011
        }
    }
    
    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    
    for (i = 0; i < ROW; ++i) {
        py0 = &y[i*nb];
        ibegin = IA[i]; iend = IA[i+1];
        for (k = ibegin; k < iend; ++k) {
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
    
    if (alpha != 1.0) {
        for (i = size; i--; ++i) {
            y[i] *= alpha;
        }
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
 * \author Zhou Zhiyang
 * \date   2010/10/25  
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
    
    if (alpha == 0.0) {
        return; // Nothting to compute
    }
    
    //-------------------------------------------------
    //   y = (1.0/alpha)*y
    //-------------------------------------------------
    
    if (alpha != 1.0) {
        temp = 1.0 / alpha;
        for (i = size; i-- ; ) y[i] *= temp; // modified by Xiaozhe, 03/11/2011
    }
    
    //-----------------------------------------------------------------
    //   y += A*x (Core Computation)
    //   each non-zero block elements are stored in row-major order
    //-----------------------------------------------------------------
    
    switch (nb) {

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
    
    if (alpha != 1.0) {
        for (i = size; i--; ) {
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
 * \param A      Pointer to the dBSR matrix
 * \param x      Pointer to the array x
 * \param y      Pointer to the array y
 *
 * \author Zhou Zhiyang
 * \date   2010/10/25 
 *
 * \note Works for general nb (Xiaozhe) 
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
    
    switch (nb) {

    case 3: 
        {
            for (i = 0; i < ROW; ++i) {
                py0 = &y[i*3];
                ibegin = IA[i]; iend = IA[i+1];  // modified by Xiaozhe, 03/11/2011
                for (k = ibegin; k < iend; ++k) {
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
            for (i = 0; i < ROW; ++i) {
                py0 = &y[i*5];
                ibegin = IA[i]; iend = IA[i+1];  // modified by Xiaozhe, 03/11/2011
                for (k = ibegin; k < iend; ++k) {
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
            for (i = 0; i < ROW; ++i) {
                py0 = &y[i*7];
                ibegin = IA[i]; iend = IA[i+1];  // modified by Xiaozhe, 03/11/2011
                for (k = ibegin; k < iend; ++k) {
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
            for (i = 0; i < ROW; ++i) {
                py0 = &y[i*nb];
                ibegin = IA[i]; iend = IA[i+1]; // modified by Xiaozhe, 03/11/2011
                for (k = ibegin; k < iend; ++k) {
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
 * \author Chunsheng FENG ,Xiaoqiang YUE, Xiaozhe Hu
 * \date   08/08/2011
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *            Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void fasp_blas_dbsr_rap (dBSRmat *R, 
                         dBSRmat *A, 
                         dBSRmat *P, 
                         dBSRmat *B)
{
    const int row=R->ROW, col=P->COL,nb=A->nb, nb2=A->nb*A->nb;
    unsigned int nB=A->NNZ;
    int i,i1,j,jj,k,length;    
    int begin_row,end_row,begin_rowA,end_rowA,begin_rowR,end_rowR;
    int istart,iistart,count;
    
    double *rj=R->val, *aj=A->val, *pj=P->val, *acj;
    int *ir=R->IA, *ia=A->IA, *ip=P->IA, *iac;
    int *jr=R->JA, *ja=A->JA, *jp=P->JA, *jac;
    
    int *index=(int *)fasp_mem_calloc(A->COL,sizeof(int));
    
    double *smat_tmp=(double *)fasp_mem_calloc(nb2,sizeof(double));
    
    int *iindex=(int *)fasp_mem_calloc(col,sizeof(int));    
    
    for (i=0; i<A->COL; ++i) index[i] = -2;
    
    memcpy(iindex,index,col*sizeof(int));
    
    jac=(int*)fasp_mem_calloc(nB,sizeof(int));    
    
    iac=(int*)fasp_mem_calloc(row+1,sizeof(int));    
    
    double *temp=(double*)fasp_mem_calloc(A->COL*nb2,sizeof(double));
    
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
            jac=(int*)fasp_mem_realloc(jac, nB*sizeof(int));
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
    
    jac=(int*)fasp_mem_realloc(jac,(iac[row])*sizeof(int));
    
    acj=(double*)fasp_mem_calloc(iac[row]*nb2,sizeof(double));
    
    int *BTindex=(int*)fasp_mem_calloc(col,sizeof(int));
    
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
