/*! \file blas_bcsr.c
 *
 *  \brief BLAS operations for block_dCSRmat matrices
 */

#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_bdcsr_aAxpy (const REAL alpha, block_dCSRmat *A,
 *                                 REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  REAL factor a
 * \param A      Pointer to block_dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Xiaozhe Hu
 * \date   06/04/2010
 */
void fasp_blas_bdcsr_aAxpy (const REAL alpha,
                            block_dCSRmat *A,
                            REAL *x,
                            REAL *y)
{
    // information of A
    INT brow = A->brow;
    
    // local variables
    register dCSRmat *A11, *A12, *A21, *A22;
    register dCSRmat *A13, *A23, *A31, *A32, *A33;
    
    unsigned INT row1, col1;
    unsigned INT row2, col2;
    
    register REAL *x1, *x2, *y1, *y2;
    register REAL *x3, *y3;
    
    INT i,j;
    INT start_row;
    INT start_col;
    
    switch (brow) {
            
        case 2:
            A11 = A->blocks[0];
            A12 = A->blocks[1];
            A21 = A->blocks[2];
            A22 = A->blocks[3];
            
            row1 = A11->row;
            col1 = A11->col;
            
            x1 = x;
            x2 = &(x[col1]);
            y1 = y;
            y2 = &(y[row1]);
            
            // y1 = alpha*A11*x1 + alpha*A12*x2 + y1
            fasp_blas_dcsr_aAxpy(alpha, A11, x1, y1);
            fasp_blas_dcsr_aAxpy(alpha, A12, x2, y1);
            
            // y2 = alpha*A21*x1 + alpha*A22*x2 + y2
            fasp_blas_dcsr_aAxpy(alpha, A21, x1, y2);
            fasp_blas_dcsr_aAxpy(alpha, A22, x2, y2);
            
            break;
            
        case 4:
            A11 = A->blocks[0];
            A12 = A->blocks[1];
            A13 = A->blocks[2];
            A21 = A->blocks[3];
            A22 = A->blocks[4];
            A23 = A->blocks[5];
            A31 = A->blocks[6];
            A32 = A->blocks[7];
            A33 = A->blocks[8];
            
            row1 = A11->row;
            col1 = A11->col;
            row2 = A22->row;
            col2 = A22->col;
            
            x1 = x;
            x2 = &(x[col1]);
            x3 = &(x[col1+col2]);
            y1 = y;
            y2 = &(y[row1]);
            y3 = &(y[row1+row2]);
            
            // y1 = alpha*A11*x1 + alpha*A12*x2 + alpha*A13*x3 + y1
            fasp_blas_dcsr_aAxpy(alpha, A11, x1, y1);
            fasp_blas_dcsr_aAxpy(alpha, A12, x2, y1);
            fasp_blas_dcsr_aAxpy(alpha, A13, x3, y1);
            
            // y2 = alpha*A21*x1 + alpha*A22*x2 + alpha*A23*x3 + y2
            fasp_blas_dcsr_aAxpy(alpha, A21, x1, y2);
            fasp_blas_dcsr_aAxpy(alpha, A22, x2, y2);
            fasp_blas_dcsr_aAxpy(alpha, A23, x3, y2);
            
            // y3 = alpha*A31*x1 + alpha*A32*x2 + alpha*A33*x3 + y2
            fasp_blas_dcsr_aAxpy(alpha, A31, x1, y3);
            fasp_blas_dcsr_aAxpy(alpha, A32, x2, y3);
            fasp_blas_dcsr_aAxpy(alpha, A33, x3, y3);
            
            break;
            
        default:
            
            start_row = 0;
            start_col = 0;
            
            for (i=0; i<brow; i++) {
                
                for (j=0; j<brow; j++){
                    
                    if (A->blocks[i*brow+j]){
                        fasp_blas_dcsr_aAxpy(alpha, A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
                    }
                    start_col = start_col + A->blocks[j*brow+j]->col;
                }
                
                start_row = start_row + A->blocks[i*brow+i]->row;
                start_col = 0;
            }
            
            break;
            
    } // end of switch
    
}

/**
 * \fn void fasp_blas_bdcsr_mxv (block_dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A      Pointer to block_dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Chensong Zhang
 * \date   04/27/2013
 */
void fasp_blas_bdcsr_mxv (block_dCSRmat *A,
                          REAL *x,
                          REAL *y)
{
    // information of A
    INT brow = A->brow;
    
    // local variables
    register dCSRmat *A11, *A12, *A21, *A22;
    register dCSRmat *A13, *A23, *A31, *A32, *A33;
    
    unsigned INT row1, col1;
    unsigned INT row2, col2;
    
    register REAL *x1, *x2, *y1, *y2;
    register REAL *x3, *y3;
    
    INT i,j;
    INT start_row;
    INT start_col;
    
    switch (brow) {
            
        case 2:
            A11 = A->blocks[0];
            A12 = A->blocks[1];
            A21 = A->blocks[2];
            A22 = A->blocks[3];
            
            row1 = A11->row;
            col1 = A11->col;
            
            x1 = x;
            x2 = &(x[col1]);
            y1 = y;
            y2 = &(y[row1]);
            
            // y1 = A11*x1 + A12*x2
            fasp_blas_dcsr_mxv(A11, x1, y1);
            fasp_blas_dcsr_aAxpy(1.0, A12, x2, y1);
            
            // y2 = A21*x1 + A22*x2
            fasp_blas_dcsr_mxv(A21, x1, y2);
            fasp_blas_dcsr_aAxpy(1.0, A22, x2, y2);
            
            break;
            
        case 4:
            A11 = A->blocks[0];
            A12 = A->blocks[1];
            A13 = A->blocks[2];
            A21 = A->blocks[3];
            A22 = A->blocks[4];
            A23 = A->blocks[5];
            A31 = A->blocks[6];
            A32 = A->blocks[7];
            A33 = A->blocks[8];
            
            row1 = A11->row;
            col1 = A11->col;
            row2 = A22->row;
            col2 = A22->col;
            
            x1 = x;
            x2 = &(x[col1]);
            x3 = &(x[col1+col2]);
            y1 = y;
            y2 = &(y[row1]);
            y3 = &(y[row1+row2]);
            
            // y1 = A11*x1 + A12*x2 + A13*x3 + y1
            fasp_blas_dcsr_mxv(A11, x1, y1);
            fasp_blas_dcsr_aAxpy(1.0, A12, x2, y1);
            fasp_blas_dcsr_aAxpy(1.0, A13, x3, y1);
            
            // y2 = A21*x1 + A22*x2 + A23*x3 + y2
            fasp_blas_dcsr_mxv(A21, x1, y2);
            fasp_blas_dcsr_aAxpy(1.0, A22, x2, y2);
            fasp_blas_dcsr_aAxpy(1.0, A23, x3, y2);
            
            // y3 = A31*x1 + A32*x2 + A33*x3 + y2
            fasp_blas_dcsr_mxv(A31, x1, y3);
            fasp_blas_dcsr_aAxpy(1.0, A32, x2, y3);
            fasp_blas_dcsr_aAxpy(1.0, A33, x3, y3);
            
            break;
            
        default:
            
            start_row = 0;
            start_col = 0;
            
            for (i=0; i<brow; i++) {
                
                for (j=0; j<brow; j++){
                    
                    if (j==0) {
                        if (A->blocks[i*brow+j]){
                            fasp_blas_dcsr_mxv(A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
                        }
                    }
                    else {
                        if (A->blocks[i*brow+j]){
                            fasp_blas_dcsr_aAxpy(1.0, A->blocks[i*brow+j], &(x[start_col]), &(y[start_row]));
                        }
                    }
                    start_col = start_col + A->blocks[j*brow+j]->col;
                }
                
                start_row = start_row + A->blocks[i*brow+i]->row;
                start_col = 0;
            }
            
            break;
            
    } // end of switch
    
}

/**
 * \fn void fasp_blas_bdbsr_aAxpy (const REAL alpha, block_BSR *A, 
 *                                 REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  REAL factor a
 * \param A      Pointer to block_BSR matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Xiaozhe Hu
 * \date   11/11/2010
 */
void fasp_blas_bdbsr_aAxpy (const REAL alpha,
                            block_BSR *A,
                            REAL *x,
                            REAL *y)
{
    register dBSRmat *Arr = &(A->ResRes);
    register dCSRmat *Arw = &(A->ResWel);
    register dCSRmat *Awr = &(A->WelRes);
    register dCSRmat *Aww = &(A->WelWel);
    
    unsigned INT Nr = Arw->row;
    
    register REAL *xr = x;
    register REAL *xw = &(x[Nr]);
    register REAL *yr = y;
    register REAL *yw = &(y[Nr]);
    
    // yr = alpha*Arr*xr + alpha*Arw*xw + yr
    fasp_blas_dbsr_aAxpy(alpha, Arr, xr, yr);
    fasp_blas_dcsr_aAxpy(alpha, Arw, xw, yr);
    
    // yw = alpha*Awr*xr + alpha*Aww*xw + yw
    fasp_blas_dcsr_aAxpy(alpha, Awr, xr, yw);
    fasp_blas_dcsr_aAxpy(alpha, Aww, xw, yw);
}

/**
 * \fn void fasp_blas_bdbsr_mxv (block_BSR *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A      Pointer to block_BSR matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Xiaozhe Hu
 * \date   11/11/2010
 */
void fasp_blas_bdbsr_mxv (block_BSR *A,
                          REAL *x,
                          REAL *y)
{
    register dBSRmat *Arr = &(A->ResRes);
    register dCSRmat *Arw = &(A->ResWel);
    register dCSRmat *Awr = &(A->WelRes);
    register dCSRmat *Aww = &(A->WelWel);
    
    unsigned INT Nr = Arw->row;
    
    register REAL *xr = x;
    register REAL *xw = &(x[Nr]);
    register REAL *yr = y;
    register REAL *yw = &(y[Nr]);
    
    // yr = alpha*Arr*xr + alpha*Arw*xw + yr
    fasp_blas_dbsr_mxv(Arr, xr, yr);
    fasp_blas_dcsr_aAxpy(1.0, Arw, xw, yr); 
    
    // yw = alpha*Awr*xr + alpha*Aww*xw + yw
    fasp_blas_dcsr_mxv(Awr, xr, yw); 
    fasp_blas_dcsr_aAxpy(1.0, Aww, xw, yw); 
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
