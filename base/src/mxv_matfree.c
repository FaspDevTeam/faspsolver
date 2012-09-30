/*! \file mxv_matfree.c
 *  \brief Matrix-free version of itsolver.
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/**
 * \fn void fasp_blas_bdcsr_mxv (block_dCSRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to block_dCSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/18/2012
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
            
        case 3:
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
            
            // y1 = A11*x1 + A12*x2 + A13*x3
            fasp_blas_dcsr_mxv(A11, x1, y1);
            fasp_blas_dcsr_aAxpy(1.0, A12, x2, y1);
            fasp_blas_dcsr_aAxpy(1.0, A13, x3, y1);
            
            // y2 = A21*x1 + A22*x2 + A23*x3
            fasp_blas_dcsr_mxv(A21, x1, y2);
            fasp_blas_dcsr_aAxpy(1.0, A22, x2, y2);
            fasp_blas_dcsr_aAxpy(1.0, A23, x3, y2);
            
            // y3 = A31*x1 + A32*x2 + A33*x3
            fasp_blas_dcsr_mxv(A31, x1, y3);
            fasp_blas_dcsr_aAxpy(1.0, A32, x2, y3);
            fasp_blas_dcsr_aAxpy(1.0, A33, x3, y3);
            
            break;
            
    } // end of switch
    
}

/**
 * \fn void fasp_blas_dstr_mxv (dSTRmat *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to dSTR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/18/2012
 */
void fasp_blas_dstr_mxv (dSTRmat *A,
                         REAL *x,
                         REAL *y)
{
    int n = (A->ngrid)*(A->nc)*(A->nc);
    
    memset(y, 0, n*sizeof(REAL));
    
    fasp_blas_dstr_aAxpy(1.0, A, x, y);
}

/**
 * \fn void fasp_blas_mxv_csr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to CSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
void fasp_blas_mxv_csr (void *A,
                        REAL *x,
                        REAL *y)
{
    dCSRmat       *a = (dCSRmat *)A;
    fasp_blas_dcsr_mxv(a, x, y);
}

/**
 * \fn void fasp_blas_mxv_bsr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to BSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
void fasp_blas_mxv_bsr (void *A,
                        REAL *x,
                        REAL *y)
{
    dBSRmat       *a = (dBSRmat *)A;
    fasp_blas_dbsr_mxv(a, x, y);
}

/**
 * \fn void fasp_blas_mxv_str (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to STR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
void fasp_blas_mxv_str (void *A,
                        REAL *x,
                        REAL *y)
{
    dSTRmat       *a = (dSTRmat *)A;
    fasp_blas_dstr_mxv(a, x, y);
}

/**
 * \fn void fasp_blas_mxv_bcsr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to bCSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
void fasp_blas_mxv_bcsr (void *A,
                         REAL *x,
                         REAL *y)
{
    block_dCSRmat *a = (block_dCSRmat *)A;
    fasp_blas_bdcsr_mxv(a, x, y);
}

/**
 * \fn void fasp_blas_mxv_bbsr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to bBSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
void fasp_blas_mxv_bbsr (void *A,
                         REAL *x,
                         REAL *y)
{
    block_BSR     *a = (block_BSR *)A;
    fasp_blas_bdbsr_mxv(a, x, y);
}

/**
 * \fn void fasp_blas_mxv_csrl (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to CSRL matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
void fasp_blas_mxv_csrl (void *A,
                         REAL *x,
                         REAL *y)
{
    dCSRLmat      *a = (dCSRLmat *)A;
    fasp_blas_dcsrl_mxv(a, x, y);
}

/**
 * \fn void itsolver_init (int matrix_format, mxv_matfree *mf, void *A)
 *
 * \brief init itsovler.
 *
 * \param matrix_format    matrix format
 * \param mf               mxv_matfree action for iterative solvers
 * \param A                void pointer to matrix
 *
 * \author Feiteng Huang
 * \date   09/18/2012
 */
void itsolver_init (int matrix_format,
                    mxv_matfree *mf,
                    void *A)
{
    mf->data = A;
    
    switch ( matrix_format ) {
            
        case MAT_CSR:
            mf->fct = fasp_blas_mxv_csr;
            break;
            
        case MAT_BSR:
            mf->fct = fasp_blas_mxv_bsr;
            break;
            
        case MAT_STR:
            mf->fct = fasp_blas_mxv_str;
            break;
            
        case MAT_bCSR:
            mf->fct = fasp_blas_mxv_bcsr;
            break;
            
        case MAT_bBSR:
            mf->fct = fasp_blas_mxv_bbsr;
            break;
            
        case MAT_CSRL:
            mf->fct = fasp_blas_mxv_csrl;
            break;
            
        default:
            printf("### ERROR: Wrong matrix format!\n");
            exit(ERROR_DATA_STRUCTURE);
            
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/