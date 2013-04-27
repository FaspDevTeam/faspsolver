/*! \file mxv_matfree.c
 *  \brief Matrix-free version of itsolver.
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

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