/*! \file  BlaSpmvMatFree.inl
 *
 *  \brief BLAS operations when matrix is implicit or its format is not specified
 *
 *  \note  This file contains Level-1 (Bla) functions, which are used in:
 *         SolMatFree.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void fasp_blas_mxv_csr (const void *A, const REAL *x, REAL *y)
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
static void fasp_blas_mxv_csr (const void *A,
                               const REAL *x,
                               REAL       *y)
{
    const dCSRmat *a = (const dCSRmat *)A;
    fasp_blas_dcsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_bsr (const void *A, const REAL *x, REAL *y)
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
static void fasp_blas_mxv_bsr (const void *A,
                               const REAL *x,
                               REAL       *y)
{
    const dBSRmat *a = (const dBSRmat *)A;
    fasp_blas_dbsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_str (const void *A, const REAL *x, REAL *y)
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
static void fasp_blas_mxv_str (const void *A,
                               const REAL *x,
                               REAL       *y)
{
    const dSTRmat *a = (const dSTRmat *)A;
    fasp_blas_dstr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_blc (const void *A, const REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to BLC matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_blc (const void *A,
                               const REAL *x,
                               REAL       *y)
{
    const dBLCmat *a = (const dBLCmat *)A;
    fasp_blas_dblc_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_csrl (const void *A, const REAL *x, REAL *y)
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
static void fasp_blas_mxv_csrl (const void *A,
                                const REAL *x,
                                REAL       *y)
{
    const dCSRLmat *a = (const dCSRLmat *)A;
    fasp_blas_dcsrl_mxv(a, x, y);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
