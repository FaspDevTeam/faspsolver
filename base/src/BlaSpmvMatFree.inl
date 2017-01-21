/*! \file BlaSpmvMatFree.inl
 *
 *  \brief BLAS2 operations when matrix is implicit or its format is not specified
 *
 *  \note This file contains Level-1 (Bla) functions, which are used in
 *        SolMatFree.c
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void fasp_blas_mxv_csr (void *A, REAL *x, REAL *y)
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
static void fasp_blas_mxv_csr (void *A,
                               REAL *x,
                               REAL *y)
{
    dCSRmat *a = (dCSRmat *)A;
    fasp_blas_dcsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_bsr (void *A, REAL *x, REAL *y)
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
static void fasp_blas_mxv_bsr (void *A,
                               REAL *x,
                               REAL *y)
{
    dBSRmat *a = (dBSRmat *)A;
    fasp_blas_dbsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_str (void *A, REAL *x, REAL *y)
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
static void fasp_blas_mxv_str (void *A,
                               REAL *x,
                               REAL *y)
{
    dSTRmat *a = (dSTRmat *)A;
    fasp_blas_dstr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_blc (void *A, REAL *x, REAL *y)
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
static void fasp_blas_mxv_blc (void *A,
                               REAL *x,
                               REAL *y)
{
    dBLCmat *a = (dBLCmat *)A;
    fasp_blas_dblc_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_csrl (void *A, REAL *x, REAL *y)
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
static void fasp_blas_mxv_csrl (void *A,
                                REAL *x,
                                REAL *y)
{
    dCSRLmat *a = (dCSRLmat *)A;
    fasp_blas_dcsrl_mxv(a, x, y);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
