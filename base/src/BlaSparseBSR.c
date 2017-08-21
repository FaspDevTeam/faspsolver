/*! \file  BlaSparseBSR.c
 *
 *  \brief Sparse matrix operations for dBSRmat matrices
 *
 *  \note  This file contains Level-1 (Bla) functions. It requires:
 *         AuxArray.c, AuxMemory.c, AuxThreads.c, BlaSmallMat.c,
 *         and BlaSmallMatInv.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2009--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 */

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dBSRmat fasp_dbsr_create (const INT ROW, const INT COL, const INT NNZ,
 *                               const INT nb, const INT storage_manner)
 *
 * \brief Create BSR sparse matrix data memory space
 *
 * \param ROW             Number of rows of block
 * \param COL             Number of columns of block
 * \param NNZ             Number of nonzero blocks
 * \param nb              Dimension of each block
 * \param storage_manner  Storage manner for each sub-block
 *
 * \return A              The new dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
dBSRmat fasp_dbsr_create (const INT  ROW,
                          const INT  COL,
                          const INT  NNZ,
                          const INT  nb,
                          const INT  storage_manner)
{
    dBSRmat A;
    
    if ( ROW > 0 ) {
        A.IA = (INT*)fasp_mem_calloc(ROW+1, sizeof(INT));
    }
    else {
        A.IA = NULL;
    }
    
    if ( NNZ > 0 ) {
        A.JA = (INT*)fasp_mem_calloc(NNZ ,sizeof(INT));
    }
    else {
        A.JA = NULL;
    }
    
    if ( nb > 0 && NNZ > 0) {
        A.val = (REAL*)fasp_mem_calloc(NNZ*nb*nb, sizeof(REAL));
    }
    else {
        A.val = NULL;
    }
    
    A.storage_manner = storage_manner;
    A.ROW = ROW;
    A.COL = COL;
    A.NNZ = NNZ;
    A.nb  = nb;
    
    return A;
}

/**
 * \fn void fasp_dbsr_alloc (const INT ROW, const INT COL, const INT NNZ,
 *                           const INT nb, const INT storage_manner, dBSRmat *A)
 *
 * \brief Allocate memory space for BSR format sparse matrix
 *
 * \param ROW             Number of rows of block
 * \param COL             Number of columns of block
 * \param NNZ             Number of nonzero blocks
 * \param nb              Dimension of each block
 * \param storage_manner  Storage manner for each sub-block
 * \param A               Pointer to new dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
void fasp_dbsr_alloc (const INT  ROW,
                      const INT  COL,
                      const INT  NNZ,
                      const INT  nb,
                      const INT  storage_manner,
                      dBSRmat   *A)
{
    if ( ROW > 0 ) {
        A->IA = (INT*)fasp_mem_calloc(ROW+1, sizeof(INT));
    }
    else {
        A->IA = NULL;
    }
    
    if ( NNZ > 0 ) {
        A->JA = (INT*)fasp_mem_calloc(NNZ, sizeof(INT));
    }
    else {
        A->JA = NULL;
    }
    
    if ( nb > 0 ) {
        A->val = (REAL*)fasp_mem_calloc(NNZ*nb*nb, sizeof(REAL));
    }
    else {
        A->val = NULL;
    }
    
    A->storage_manner = storage_manner;
    A->ROW = ROW;
    A->COL = COL;
    A->NNZ = NNZ;
    A->nb  = nb;
    
    return;
}


/**
 * \fn void fasp_dbsr_free (dBSRmat *A)
 *
 * \brief Free memory space for BSR format sparse matrix
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
void fasp_dbsr_free (dBSRmat *A)
{
    if (A==NULL) return;
    
    fasp_mem_free(A->IA);
    fasp_mem_free(A->JA);
    fasp_mem_free(A->val);
    
    A->ROW = 0;
    A->COL = 0;
    A->NNZ = 0;
    A->nb  = 0;
    A->storage_manner = 0;
}

/**
 * \fn void fasp_dbsr_cp (const dBSRmat *A, dBSRmat *B)
 *
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   Pointer to the dBSRmat matrix
 * \param B   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   08/07/2011
 */
void fasp_dbsr_cp (const dBSRmat *A,
                   dBSRmat       *B)
{
    B->ROW = A->ROW;
    B->COL = A->COL;
    B->NNZ = A->NNZ;
    B->nb  = A->nb;
    B->storage_manner = A->storage_manner;
    
    memcpy(B->IA,A->IA,(A->ROW+1)*sizeof(INT));
    memcpy(B->JA,A->JA,(A->NNZ)*sizeof(INT));
    memcpy(B->val,A->val,(A->NNZ)*(A->nb)*(A->nb)*sizeof(REAL));
}

/**
 * \fn INT fasp_dbsr_trans (const dBSRmat *A, dBSRmat *AT)
 *
 * \brief Find A^T from given dBSRmat matrix A
 *
 * \param A   Pointer to the dBSRmat matrix
 * \param AT  Pointer to the transpose of dBSRmat matrix A
 *
 * \author Chunsheng FENG
 * \date   2011/06/08
 *
 * Modified by Xiaozhe Hu (08/06/2011)
 */
INT fasp_dbsr_trans (const dBSRmat *A,
                     dBSRmat       *AT)
{
    const INT n = A->ROW, m = A->COL, nnz = A->NNZ, nb = A->nb;
    
    INT status = FASP_SUCCESS;
    INT i,j,k,p,inb,jnb,nb2;
    
    AT->ROW = m;
    AT->COL = n;
    AT->NNZ = nnz;
    AT->nb  = nb;
    AT->storage_manner = A->storage_manner;
    
    AT->IA  = (INT*)fasp_mem_calloc(m+1,sizeof(INT));
    AT->JA  = (INT*)fasp_mem_calloc(nnz,sizeof(INT));
    nb2     = nb*nb;
    
    if (A->val) {
        AT->val = (REAL*)fasp_mem_calloc(nnz*nb2,sizeof(REAL));
    }
    else {
        AT->val = NULL;
    }
    
    // first pass: find the number of nonzeros in the first m-1 columns of A
    // Note: these numbers are stored in the array AT.IA from 1 to m-1
    fasp_iarray_set(m+1, AT->IA, 0);
    
    for ( j=0; j<nnz; ++j ) {
        i=A->JA[j]; // column number of A = row number of A'
        if (i<m-1) AT->IA[i+2]++;
    }
    
    for ( i=2; i<=m; ++i ) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if ( A->val ) {
        for ( i=0; i<n; ++i ) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for ( p=ibegin; p<iend1; p++ ) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                for ( inb=0; inb<nb; inb++ )
                    for ( jnb=0; jnb<nb; jnb++ )
                        AT->val[nb2*k + inb*nb + jnb] = A->val[nb2*p + jnb*nb + inb];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
        
    }
    else {
        for ( i=0; i<n; ++i ) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for ( p=ibegin; p<iend1; p++ ) {
                j=A->JA[p]+1;
                k=AT->IA[j];
                AT->JA[k]=i;
                AT->IA[j]=k+1;
            } // end for p
        } // end of i
        
    } // end if
    
    return (status);
}

/**
 * \fn SHORT fasp_dbsr_getblk (const dBSRmat *A, const INT *Is, const INT *Js,
 *                             const INT m, const INT n, dBSRmat *B)
 *
 * \brief Get a sub BSR matrix of A with specified rows and columns.
 *
 * \param A     Pointer to dBSRmat BSR matrix
 * \param B     Pointer to dBSRmat BSR matrix
 * \param Is    Pointer to selected rows
 * \param Js    Pointer to selected columns
 * \param m     Number of selected rows
 * \param n     Number of selected columns
 *
 * \return      FASP_SUCCESS if succeeded, otherwise return error information.
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   12/25/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
SHORT fasp_dbsr_getblk (const dBSRmat  *A,
                        const INT      *Is,
                        const INT      *Js,
                        const INT       m,
                        const INT       n,
                        dBSRmat        *B)
{
    INT   status = FASP_SUCCESS;
    INT   i,j,k,nnz=0;
    INT  *col_flag;
    SHORT use_openmp = FALSE;
    
    const INT nb = A->nb;
    const INT nb2=nb*nb;
    
#ifdef _OPENMP
    INT myid, mybegin, stride_i, myend, nthreads;
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    // create colum flags
    col_flag = (INT*)fasp_mem_calloc(A->COL,sizeof(INT));
    
    B->ROW=m; B->COL=n; B->nb=nb; B->storage_manner=A->storage_manner;
    
    B->IA = (INT*)fasp_mem_calloc(m+1,sizeof(INT));
    
    if ( use_openmp ) {
#ifdef _OPENMP
        stride_i = n/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            mybegin = myid*stride_i;
            if ( myid < nthreads-1 )  myend = mybegin+stride_i;
            else myend = n;
            for ( i = mybegin; i < myend; ++i ) {
                col_flag[Js[i]]=i+1;
            }
        }
#endif
    }
    else {
        for ( i=0; i<n; ++i ) col_flag[Js[i]]=i+1;
    }
    
    // first pass: count nonzeros for sub matrix
    B->IA[0] = 0;
    for ( i=0; i<m; ++i ) {
        for ( k=A->IA[Is[i]]; k<A->IA[Is[i]+1]; ++k ) {
            j=A->JA[k];
            if (col_flag[j]>0) nnz++;
        } /* end for k */
        B->IA[i+1] = nnz;
    } /* end for i */
    B->NNZ = nnz;
    
    // allocate
    B->JA  = (INT*)fasp_mem_calloc(nnz,sizeof(INT));
    B->val = (REAL*)fasp_mem_calloc(nnz*nb2,sizeof(REAL));
    
    // second pass: copy data to B
    // TODO: No need to do the following loop, need to be modified!!  --Xiaozhe
    nnz = 0;
    for ( i=0; i<m; ++i)  {
        for ( k=A->IA[Is[i]]; k<A->IA[Is[i]+1]; ++k ) {
            j = A->JA[k];
            if ( col_flag[j] > 0 ) {
                B->JA[nnz]=col_flag[j]-1;
                memcpy(B->val+nnz*nb2, A->val+k*nb2, nb2*sizeof(REAL));
                nnz++;
            }
        } /* end for k */
    } /* end for i */
    
    fasp_mem_free(col_flag);
    
    return(status);
}

/*!
 * \fn SHORT fasp_dbsr_diagpref (dBSRmat *A)
 *
 * \brief Reorder the column and data arrays of a square BSR matrix,
 *        so that the first entry in each row is the diagonal one.
 *
 * \param A   Pointer to the BSR matrix
 *
 * \author Xiaozhe Hu
 * \date   03/10/2011
 *
 * \author Chunsheng Feng, Zheng Li
 * \date   09/02/2012
 *
 * \note Reordering is done in place.
 */
SHORT fasp_dbsr_diagpref (dBSRmat *A)
{
    SHORT         status = FASP_SUCCESS;
    const INT     num_rowsA = A -> ROW;
    const INT     num_colsA = A -> COL;
    const INT     nb        = A->nb;
    const INT     nb2       = nb*nb;
    
    const INT    *A_i       = A -> IA;
    INT          *A_j       = A -> JA;
    REAL         *A_data    = A -> val;
    
    INT   i, j, tempi, row_size;
    
#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend, ibegin, iend;
    INT nthreads = fasp_get_num_threads();
#endif
    
    /* the matrix should be square */
    if (num_rowsA != num_colsA) return ERROR_INPUT_PAR;
    
#ifdef _OPENMP
    if (num_rowsA > OPENMP_HOLDS) {
        REAL *tempd = (REAL*)fasp_mem_calloc(nb2*nthreads, sizeof(REAL));
#pragma omp parallel for private (myid,mybegin,myend,i,j,tempi,ibegin,iend)
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, num_rowsA, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                ibegin = A_i[i+1]; iend = A_i[i];
                for (j = ibegin; j < iend; j ++) {
                    if (A_j[j] == i) {
                        if (j != ibegin) {
                            // swap index
                            tempi  = A_j[ibegin];
                            A_j[ibegin] = A_j[j];
                            A_j[j] = tempi;
                            // swap block
                            memcpy(tempd+myid*nb2,    A_data+ibegin*nb2, (nb2)*sizeof(REAL));
                            memcpy(A_data+ibegin*nb2, A_data+j*nb2,      (nb2)*sizeof(REAL));
                            memcpy(A_data+j*nb2,      tempd+myid*nb2,    (nb2)*sizeof(REAL));
                        }
                        break;
                    }
                    /* diagonal element is missing */
                    if (j == iend-1) {
                        status = -2;
                        break;
                    }
                }
            }
        }
        fasp_mem_free(tempd);
    }
    else {
#endif
        REAL *tempd = (REAL*)fasp_mem_calloc(nb2, sizeof(REAL));
        for (i = 0; i < num_rowsA; i ++) {
            row_size = A_i[i+1] - A_i[i];
            for (j = 0; j < row_size; j ++) {
                if (A_j[j] == i) {
                    if (j != 0) {
                        // swap index
                        tempi  = A_j[0];
                        A_j[0] = A_j[j];
                        A_j[j] = tempi;
                        // swap block
                        memcpy(tempd, A_data, (nb2)*sizeof(REAL));
                        memcpy(A_data, A_data+j*nb2, (nb2)*sizeof(REAL));
                        memcpy(A_data+j*nb2, tempd, (nb2)*sizeof(REAL));
                    }
                    break;
                }
                /* diagonal element is missing */
                if (j == row_size-1) return -2;
            }
            A_j    += row_size;
            A_data += row_size*nb2;
        }
        fasp_mem_free(tempd);
#ifdef _OPENMP
    }
#endif
    
    if (status < 0) return status;
    else            return FASP_SUCCESS;
}

/**
 * \fn dvector fasp_dbsr_getdiaginv (const dBSRmat *A)
 *
 * \brief Get D^{-1} of matrix A
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   02/19/2013
 *
 * \note Works for general nb (Xiaozhe)
 */
dvector fasp_dbsr_getdiaginv (const dBSRmat *A)
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb2;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;
    
    dvector diaginv;
    
    INT i,k;
    
    // Variables for OpenMP
    SHORT nthreads = 1, use_openmp = FALSE;
    INT   myid, mybegin, myend;
    
#ifdef _OPENMP
    if ( ROW > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    // allocate memory
    diaginv.row = size;
    diaginv.val = (REAL *)fasp_mem_calloc(size, sizeof(REAL));
    
    // get all the diagonal sub-blocks
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend, k)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) {
                for (k = IA[i]; k < IA[i+1]; ++k) {
                    if (JA[k] == i)
                        memcpy(diaginv.val+i*nb2, val+k*nb2, nb2*sizeof(REAL));
                }
            }
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {
            for (k = IA[i]; k < IA[i+1]; ++k) {
                if (JA[k] == i)
                    memcpy(diaginv.val+i*nb2, val+k*nb2, nb2*sizeof(REAL));
            }
        }
    }
    // compute the inverses of all the diagonal sub-blocks
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
            if (nb > 1) {
                for (i = mybegin; i < myend; ++i) {
                    fasp_smat_inv(diaginv.val+i*nb2, nb);
                }
            }
            else {
                for (i = mybegin; i < myend; ++i) {
                    // zero-diagonal should be tested previously
                    diaginv.val[i] = 1.0 / diaginv.val[i];
                }
            }
        }
    }
    else {
        if (nb > 1) {
            for (i = 0; i < ROW; ++i) {
                fasp_smat_inv(diaginv.val+i*nb2, nb);
            }
        }
        else {
            for (i = 0; i < ROW; ++i) {
                // zero-diagonal should be tested previously
                diaginv.val[i] = 1.0 / diaginv.val[i];
            }
        }
    }
    
    return (diaginv);
}

// TODO: Need to clean up fasp_dbsr_diaginv* functions

/**
 * \fn dBSRmat fasp_dbsr_diaginv (const dBSRmat *A)
 *
 * \brief Compute B := D^{-1}*A, where 'D' is the block diagonal part of A.
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Zhiyang Zhou
 * \date   2010/10/26
 *
 * \note Works for general nb (Xiaozhe)
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/25/2012
 */
dBSRmat fasp_dbsr_diaginv (const dBSRmat *A)
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     COL = A->COL;
    const INT     NNZ = A->NNZ;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb2;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;
    
    dBSRmat B;
    INT    *IAb  = NULL;
    INT    *JAb  = NULL;
    REAL   *valb = NULL;
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    REAL *diaginv = NULL;
    
    INT i,j,k,m,l;
    
    // Variables for OpenMP
    SHORT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;
    
#ifdef _OPENMP
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    memcpy(IAb, IA, (ROW+1)*sizeof(INT));
    memcpy(JAb, JA, NNZ*sizeof(INT));
    
    // allocate memory
    diaginv = (REAL *)fasp_mem_calloc(size, sizeof(REAL));
    
    // get all the diagonal sub-blocks
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend, k)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) {
                for (k = IA[i]; k < IA[i+1]; ++k) {
                    if (JA[k] == i)
                        memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(REAL));
                }
            }
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {
            for (k = IA[i]; k < IA[i+1]; ++k) {
                if (JA[k] == i)
                    memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(REAL));
            }
        }
    }
    
    // compute the inverses of all the diagonal sub-blocks
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
            if (nb > 1) {
                for (i = mybegin; i < myend; ++i) {
                    fasp_smat_inv(diaginv+i*nb2, nb);
                }
            }
            else {
                for (i = mybegin; i < myend; ++i) {
                    // zero-diagonal should be tested previously
                    diaginv[i] = 1.0 / diaginv[i];
                }
            }
        }
    }
    else {
        if (nb > 1) {
            for (i = 0; i < ROW; ++i) {
                fasp_smat_inv(diaginv+i*nb2, nb);
            }
        }
        else {
            for (i = 0; i < ROW; ++i) {
                // zero-diagonal should be tested previously
                diaginv[i] = 1.0 / diaginv[i];
            }
        }
    }
    
    // compute D^{-1}*A
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, k, m, j, l)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) {
                for (k = IA[i]; k < IA[i+1]; ++k) {
                    m = k*nb2;
                    j = JA[k];
                    if (j == i) {
                        // Identity sub-block
                        memset(valb+m, 0X0, nb2*sizeof(REAL));
                        for (l = 0; l < nb; l ++) valb[m+l*nb+l] = 1.0;
                    }
                    else {
                        fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                    }
                }
            }
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {
            for (k = IA[i]; k < IA[i+1]; ++k) {
                m = k*nb2;
                j = JA[k];
                if (j == i) {
                    // Identity sub-block
                    memset(valb+m, 0X0, nb2*sizeof(REAL));
                    for (l = 0; l < nb; l ++) valb[m+l*nb+l] = 1.0;
                }
                else {
                    fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                }
            }
        }
    }
    
    fasp_mem_free(diaginv);
    
    return (B);
}

/**
 * \fn dBSRmat fasp_dbsr_diaginv2 (const dBSRmat *A, REAL *diaginv)
 *
 * \brief Compute B := D^{-1}*A, where 'D' is the block diagonal part of A.
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param diaginv  Pointer to the inverses of all the diagonal blocks
 *
 * \author Zhiyang Zhou
 * \date  2010/11/07
 *
 * \note Works for general nb (Xiaozhe)
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/25/2012
 */
dBSRmat fasp_dbsr_diaginv2 (const dBSRmat *A,
                            REAL          *diaginv)
{
    // members of A
    const INT ROW = A->ROW;
    const INT COL = A->COL;
    const INT NNZ = A->NNZ;
    const INT nb  = A->nb, nbp1 = nb+1;
    const INT nb2 = nb*nb;
    
    INT    *IA  = A->IA;
    INT    *JA  = A->JA;
    REAL   *val = A->val;
    
    dBSRmat B;
    INT    *IAb  = NULL;
    INT    *JAb  = NULL;
    REAL   *valb = NULL;
    
    INT i,k,m,l,ibegin,iend;
    
    // Variables for OpenMP
    SHORT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;
    
#ifdef _OPENMP
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    memcpy(IAb, IA, (ROW+1)*sizeof(INT));
    memcpy(JAb, JA, NNZ*sizeof(INT));
    
    // compute D^{-1}*A
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private (myid, i, mybegin, myend, ibegin, iend, k, m, l)
#endif
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
            for (i = mybegin; i < myend; ++i) {
                ibegin = IA[i]; iend = IA[i+1];
                for (k = ibegin; k < iend; ++k) {
                    m = k*nb2;
                    if (JA[k] != i) {
                        fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                    }
                    else {
                        // Identity sub-block
                        memset(valb+m, 0X0, nb2*sizeof(REAL));
                        for (l = 0; l < nb; l ++) valb[m+l*nbp1] = 1.0;
                    }
                }
            }
        }
    }
    else {
        // compute D^{-1}*A
        for (i = 0; i < ROW; ++i) {
            ibegin = IA[i]; iend = IA[i+1];
            for (k = ibegin; k < iend; ++k) {
                m = k*nb2;
                if (JA[k] != i) {
                    fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                }
                else {
                    // Identity sub-block
                    memset(valb+m, 0X0, nb2*sizeof(REAL));
                    for (l = 0; l < nb; l ++) valb[m+l*nbp1] = 1.0;
                }
            }
        }
    }
    
    return (B);
}

/**
 * \fn dBSRmat fasp_dbsr_diaginv3 (const dBSRmat *A, REAL *diaginv)
 *
 * \brief Compute B := D^{-1}*A, where 'D' is the block diagonal part of A.
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param diaginv  Pointer to the inverses of all the diagonal blocks
 *
 * \return BSR matrix after diagonal scaling
 *
 * \author Xiaozhe Hu
 * \date   12/25/2010
 *
 * \note Works for general nb (Xiaozhe)
 *
 * Modified by Xiaozhe Hu on 05/26/2012
 */
dBSRmat fasp_dbsr_diaginv3 (const dBSRmat *A,
                            REAL          *diaginv)
{
    dBSRmat B;
    // members of A
    INT     ROW = A->ROW;
    INT     ROW_plus_one = ROW+1;
    INT     COL = A->COL;
    INT     NNZ = A->NNZ;
    INT     nb  = A->nb;
    INT    *IA  = A->IA;
    INT    *JA  = A->JA;
    REAL   *val = A->val;
    
    INT    *IAb  = NULL;
    INT    *JAb  = NULL;
    REAL   *valb = NULL;
    
    INT     nb2  = nb*nb;
    INT     i,j,k,m;
    
    SHORT   use_openmp = FALSE;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, stride_i, nthreads;
    if ( ROW > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    fasp_iarray_cp(ROW_plus_one, IA, IAb);
    fasp_iarray_cp(NNZ, JA, JAb);
    
    switch (nb) {
            
        case 2:
            // main loop
            if (use_openmp) {
#ifdef _OPENMP
                stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
                {
                    myid = omp_get_thread_num();
                    mybegin = myid*stride_i;
                    if (myid < nthreads-1)  myend = mybegin+stride_i;
                    else myend = ROW;
                    for (i=mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        
                        k = IA[i];
                        m = k*4;
                        memcpy(diaginv+i*4, val+m, 4*sizeof(REAL));
                        fasp_smat_identity_nc2(valb+m);
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc2(diaginv+i*4);
                        // compute D^{-1}*A
                        for (k = IA[i]+1; k < IA[i+1]; ++k)
                        {
                            m = k*4;
                            j = JA[k];
                            fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
                        }
                    }// end of main loop
                }
#endif
            }
            else {
                // main loop
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    k = IA[i];
                    m = k*4;
                    memcpy(diaginv+i*4, val+m, 4*sizeof(REAL));
                    fasp_smat_identity_nc2(valb+m);
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc2(diaginv+i*4);
                    // compute D^{-1}*A
                    for (k = IA[i]+1; k < IA[i+1]; ++k) {
                        m = k*4;
                        j = JA[k];
                        fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
                    }
                }// end of main loop
            }
            
            break;
            
        case 3:
            // main loop
            if (use_openmp) {
#ifdef _OPENMP
                stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
                {
                    myid = omp_get_thread_num();
                    mybegin = myid*stride_i;
                    if (myid < nthreads-1)  myend = mybegin+stride_i;
                    else myend = ROW;
                    for (i=mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            if (JA[k] == i) {
                                m = k*9;
                                memcpy(diaginv+i*9, val+m, 9*sizeof(REAL));
                                fasp_smat_identity_nc3(valb+m);
                            }
                        }
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc3(diaginv+i*9);
                        // compute D^{-1}*A
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            m = k*9;
                            j = JA[k];
                            if (j != i) fasp_blas_smat_mul_nc3(diaginv+i*9, val+m, valb+m);
                        }
                    }// end of main loop
                }
#endif
            }
            
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        if (JA[k] == i) {
                            m = k*9;
                            memcpy(diaginv+i*9, val+m, 9*sizeof(REAL));
                            fasp_smat_identity_nc3(valb+m);
                        }
                    }
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc3(diaginv+i*9);
                    
                    // compute D^{-1}*A
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        m = k*9;
                        j = JA[k];
                        if (j != i) fasp_blas_smat_mul_nc3(diaginv+i*9, val+m, valb+m);
                    }
                }// end of main loop
            }
            
            break;
            
        case 5:
            // main loop
            if (use_openmp) {
#ifdef _OPENMP
                stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
                {
                    myid = omp_get_thread_num();
                    mybegin = myid*stride_i;
                    if (myid < nthreads-1)  myend = mybegin+stride_i;
                    else myend = ROW;
                    for (i=mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            if (JA[k] == i) {
                                m = k*25;
                                memcpy(diaginv+i*25, val+m, 25*sizeof(REAL));
                                fasp_smat_identity_nc5(valb+m);
                            }
                        }
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc5(diaginv+i*25);
                        
                        // compute D^{-1}*A
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            m = k*25;
                            j = JA[k];
                            if (j != i) fasp_blas_smat_mul_nc5(diaginv+i*25, val+m, valb+m);
                        }
                    }// end of main loop
                }
#endif
            }
            
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        if (JA[k] == i)
                        {
                            m = k*25;
                            memcpy(diaginv+i*25, val+m, 25*sizeof(REAL));
                            fasp_smat_identity_nc5(valb+m);
                        }
                    }
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc5(diaginv+i*25);
                    
                    // compute D^{-1}*A
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        m = k*25;
                        j = JA[k];
                        if (j != i) fasp_blas_smat_mul_nc5(diaginv+i*25, val+m, valb+m);
                    }
                }// end of main loop
            }
            
            break;
            
        case 7:
            // main loop
            if (use_openmp) {
#ifdef _OPENMP
                stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
                {
                    myid = omp_get_thread_num();
                    mybegin = myid*stride_i;
                    if (myid < nthreads-1)  myend = mybegin+stride_i;
                    else myend = ROW;
                    for (i=mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            if (JA[k] == i) {
                                m = k*49;
                                memcpy(diaginv+i*49, val+m, 49*sizeof(REAL));
                                fasp_smat_identity_nc7(valb+m);
                            }
                        }
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc7(diaginv+i*49);
                        
                        // compute D^{-1}*A
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            m = k*49;
                            j = JA[k];
                            if (j != i) fasp_blas_smat_mul_nc7(diaginv+i*49, val+m, valb+m);
                        }
                    }// end of main loop
                }
#endif
            }
            
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        if (JA[k] == i) {
                            m = k*49;
                            memcpy(diaginv+i*49, val+m, 49*sizeof(REAL));
                            fasp_smat_identity_nc7(valb+m);
                        }
                    }
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc7(diaginv+i*49);
                    
                    // compute D^{-1}*A
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        m = k*49;
                        j = JA[k];
                        if (j != i) fasp_blas_smat_mul_nc7(diaginv+i*49, val+m, valb+m);
                    }
                }// end of main loop
            }
            
            break;
            
        default:
            // main loop
            if (use_openmp) {
#ifdef _OPENMP
                stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
                {
                    myid = omp_get_thread_num();
                    mybegin = myid*stride_i;
                    if (myid < nthreads-1)  myend = mybegin+stride_i;
                    else myend = ROW;
                    for (i=mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            if (JA[k] == i) {
                                m = k*nb2;
                                memcpy(diaginv+i*nb2, val+m, nb2*sizeof(REAL));
                                fasp_smat_identity(valb+m, nb, nb2);
                            }
                        }
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv(diaginv+i*nb2, nb);
                        
                        // compute D^{-1}*A
                        for (k = IA[i]; k < IA[i+1]; ++k) {
                            m = k*nb2;
                            j = JA[k];
                            if (j != i) fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                        }
                    }// end of main loop
                }
#endif
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        if (JA[k] == i) {
                            m = k*nb2;
                            memcpy(diaginv+i*nb2, val+m, nb2*sizeof(REAL));
                            fasp_smat_identity(valb+m, nb, nb2);
                        }
                    }
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv(diaginv+i*nb2, nb);
                    
                    // compute D^{-1}*A
                    for (k = IA[i]; k < IA[i+1]; ++k) {
                        m = k*nb2;
                        j = JA[k];
                        if (j != i) fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                    }
                } // end of main loop
            }
            
            break;
    }
    
    return (B);
}

/**
 * \fn dBSRmat fasp_dbsr_diaginv4 (const dBSRmat *A, REAL *diaginv)
 *
 * \brief Compute B := D^{-1}*A, where 'D' is the block diagonal part of A.
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param diaginv  Pointer to the inverses of all the diagonal blocks
 *
 * \return BSR matrix after diagonal scaling
 *
 * \note Works for general nb (Xiaozhe)
 *
 * \note A is pre-ordered that the first block of each row is the diagonal block!
 *
 * \author Xiaozhe Hu
 * \date 03/12/2011
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/26/2012
 */
dBSRmat fasp_dbsr_diaginv4 (const dBSRmat *A,
                            REAL          *diaginv)
{
    // members of A
    const INT     ROW = A->ROW;
    const INT     COL = A->COL;
    const INT     NNZ = A->NNZ;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;
    REAL         *val = A->val;
    
    dBSRmat B;
    INT    *IAb  = NULL;
    INT    *JAb  = NULL;
    REAL   *valb = NULL;
    
    INT i,k,m;
    INT ibegin, iend;
    
    // Variables for OpenMP
    SHORT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;
    
#ifdef _OPENMP
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    memcpy(IAb, IA, (ROW+1)*sizeof(INT));
    memcpy(JAb, JA, NNZ*sizeof(INT));
    
    switch (nb) {
            
        case 2:
            if (use_openmp) {
#ifdef _openmp
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, m, k)
#endif
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
                    for (i = mybegin; i < myend; ++i) {
                        ibegin = IA[i]; iend = IA[i+1];
                        // get the diagonal sub-blocks (It is the first block of each row)
                        m = ibegin*4;
                        memcpy(diaginv+i*4, val+m, 4*sizeof(REAL));
                        fasp_smat_identity_nc2(valb+m);
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc2(diaginv+i*4); // fixed by Zheng Li
                        
                        // compute D^{-1}*A
                        for (k = ibegin+1; k < iend; ++k) {
                            m = k*4;
                            fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
                        }
                    }
                }// end of main loop
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    ibegin = IA[i]; iend = IA[i+1];
                    // get the diagonal sub-blocks (It is the first block of each row)
                    m = ibegin*4;
                    memcpy(diaginv+i*4, val+m, 4*sizeof(REAL));
                    fasp_smat_identity_nc2(valb+m);
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc2(diaginv+i*4); // fixed by Zheng Li
                    
                    // compute D^{-1}*A
                    for (k = ibegin+1; k < iend; ++k) {
                        m = k*4;
                        fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
                    }
                } // end of main loop
            }
            
            break;
            
        case 3:
            if (use_openmp) {
#ifdef _openmp
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, m, k)
#endif
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
                    for (i = mybegin; i < myend; ++i) {
                        ibegin = IA[i]; iend = IA[i+1];
                        // get the diagonal sub-blocks (It is the first block of each row)
                        m = ibegin*9;
                        memcpy(diaginv+i*9, val+m, 9*sizeof(REAL));
                        fasp_smat_identity_nc3(valb+m);
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc3(diaginv+i*9);
                        // compute D^{-1}*A
                        for (k = ibegin+1; k < iend; ++k) {
                            m = k*9;
                            fasp_blas_smat_mul_nc3(diaginv+i*9, val+m, valb+m);
                        }
                    }
                }// end of main loop
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    ibegin = IA[i]; iend = IA[i+1];
                    // get the diagonal sub-blocks (It is the first block of each row)
                    m = ibegin*9;
                    memcpy(diaginv+i*9, val+m, 9*sizeof(REAL));
                    fasp_smat_identity_nc3(valb+m);
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc3(diaginv+i*9);
                    
                    // compute D^{-1}*A
                    for (k = ibegin+1; k < iend; ++k) {
                        m = k*9;
                        fasp_blas_smat_mul_nc3(diaginv+i*9, val+m, valb+m);
                    }
                }// end of main loop
            }
            
            break;
            
        case 5:
            if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, m, k)
#endif
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
                    for (i = mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        ibegin = IA[i]; iend = IA[i+1];
                        m = ibegin*25;
                        memcpy(diaginv+i*25, val+m, 25*sizeof(REAL));
                        fasp_smat_identity_nc5(valb+m);
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc5(diaginv+i*25);
                        
                        // compute D^{-1}*A
                        for (k = ibegin+1; k < iend; ++k) {
                            m = k*25;
                            fasp_blas_smat_mul_nc5(diaginv+i*25, val+m, valb+m);
                        }
                    }
                }
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    ibegin = IA[i]; iend = IA[i+1];
                    m = ibegin*25;
                    memcpy(diaginv+i*25, val+m, 25*sizeof(REAL));
                    fasp_smat_identity_nc5(valb+m);
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc5(diaginv+i*25);
                    
                    // compute D^{-1}*A
                    for (k = ibegin+1; k < iend; ++k) {
                        m = k*25;
                        fasp_blas_smat_mul_nc5(diaginv+i*25, val+m, valb+m);
                    }
                }// end of main loop
            }
            break;
            
        case 7:
            if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, i, mybegin, myend, ibegin, iend, m, k)
#endif
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
                    for (i = mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        ibegin = IA[i]; iend = IA[i+1];
                        m = ibegin*49;
                        memcpy(diaginv+i*49, val+m, 49*sizeof(REAL));
                        fasp_smat_identity_nc7(valb+m);
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv_nc7(diaginv+i*49);
                        
                        // compute D^{-1}*A
                        for (k = ibegin+1; k < iend; ++k) {
                            m = k*49;
                            fasp_blas_smat_mul_nc7(diaginv+i*49, val+m, valb+m);
                        }
                    }
                }// end of main loop
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    ibegin = IA[i]; iend = IA[i+1];
                    m = ibegin*49;
                    memcpy(diaginv+i*49, val+m, 49*sizeof(REAL));
                    fasp_smat_identity_nc7(valb+m);
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv_nc7(diaginv+i*49);
                    
                    // compute D^{-1}*A
                    for (k = ibegin+1; k < iend; ++k) {
                        m = k*49;
                        fasp_blas_smat_mul_nc7(diaginv+i*49, val+m, valb+m);
                    }
                }// end of main loop
            }
            
            break;
            
        default:
            if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, m, k)
#endif
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, ROW, &mybegin, &myend);
                    for (i = mybegin; i < myend; ++i) {
                        // get the diagonal sub-blocks
                        ibegin = IA[i]; iend = IA[i+1];
                        m = ibegin*nb2;
                        memcpy(diaginv+i*nb2, val+m, nb2*sizeof(REAL));
                        fasp_smat_identity(valb+m, nb, nb2);
                        
                        // compute the inverses of the diagonal sub-blocks
                        fasp_smat_inv(diaginv+i*nb2, nb);
                        
                        // compute D^{-1}*A
                        for (k = ibegin+1; k < iend; ++k) {
                            m = k*nb2;
                            fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                        }
                    }
                }// end of main loop
            }
            else {
                for (i = 0; i < ROW; ++i) {
                    // get the diagonal sub-blocks
                    ibegin = IA[i]; iend = IA[i+1];
                    m = ibegin*nb2;
                    memcpy(diaginv+i*nb2, val+m, nb2*sizeof(REAL));
                    fasp_smat_identity(valb+m, nb, nb2);
                    
                    // compute the inverses of the diagonal sub-blocks
                    fasp_smat_inv(diaginv+i*nb2, nb);
                    
                    // compute D^{-1}*A
                    for (k = ibegin+1; k < iend; ++k) {
                        m = k*nb2;
                        fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                    }
                } // end of main loop
            }
            
            break;
    }
    
    return (B);
}

/*!
 * \fn void fasp_dbsr_getdiag (INT n, const dBSRmat *A, REAL *diag)
 *
 * \brief Abstract the diagonal blocks of a BSR matrix.
 *
 * \param n     Number of blocks to get
 * \param A     Pointer to the 'dBSRmat' type matrix
 * \param diag  Pointer to array which stores the diagonal blocks in row by row manner
 *
 * \author Zhiyang Zhou
 * \date   2010/10/26
 *
 * \note Works for general nb (Xiaozhe)
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/25/2012
 */
void fasp_dbsr_getdiag (INT            n,
                        const dBSRmat *A,
                        REAL          *diag )
{
    const INT nb2 = A->nb*A->nb;
    
    INT i,k;
    
    if ( n==0 || n>A->ROW || n>A->COL ) n = MIN(A->ROW,A->COL);
    
#ifdef _OPENMP
#pragma omp parallel for private(i,k) if(n>OPENMP_HOLDS)
#endif
    for (i = 0; i < n; ++i) {
        for (k = A->IA[i]; k < A->IA[i+1]; ++k) {
            if (A->JA[k] == i) {
                memcpy(diag+i*nb2, A->val+k*nb2, nb2*sizeof(REAL));
                break;
            }
        }
    }
}

/**
 * \fn dBSRmat fasp_dbsr_diagLU (const dBSRmat *A, REAL *DL, REAL *DU)
 *
 * \brief Compute B := DL*A*DU. We decompose each diagonal block of A into LDU form
 *        and DL = diag(L^{-1}) and DU = diag(U^{-1}).
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param DL       Pointer to the diag(L^{-1})
 * \param DU       Pointer to the diag(U^{-1})
 *
 * \return BSR matrix after scaling
 *
 * \author Xiaozhe Hu
 * \date   04/02/2014
 */
dBSRmat fasp_dbsr_diagLU (const dBSRmat *A,
                          REAL          *DL,
                          REAL          *DU)
{
    
    // members of A
    const INT  ROW = A->ROW;
    const INT  ROW_plus_one = ROW+1;
    const INT  COL = A->COL;
    const INT  NNZ = A->NNZ;
    const INT  nb  = A->nb;
    
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;
    const REAL *val = A->val;
    
    INT  *IAb  = NULL;
    INT  *JAb  = NULL;
    REAL *valb = NULL;
    
    INT nb2  = nb*nb;
    INT i, j, k;
    
    // Create a dBSRmat 'B'
    dBSRmat B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    fasp_iarray_cp(ROW_plus_one, IA, IAb);
    fasp_iarray_cp(NNZ, JA, JAb);
    
    // work array
    REAL *temp = (REAL *)fasp_mem_calloc(nb2, sizeof(REAL));
    
    // get DL and DU
    switch (nb) {
            
        case 2:
            
            for (i=0; i<ROW; i++){
                
                for (j=IA[i]; j<IA[i+1]; j++){
                    
                    if (JA[j] == i){
                        
                        temp[0] = val[j*nb2];
                        temp[1] = val[j*nb2+1];
                        temp[2] = val[j*nb2+2];
                        temp[3] = val[j*nb2+3];
                        
                        // form DL
                        DL[i*nb2]   = 1.0;
                        DL[i*nb2+1] = 0.0;
                        DL[i*nb2+2] = -temp[2]/temp[0];
                        DL[i*nb2+3] = 1.0;
                        //DL[i*nb2+2] = -temp[2]/(temp[0]*s);
                        //DL[i*nb2+3] = 1.0/s;
                        
                        // form DU
                        DU[i*nb2]   = 1.0;
                        DU[i*nb2+1] = -temp[1]/temp[0];
                        DU[i*nb2+2] = 0.0;
                        DU[i*nb2+3] = 1.0;
                        
                        break;
                        
                    } // end of if (JA[j] == i)
                    
                } // end of for (j=IA[i]; j<IA[i+1]; j++)
                
            } // end of for (i=0; i<ROW; i++)
            
            break;
            
        case 3:
            
            for (i=0; i<ROW; i++){
                
                for (j=IA[i]; j<IA[i+1]; j++){
                    
                    if (JA[j] == i){
                        
                        temp[0] = val[j*nb2];
                        temp[1] = val[j*nb2+1];
                        temp[2] = val[j*nb2+2];
                        temp[3] = val[j*nb2+3];
                        temp[4] = val[j*nb2+4];
                        temp[5] = val[j*nb2+5];
                        temp[6] = val[j*nb2+6];
                        temp[7] = val[j*nb2+7];
                        temp[8] = val[j*nb2+8];
                        
                        // some auxiliry variables
                        REAL s22 = temp[4] - ((temp[1]*temp[3])/temp[0]);
                        REAL s23 = temp[5] - ((temp[2]*temp[3])/temp[0]);
                        REAL s32 = temp[7] - ((temp[1]*temp[6])/temp[0]);
                        
                        // form DL
                        DL[i*nb2]   = 1.0;
                        DL[i*nb2+1] = 0.0;
                        DL[i*nb2+2] = 0.0;
                        DL[i*nb2+3] = -temp[3]/temp[0];
                        DL[i*nb2+4] = 1.0;
                        DL[i*nb2+5] = 0.0;
                        DL[i*nb2+6] = -temp[6]/temp[0] + (temp[3]/temp[0])*(s32/s22);
                        DL[i*nb2+7] = -s32/s22;
                        DL[i*nb2+8] = 1.0;
                        
                        // form DU
                        DU[i*nb2]   = 1.0;
                        DU[i*nb2+1] = -temp[1]/temp[0];
                        DU[i*nb2+2] = -temp[2]/temp[0] + (temp[1]/temp[0])*(s23/s22);
                        DU[i*nb2+3] = 0.0;
                        DU[i*nb2+4] = 1.0;
                        DU[i*nb2+5] = -s23/s22;
                        DU[i*nb2+6] = 0.0;
                        DU[i*nb2+7] = 0.0;
                        DU[i*nb2+8] = 1.0;
                        
                        break;
                        
                    } // end of if (JA[j] == i)
                    
                } // end of for (j=IA[i]; j<IA[i+1]; j++)
                
            } // end of for (i=0; i<ROW; i++)
            
            break;
            
        default:
            printf("### ERROR: Only works for nb = 2 or 3 now!");
            break;
            
    } // end of switch
    
    // compute B = DL*A*DU
    switch (nb) {
            
        case 2:
            
            for (i=0; i<ROW; i++){
                
                for (j=IA[i]; j<IA[i+1]; j++){
                    
                    k = JA[j];
                    
                    // left multiply DL
                    fasp_blas_smat_mul_nc2(DL+i*nb2, val+j*nb2, temp);
                    
                    // right multiply DU
                    fasp_blas_smat_mul_nc2(temp, DU+k*nb2, valb+j*nb2);
                    
                    // for diagonal block, set it to be diagonal matrix
                    if (JA[j] == i){
                        
                        valb[j*nb2+1] = 0.0;
                        valb[j*nb2+2] = 0.0;
                        
                    } // end if (JA[j] == i)
                    
                    
                } // end for (j=IA[i]; j<IA[i+1]; j++)
                
            } // end of for (i=0; i<ROW; i++)
            
            break;
            
        case 3:
            
            for (i=0; i<ROW; i++){
                
                for (j=IA[i]; j<IA[i+1]; j++){
                    
                    k = JA[j];
                    
                    // left multiply DL
                    fasp_blas_smat_mul_nc3(DL+i*nb2, val+j*nb2, temp);
                    
                    // right multiply DU
                    fasp_blas_smat_mul_nc3(temp, DU+k*nb2, valb+j*nb2);
                    
                    // for diagonal block, set it to be diagonal matrix
                    if (JA[j] == i){
                        
                        valb[j*nb2+1] = 0.0;
                        valb[j*nb2+2] = 0.0;
                        valb[j*nb2+3] = 0.0;
                        valb[j*nb2+5] = 0.0;
                        valb[j*nb2+6] = 0.0;
                        valb[j*nb2+7] = 0.0;
                        if (ABS(valb[j*nb2+4]) < SMALLREAL) valb[j*nb2+4] = SMALLREAL;
                        if (ABS(valb[j*nb2+8]) < SMALLREAL) valb[j*nb2+8] = SMALLREAL;
                        
                    } // end if (JA[j] == i)
                    
                } // end for (j=IA[i]; j<IA[i+1]; j++)
                
            } // end of for (i=0; i<ROW; i++)
            
            break;
            
        default:
            printf("### ERROR: Only works for nb = 2 or 3 now!");
            break;
            
    }
    
    // return
    return B;
    
}

/**
 * \fn dBSRmat fasp_dbsr_diagLU2 (dBSRmat *A, REAL *DL, REAL *DU)
 *
 * \brief Compute B := DL*A*DU. We decompose each diagonal block of A into LDU form
 *        and DL = diag(L^{-1}) and DU = diag(U^{-1}).
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param DL       Pointer to the diag(L^{-1})
 * \param DU       Pointer to the diag(U^{-1})
 *
 * \return BSR matrix after scaling
 *
 * \author Zheng Li, Xiaozhe Hu
 * \date   06/17/2014
 */
dBSRmat fasp_dbsr_diagLU2 (dBSRmat *A,
                           REAL    *DL,
                           REAL    *DU)
{
    
    // members of A
    INT  ROW = A->ROW;
    INT  ROW_plus_one = ROW+1;
    INT  COL = A->COL;
    INT  NNZ = A->NNZ;
    INT  nb  = A->nb;
    INT  *IA  = A->IA;
    INT  *JA  = A->JA;
    REAL *val = A->val;
    
    INT  *IAb  = NULL;
    INT  *JAb  = NULL;
    REAL *valb = NULL;
    
    INT nb2  = nb*nb;
    INT i,j,k;
    
    REAL sqt3, sqt4, sqt8;
    
    // Create a dBSRmat 'B'
    dBSRmat B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    REAL *temp = (REAL *)fasp_mem_calloc(nb*nb, sizeof(REAL));
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    fasp_iarray_cp(ROW_plus_one, IA, IAb);
    fasp_iarray_cp(NNZ, JA, JAb);
    
    // get DL and DU
    switch (nb) {
        case 2:
            for (i=0; i<ROW; i++){
                for (j=IA[i]; j<IA[i+1]; j++){
                    if (JA[j] == i){
                        REAL temp0 = val[j*nb2];
                        REAL temp1 = val[j*nb2+1];
                        REAL temp2 = val[j*nb2+2];
                        REAL temp3 = val[j*nb2+3];
                        
                        if (ABS(temp3) < SMALLREAL) temp3 = SMALLREAL;
                        
                        sqt3 = sqrt(ABS(temp3));
                        
                        // form DL
                        DL[i*nb2]   = 1.0;
                        DL[i*nb2+1] = 0.0;
                        DL[i*nb2+2] = -temp2/temp0/sqt3;
                        DL[i*nb2+3] = 1.0/sqt3;
                        
                        // form DU
                        DU[i*nb2]   = 1.0;
                        DU[i*nb2+1] = -temp1/temp0/sqt3;
                        DU[i*nb2+2] = 0.0;
                        DU[i*nb2+3] = 1.0/sqt3;
                        break;
                        
                    }
                }
            }
            
            break;
            
        case 3:
            for (i=0; i<ROW; i++){
                for (j=IA[i]; j<IA[i+1]; j++){
                    if (JA[j] == i){
                        REAL temp0 = val[j*nb2];
                        REAL temp1 = val[j*nb2+1];
                        REAL temp2 = val[j*nb2+2];
                        REAL temp3 = val[j*nb2+3];
                        REAL temp4 = val[j*nb2+4];
                        REAL temp5 = val[j*nb2+5];
                        REAL temp6 = val[j*nb2+6];
                        REAL temp7 = val[j*nb2+7];
                        REAL temp8 = val[j*nb2+8];
                        
                        if (ABS(temp4) < SMALLREAL)  temp4 = SMALLREAL;
                        if (ABS(temp8) < SMALLREAL)  temp8 = SMALLREAL;
                        
                        sqt4 = sqrt(ABS(temp4));
                        sqt8 = sqrt(ABS(temp8));
                        
                        // some auxiliary variables
                        REAL s22 = temp4 - ((temp1*temp3)/temp0);
                        REAL s23 = temp5 - ((temp2*temp3)/temp0);
                        REAL s32 = temp7 - ((temp1*temp6)/temp0);
                        
                        // form DL
                        DL[i*nb2]   = 1.0;
                        DL[i*nb2+1] = 0.0;
                        DL[i*nb2+2] = 0.0;
                        DL[i*nb2+3] = -temp3/temp0/sqt4;
                        DL[i*nb2+4] = 1.0/sqt4;
                        DL[i*nb2+5] = 0.0;
                        DL[i*nb2+6] = (-temp6/temp0 + (temp3/temp0)*(s32/s22))/sqt8;
                        DL[i*nb2+7] = -s32/s22/sqt8;
                        DL[i*nb2+8] = 1.0/sqt8;
                        
                        // form DU
                        DU[i*nb2]   = 1.0;
                        DU[i*nb2+1] = -temp1/temp0/sqt4;
                        DU[i*nb2+2] = (-temp2/temp0 + (temp1/temp0)*(s23/s22))/sqt8;
                        DU[i*nb2+3] = 0.0;
                        DU[i*nb2+4] = 1.0/sqt4;
                        DU[i*nb2+5] = -s23/s22/sqt8;
                        DU[i*nb2+6] = 0.0;
                        DU[i*nb2+7] = 0.0;
                        DU[i*nb2+8] = 1.0/sqt8;
                        
                        break;
                        
                    }
                }
            }
            
            break;
            
        default:
            printf("### ERROR: Only works for nb = 2 or 3 now!");
            break;
            
    } // end of switch
    
    // compute B = DL*A*DU
    switch (nb) {
            
        case 2:
            for (i=0; i<ROW; i++){
                for (j=IA[i]; j<IA[i+1]; j++){
                    k = JA[j];
                    // left multiply DL
                    fasp_blas_smat_mul_nc2(DL+i*nb2, val+j*nb2, temp);
                    // right multiply DU
                    fasp_blas_smat_mul_nc2(temp, DU+k*nb2, valb+j*nb2);
                    // for diagonal block, set it to be diagonal matrix
                    if (JA[j] == i){
                        valb[j*nb2+1] = 0.0;
                        valb[j*nb2+2] = 0.0;
                        if (ABS(valb[j*nb2+3]) < SMALLREAL) valb[j*nb2+3] = SMALLREAL;
                    }
                }
            }
            
            break;
            
        case 3:
            for (i=0; i<ROW; i++){
                for (j=IA[i]; j<IA[i+1]; j++){
                    k = JA[j];
                    // left multiply DL
                    fasp_blas_smat_mul_nc3(DL+i*nb2, val+j*nb2, temp);
                    // right multiply DU
                    fasp_blas_smat_mul_nc3(temp, DU+k*nb2, valb+j*nb2);
                    // for diagonal block, set it to be diagonal matrix
                    if (JA[j] == i){
                        valb[j*nb2+1] = 0.0;
                        valb[j*nb2+2] = 0.0;
                        valb[j*nb2+3] = 0.0;
                        valb[j*nb2+5] = 0.0;
                        valb[j*nb2+6] = 0.0;
                        valb[j*nb2+7] = 0.0;
                        if (ABS(valb[j*nb2+4]) < SMALLREAL) valb[j*nb2+4] = SMALLREAL;
                        if (ABS(valb[j*nb2+8]) < SMALLREAL) valb[j*nb2+8] = SMALLREAL;
                    }
                }
            }
            break;
            
        default:
            printf("### ERROR: Only works for nb = 2 or 3 now!");
            break;
    }
    
    // return
    return B;
    
}

/**
 * \fn dBSRmat fasp_dbsr_perm (const dBSRmat *A, const INT *P)
 *
 * \brief Apply permutation of A, i.e. Aperm=PAP' by the orders given in P
 *
 * \param A  Pointer to the original dBSRmat matrix
 * \param P  Pointer to the given ordering
 *
 * \return   The new ordered dBSRmat matrix if succeed, NULL if fail
 *
 * \author Zheng Li
 * \date   24/9/2015
 *
 * \note   P[i] = k means k-th row and column become i-th row and column!
 */
dBSRmat fasp_dbsr_perm (const dBSRmat *A,
                        const INT     *P)
{
    const INT   n = A->ROW, nnz = A->NNZ;
    const INT  *ia= A->IA, *ja = A->JA;
    const REAL *Aval = A->val;
    const INT   nb = A->nb, nb2 = nb*nb;
    const INT   manner = A->storage_manner;
    SHORT       nthreads = 1, use_openmp = FALSE;
    
    INT i,j,k,jaj,i1,i2,start,jj;
    
#ifdef _OPENMP
    if ( MIN(n, nnz) > OPENMP_HOLDS ) {
        use_openmp = 0;//TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    dBSRmat Aperm = fasp_dbsr_create(n,n,nnz,nb,manner);
    
    // form the transpose of P
    INT *Pt = (INT*)fasp_mem_calloc(n,sizeof(INT));
    
    if (use_openmp) {
        INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i)
#endif
        for (myid=0; myid<nthreads; ++myid) {
            fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) Pt[P[i]] = i;
        }
    }
    else {
        for (i=0; i<n; ++i) Pt[P[i]] = i;
    }
    
    // compute IA of P*A (row permutation)
    Aperm.IA[0] = 0;
    for (i=0; i<n; ++i) {
        k = P[i];
        Aperm.IA[i+1] = Aperm.IA[i]+(ia[k+1]-ia[k]);
    }
    
    // perform actual P*A
    if (use_openmp) {
        INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, i1, i2, k, start, j, jaj, jj)
#endif
        for (myid=0; myid<nthreads; ++myid) {
            fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) {
                i1 = Aperm.IA[i]; i2 = Aperm.IA[i+1]-1;
                k = P[i];
                start = ia[k];
                for (j=i1; j<=i2; ++j) {
                    jaj = start+j-i1;
                    Aperm.JA[j] = ja[jaj];
                    for (jj=0; jj<nb2;++jj)
                        Aperm.val[j*nb2+jj] = Aval[jaj*nb2+jj];
                }
            }
        }
    }
    else {
        for (i=0; i<n; ++i) {
            i1 = Aperm.IA[i]; i2 = Aperm.IA[i+1]-1;
            k = P[i];
            start = ia[k];
            for (j=i1; j<=i2; ++j) {
                jaj = start+j-i1;
                Aperm.JA[j] = ja[jaj];
                for (jj=0; jj<nb2;++jj)
                    Aperm.val[j*nb2+jj] = Aval[jaj*nb2+jj];
            }
        }
    }
    
    // perform P*A*P' (column permutation)
    if (use_openmp) {
        INT myid, mybegin, myend;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, k, j)
#endif
        for (myid=0; myid<nthreads; ++myid) {
            fasp_get_start_end(myid, nthreads, nnz, &mybegin, &myend);
            for (k=mybegin; k<myend; ++k) {
                j = Aperm.JA[k];
                Aperm.JA[k] = Pt[j];
            }
        }
    }
    else {
        for (k=0; k<nnz; ++k) {
            j = Aperm.JA[k];
            Aperm.JA[k] = Pt[j];
        }
    }
    
    fasp_mem_free(Pt);
    
    return(Aperm);
}

/**
 * \fn INT fasp_dbsr_merge_col (const dBSRmat *A)
 *
 * \brief Check and merge some same col index in one row.
 *
 * \param A  Pointer to the original dCSRmat matrix
 *
 * \return   The new merged dCSRmat matrix
 *
 * \author Chunsheng Feng
 * \date   30/07/2017
 */
INT fasp_dbsr_merge_col (dBSRmat *A)
{
    INT           count = 0;
    const INT     num_rowsA = A -> ROW;
    const INT     nb = A->nb;
    const INT     nb2 = nb*nb;
    INT          *A_i    = A -> IA;
    INT          *A_j    = A -> JA;
    REAL         *A_data = A -> val;
    
    INT   i, ii, j, jj, ibegin, iend, iend1;
    
#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend;
    INT nthreads = fasp_get_num_threads();
#endif
    
#ifdef _OPENMP
    if (num_rowsA > OPENMP_HOLDS) {
#pragma omp parallel for private (myid,mybegin,myend,i,ii,j,jj,ibegin,iend,iend1)
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, num_rowsA, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                ibegin = A_i[i]; iend = A_i[i+1]; iend1 =  iend-1;
                for (j = ibegin; j < iend1; j ++) {
                    if (A_j[j] > -1) {
                        for (jj = j+1; jj < iend; jj ++) {
                            if (A_j[j] == A_j[jj]) {
                                // add jj col to j
                                for ( ii=0; ii <nb2; ii++ )
                                    A_data[j*nb2 +ii] += A_data[ jj*nb2+ii];
                                A_j[jj] = -1;
                                count ++;
                            }
                        }
                    }
                }
            }
        }
    }
    else {
#endif
        for (i = 0; i < num_rowsA; i ++) {
            ibegin = A_i[i]; iend = A_i[i+1]; iend1 =  iend-1;
            for (j = ibegin; j < iend1; j ++) {
                if (A_j[j] > -1) {
                    for (jj = j+1; jj < iend; jj ++) {
                        if (A_j[j] == A_j[jj]) {
                            // add jj col to j
                            for ( ii=0; ii <nb2; ii++ )
                                A_data[j*nb2 +ii] += A_data[ jj*nb2+ii];
                            printf("### WARNING: Same col indices at %d, col %d (%d %d)\n",
                                   i, A_j[j], j, jj );
                            A_j[jj] = -1;
                            count ++;
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
    }
#endif
    
    if ( count > 0 ) {
        INT *tempA_i = (INT*)fasp_mem_calloc(num_rowsA+1, sizeof(INT));
        memcpy(tempA_i, A_i, (num_rowsA+1 )*sizeof(INT));
        jj = 0; 	A_i[0] = jj;
        for (i = 0; i < num_rowsA; i ++) {
            ibegin = tempA_i[i]; iend = tempA_i[i+1];
            for (j = ibegin; j < iend; j ++) {
                if (A_j[j] > -1 ) {
                    memcpy(A_data +jj*nb2, A_data+j*nb2, (nb2)*sizeof(REAL));
                    A_j[jj] = A_j[j];
                    jj++;
                }
            }
            A_i[i+1]	= jj;
        }
        A-> NNZ = jj;
        fasp_mem_free(tempA_i);
        
        printf("### WARNING: %d col indices have been merged!\n", count);
    }
    
    return count;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
