/*! \file BlaSparseBLC.c
 *
 *  \brief Sparse matrix block operations
 *
 *  \note This file contains Level-1 (Bla) functions. It requires
 *        AuxMemory.c and BlaSparseCSR.c
 */

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_dblc_free (dBLCmat *A)
 *
 * \brief Free block CSR sparse matrix data memory space
 *
 * \param A   Pointer to the dBLCmat matrix
 *
 * \author Xiaozhe Hu
 * \date   04/18/2014
 */
void fasp_dblc_free (dBLCmat *A)
{
    if (A == NULL) return; // Nothing need to be freed!
    
    INT i;
    INT num_blocks = (A->brow)*(A->bcol);
    
    for ( i=0; i<num_blocks; i++ ) {
        fasp_dcsr_free(A->blocks[i]);
        A->blocks[i] = NULL;
    }
    
    fasp_mem_free(A->blocks);
    A->blocks = NULL;
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
    col_flag=(INT*)fasp_mem_calloc(A->COL,sizeof(INT));
    
    B->ROW=m; B->COL=n; B->nb=nb; B->storage_manner=A->storage_manner;
    
    B->IA=(INT*)fasp_mem_calloc(m+1,sizeof(INT));
    
    if (use_openmp) {
#ifdef _OPENMP
        stride_i = n/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            mybegin = myid*stride_i;
            if (myid < nthreads-1)  myend = mybegin+stride_i;
            else myend = n;
            for (i=mybegin; i < myend; ++i) {
                col_flag[Js[i]]=i+1;
            }
        }
#endif
    }
    else {
        for (i=0;i<n;++i) col_flag[Js[i]]=i+1;
    }
    
    // first pass: count nonzeros for sub matrix
    B->IA[0]=0;
    for (i=0;i<m;++i) {
        for (k=A->IA[Is[i]];k<A->IA[Is[i]+1];++k) {
            j=A->JA[k];
            if (col_flag[j]>0) nnz++;
        } /* end for k */
        B->IA[i+1]=nnz;
    } /* end for i */
    B->NNZ=nnz;
    
    // allocate
    B->JA=(INT*)fasp_mem_calloc(nnz,sizeof(INT));
    
    B->val=(REAL*)fasp_mem_calloc(nnz*nb2,sizeof(REAL));
    
    // second pass: copy data to B
    // TODO: No need to do the following loop, need to be modified!!  --Xiaozhe
    nnz = 0;
    for (i=0;i<m;++i) {
        for (k=A->IA[Is[i]];k<A->IA[Is[i]+1];++k) {
            j=A->JA[k];
            if (col_flag[j]>0) {
                B->JA[nnz]=col_flag[j]-1;
                memcpy(B->val+nnz*nb2, A->val+k*nb2, nb2*sizeof(REAL));
                nnz++;
            }
        } /* end for k */
    } /* end for i */
    
    fasp_mem_free(col_flag);
	
    return(status);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
