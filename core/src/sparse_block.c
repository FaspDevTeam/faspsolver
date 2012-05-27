/*! \file sparse_block.c
 *  \brief Functions and operation for block sparse matrices. 
 *
 */

#include <time.h>

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_dcsr_getblk (dCSRmat *A, INT *Is, INT *Js, INT m, INT n, dCSRmat *B)
 *
 * \brief Get a sub CSR matrix of A with specified rows and colums
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param B     Pointer to dCSRmat CSR matrix
 * \param Is    Pointer to selected rows
 * \param Js    Pointer to selected colums
 * \param m     Number of selected rows
 * \param n     Number of selected colums
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   12/25/2010
 *
 * Modified by Feiteng Huang on 5/25/2012: remove unnecessary loop
 */
SHORT fasp_dcsr_getblk (dCSRmat *A, 
                        INT *Is, 
                        INT *Js, 
                        INT m, 
                        INT n, 
                        dCSRmat *B)
{
    
    INT    i,j,k,nnz=0;
    SHORT  status = SUCCESS;
    
    // create colum flags
    INT *col_flag = (INT*)fasp_mem_calloc(A->col,sizeof(INT)); 
    
    B->row=m; B->col=n;
    
    // allocate memory
    B->IA  = (INT*)fasp_mem_calloc(m+1,sizeof(INT));
    B->JA  = (INT*)fasp_mem_calloc(A->nnz,sizeof(INT)); 
    B->val = (REAL*)fasp_mem_calloc(A->nnz,sizeof(REAL));
    
    for (i=0;i<n;++i) col_flag[N2C(Js[i])]=i+1;
    
    B->IA[0]=0;
    for (i=0;i<m;++i) {    
        for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
            j=A->JA[N2C(k)];
            if (col_flag[N2C(j)]>0) {
                B->JA[nnz]=col_flag[j]-1;
                B->val[nnz]=A->val[N2C(k)];
                nnz++;
            }
        } /* end for k */
        B->IA[i+1]=nnz;
    } /* end for i */
    
    B->nnz=nnz;
    B->JA=(INT*)fasp_mem_realloc(B->JA, sizeof(int)*nnz);
    B->val=(REAL*)fasp_mem_realloc(B->val, sizeof(REAL)*nnz);
    
    fasp_mem_free(col_flag);   
    
    return(status);
}

/**
 * \fn SHORT fasp_dbsr_getblk (dBSRmat *A, INT *Is, INT *Js, INT m, INT n, dBSRmat *B)
 *
 * \brief Get a sub BSR matrix of A with specified rows and columns. 
 *
 * \param A     Pointer to dBSRmat BSR matrix
 * \param B     Pointer to dBSRmat BSR matrix
 * \param Is    Pointer to selected rows
 * \param Js    Pointer to selected colums
 * \param m     Number of selected rows
 * \param n     Number of selected colums
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   12/25/2010
 */
SHORT fasp_dbsr_getblk (dBSRmat *A, 
                        INT *Is, 
                        INT *Js, 
                        INT m, 
                        INT n, 
                        dBSRmat *B)
{
    const INT nb = A->nb;
    const INT nb2=nb*nb;
    
    INT i,j,k,nnz=0;
    SHORT status = SUCCESS;
    
    // create colum flags
    INT *col_flag=(INT*)fasp_mem_calloc(A->COL,sizeof(INT)); 
    
    B->ROW=m; B->COL=n; B->nb=nb; B->storage_manner=A->storage_manner;  
    
    B->IA=(INT*)fasp_mem_calloc(m+1,sizeof(INT));
    
    for (i=0;i<n;++i) col_flag[N2C(Js[i])]=i+1;
    
    // first pass: count nonzeros for sub matrix
    B->IA[0]=0;
    for (i=0;i<m;++i){    
        for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
            j=A->JA[N2C(k)];
            if (col_flag[N2C(j)]>0) nnz++;
        } /* end for k */
        B->IA[i+1]=nnz;
    } /* end for i */
    B->NNZ=nnz;
    
    // allocate 
    B->JA=(INT*)fasp_mem_calloc(nnz,sizeof(INT)); 
    
    B->val=(REAL*)fasp_mem_calloc(nnz*nb2,sizeof(REAL));
    
    // second pass: copy data to B
    // no need to do the following loop, need to be modified!!  Xiaozhe 
    nnz = 0;
    for (i=0;i<m;++i){    
        for (k=A->IA[N2C(Is[i])];k<A->IA[N2C(Is[i])+1];++k) {
            j=A->JA[N2C(k)];
            if (col_flag[N2C(j)]>0) {
                B->JA[nnz]=col_flag[j]-1;
                memcpy(B->val+nnz*nb2, A->val+N2C(k)*nb2, nb2*sizeof(REAL));
                nnz++;
            }
        } /* end for k */
    } /* end for i */
    
    fasp_mem_free(col_flag);   
    
    return(status);
}

/**
 * \fn dCSRmat fasp_dbsr_getblk_dcsr (dBSRmat *A)
 *
 * \brief get dCSRmat block from a dBSRmat matrix 
 * 
 * \param *A   pointer to the BSR format matrix
 *
 * \return     dCSRmat matrix if succeed, NULL if fail
 *
 * \author Xiaozhe Hu
 * \date   03/16/2012 
 *
 * \note Required by the reservoir simulation package fasp4monix!
 *
 */
dCSRmat fasp_dbsr_getblk_dcsr (dBSRmat *A)
{
    // information about A
    const INT ROW = A->ROW;
    const INT COL = A->COL;
    const INT NNZ = A->NNZ;
    const SHORT nc = A->nb;
    const INT nc2 = nc*nc;
    
    REAL *val = A->val;
    INT *IA = A->IA;
    INT *JA = A->JA;
    
    // Pressure block
    dCSRmat P_csr = fasp_dcsr_create(ROW, COL, NNZ);
    REAL *Pval=P_csr.val;
    
    // local variable
    INT i,j;
    SHORT status = SUCCESS;
    
    // get pressure block
    memcpy(P_csr.JA, JA, NNZ*sizeof(INT)); 
    memcpy(P_csr.IA, IA, (ROW+1)*sizeof(INT));
    
    //for (i=NNZ, j=NNZ*nc2-nc2; i--; j-=nc2)
    for (i=NNZ, j=NNZ*nc2-nc2 + (0*nc+0); i--; j-=nc2) {
            Pval[i] = val[j];
        }
    
    // compress CSR format 
    status = fasp_dcsr_compress_inplace(&P_csr,1e-8);
    
    // return P 
    return P_csr;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
