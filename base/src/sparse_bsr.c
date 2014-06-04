/*! \file sparse_bsr.c
 *
 *  \brief Sparse matrix operations for dBSRmat matrices
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dBSRmat fasp_dbsr_create (INT ROW, INT COL, INT NNZ, INT nb, INT storage_manner)
 *
 * \brief Create BSR sparse matrix data memory space
 *
 * \param ROW             Number of rows of block
 * \param COL             Number of columns of block
 * \param NNZ             Number of nonzero blocks
 * \param nb              Dimension of exch block
 * \param storage_manner  Storage manner for each sub-block 
 *
 * \return A              The new dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010  
 */
dBSRmat fasp_dbsr_create (INT ROW, 
                          INT COL, 
                          INT NNZ, 
                          INT nb, 
                          INT storage_manner)
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
    
    A.ROW = ROW; A.COL = COL; A.NNZ = NNZ; 
    A.nb = nb; A.storage_manner = storage_manner;
    
    return A;
}

/**
 * \fn void fasp_dbsr_alloc (INT ROW, INT COL, INT NNZ, INT nb, INT storage_manner, 
 *                           dBSRmat *A)
 *
 * \brief Allocate memory space for BSR format sparse matrix
 *
 * \param ROW             Number of rows of block
 * \param COL             Number of columns of block
 * \param NNZ             Number of nonzero blocks
 * \param nb              Dimension of exch block
 * \param storage_manner  Storage manner for each sub-block 
 * \param A               Pointer to new dBSRmat matrix 
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010 
 */
void fasp_dbsr_alloc (INT ROW, 
                      INT COL, 
                      INT NNZ, 
                      INT nb, 
                      INT storage_manner, 
                      dBSRmat *A)
{    
    if ( ROW > 0 ) {
        A->IA = (INT*)fasp_mem_calloc(ROW+1, sizeof(INT));
    }
    else {
        A->IA = NULL;
    }
    
    if ( NNZ > 0 ) {
        A->JA = (INT*)fasp_mem_calloc(NNZ ,sizeof(INT));
    }
    else {
        A->JA = NULL;
    }    
    
    if ( nb > 0 ) {
        //A->val = (REAL*)fasp_mem_calloc(NNZ*nb*nb, sizeof(REAL));
		A->val = (REAL*)malloc(NNZ*nb*nb*sizeof(REAL));
    }
    else {
        A->val = NULL;
    }
    
    A->ROW = ROW; A->COL = COL; A->NNZ = NNZ; 
    A->nb = nb; A->storage_manner = storage_manner;
    
    return;
}


/**
 * \fn void fasp_dbsr_free (dBSRmat *A)
 *
 * \brief Free memeory space for BSR format sparse matrix
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
    A->nb = 0;
    A->storage_manner = 0;
}

/**
 * \fn void fasp_dbsr_null (dBSRmat *A)
 *
 * \brief Initialize sparse matrix on structured grid
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010 
 */
void fasp_dbsr_null (dBSRmat *A)
{    
    A->ROW=0;
    A->COL=0;
    A->NNZ=0;
    A->nb=0;
    A->storage_manner=0;
    A->IA=NULL;
    A->JA=NULL;
    A->val=NULL;
}

/**
 * \fn void fasp_dbsr_cp (dBSRmat *A, dBSRmat *B)
 *
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   Pointer to the dBSRmat matrix
 * \param B   Pointer to the dBSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date   08/07/2011  
 */
void fasp_dbsr_cp (dBSRmat *A, 
                   dBSRmat *B)
{    
    B->ROW=A->ROW;
    B->COL=A->COL;
    B->NNZ=A->NNZ;
    B->nb = A->nb;
    B->storage_manner = A->storage_manner;
    
    memcpy(B->IA,A->IA,(A->ROW+1)*sizeof(INT));
    memcpy(B->JA,A->JA,(A->NNZ)*sizeof(INT));
    memcpy(B->val,A->val,(A->NNZ)*(A->nb)*(A->nb)*sizeof(REAL));
}

/**
 * \fn INT fasp_dbsr_trans (dBSRmat *A, dBSRmat *AT)
 *
 * \brief Find A^T from given dBSRmat matrix A
 *
 * \param A   Pointer to the dBSRmat matrix
 * \param AT  Pointer to the transpose of dBSRmat matrix A
 *
 * \author Chunsheng FENG 
 * \date 2011/06/08
 *
 * Modified by Xiaozhe Hu (08/06/2011)
 */
INT fasp_dbsr_trans (dBSRmat *A, 
                     dBSRmat *AT)
{
    const INT n=A->ROW, m=A->COL, nnz=A->NNZ, nb=A->nb;    
    INT i,j,k,p,inb,jnb,nb2;
    INT status = SUCCESS;
    
    AT->ROW=m;
    AT->COL=n;
    AT->NNZ=nnz;
    AT->nb = nb;
    AT->storage_manner = A->storage_manner;
    
    AT->IA=(INT*)fasp_mem_calloc(m+1,sizeof(INT));
    nb2 =  A->nb*A->nb;

    AT->JA=(INT*)fasp_mem_calloc(nnz,sizeof(INT));
    
    if (A->val) {
        AT->val=(REAL*)fasp_mem_calloc(nnz*nb*nb,sizeof(REAL));     
    }
    else { 
        AT->val=NULL; 
    }
    
    // first pass: find the number of nonzeros in the first m-1 columns of A 
    // Note: these numbers are stored in the array AT.IA from 1 to m-1
    //r (i=0;i<m;++i) AT->IA[i] = 0;
    fasp_iarray_set(m+1, AT->IA, 0);

    for (j=0;j<nnz;++j) {
        i=N2C(A->JA[j]); // column number of A = row number of A'
        if (i<m-1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val) {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for (p=ibegin;p<iend1;p++) {
                j=A->JA[N2C(p)]+1;
                k=AT->IA[N2C(j)];
                AT->JA[N2C(k)]=C2N(i);
                for (inb=0;inb<nb;inb++)
                    for (jnb=0;jnb<nb;jnb++)
                        AT->val[ nb2*N2C(k) + inb*nb + jnb ] =A->val[nb2*N2C(p) + jnb*nb + inb ];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    
    }
    else {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
            for (p=ibegin;p<iend1;p++) {
                j=A->JA[N2C(p)]+1;
                k=AT->IA[N2C(j)];
                AT->JA[N2C(k)]=C2N(i);
                AT->IA[j]=k+1;
            } // end for p
        } // end of i
        
    } // end if 
    
    return (status);
}    

/*!
 * \fn SHORT fasp_dbsr_diagpref ( dBSRmat *A )
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
    SHORT         status = SUCCESS;
    const INT     num_rowsA = A -> ROW;
    const INT     num_colsA = A -> COL; 
    const INT     nb = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *A_i    = A -> IA;
    INT          *A_j    = A -> JA;
    REAL         *A_data = A -> val;
    
    INT   i, j, tempi, row_size;

#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend, ibegin, iend;
    INT nthreads = FASP_GET_NUM_THREADS();
#endif

    /* the matrix should be square */
    if (num_rowsA != num_colsA) return ERROR_INPUT_PAR;
    
    //REAL *tempd = (REAL*)fasp_mem_calloc(nb2, sizeof(REAL));
   
#ifdef _OPENMP
    if (num_rowsA > OPENMP_HOLDS) {
        REAL *tempd = (REAL*)fasp_mem_calloc(nb2*nthreads, sizeof(REAL));
//#pragma omp parallel for private (myid,mybegin,myend,i,j,tempi,tempd, ibegin,iend)
#pragma omp parallel for private (myid,mybegin,myend,i,j,tempi,ibegin,iend)
        for (myid = 0; myid < nthreads; myid++) {
            FASP_GET_START_END(myid, nthreads, num_rowsA, &mybegin, &myend);
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
                             //memcpy(tempd,     A_data+ibegin*nb2,  (nb2)*sizeof(REAL));
                             memcpy(tempd+myid*nb2,     A_data+ibegin*nb2,  (nb2)*sizeof(REAL));
                             memcpy(A_data+ibegin*nb2,  A_data+j*nb2,       (nb2)*sizeof(REAL));
                             //memcpy(A_data+j*nb2,     tempd,     (nb2)*sizeof(REAL));
                             memcpy(A_data+j*nb2,       tempd+myid*nb2,     (nb2)*sizeof(REAL));
                         }
                         break;
                    }
                    /* diagonal element is missing */
                    if (j == iend-1) {
                        status = -2;
                        break;
			//return -2;
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

    // free tempd
    //fasp_mem_free(tempd);
    
    if (status < 0) return status;
    else            return SUCCESS;
}

/**
 * \fn dvector fasp_dbsr_getdiaginv ( dBSRmat *A )
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
dvector fasp_dbsr_getdiaginv (dBSRmat *A)
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
    INT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;
    
#ifdef _OPENMP
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
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
            FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
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
            FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
            if (nb > 1) {
                for (i = mybegin; i < myend; ++i) {
                    fasp_blas_smat_inv(diaginv.val+i*nb2, nb);
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
                fasp_blas_smat_inv(diaginv.val+i*nb2, nb);
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

/**
 * \fn dBSRmat fasp_dbsr_diaginv ( dBSRmat *A )
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
dBSRmat fasp_dbsr_diaginv (dBSRmat *A)
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
    INT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;

#ifdef _OPENMP    
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
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
            FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);   
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
            FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
            if (nb > 1) {
                for (i = mybegin; i < myend; ++i) {
                    fasp_blas_smat_inv(diaginv+i*nb2, nb);
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
                fasp_blas_smat_inv(diaginv+i*nb2, nb);
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
            FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
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
 * \fn dBSRmat fasp_dbsr_diaginv2 ( dBSRmat *A, REAL *diaginv )
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
dBSRmat fasp_dbsr_diaginv2 (dBSRmat *A, 
                            REAL *diaginv)
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
    INT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;

#ifdef _OPENMP
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
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
            FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
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
 * \fn dBSRmat fasp_dbsr_diaginv3 (dBSRmat *A, REAL *diaginv)
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

dBSRmat fasp_dbsr_diaginv3 (dBSRmat *A, 
                            REAL *diaginv) 
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
    REAL *val = A->val;
    
    INT    *IAb  = NULL;
    INT    *JAb  = NULL;
    REAL *valb = NULL;
    
    INT nb2  = nb*nb;
    INT i,j,k,m;

    INT use_openmp = FALSE;

#ifdef _OPENMP  
    INT myid, mybegin, myend, stride_i, nthreads;
    if ( ROW > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
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
                    fasp_blas_smat_inv_nc2(diaginv+i*4);
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
                fasp_blas_smat_inv_nc2(diaginv+i*4);
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
                    fasp_blas_smat_inv_nc3(diaginv+i*9);
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
                fasp_blas_smat_inv_nc3(diaginv+i*9);
    
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
                    fasp_blas_smat_inv_nc5(diaginv+i*25);
    
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
                fasp_blas_smat_inv_nc5(diaginv+i*25);
    
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
                    fasp_blas_smat_inv_nc7(diaginv+i*49);
    
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
                fasp_blas_smat_inv_nc7(diaginv+i*49);
    
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
                    fasp_blas_smat_inv(diaginv+i*nb2, nb);
    
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
                fasp_blas_smat_inv(diaginv+i*nb2, nb);
    
                // compute D^{-1}*A
                for (k = IA[i]; k < IA[i+1]; ++k) {
                    m = k*nb2;
                    j = JA[k];
                    if (j != i) fasp_blas_smat_mul(diaginv+i*nb2, val+m, valb+m, nb);
                }
            }// end of main loop
        }
    
        break;
    }

    return (B);
}

/**
 * \fn dBSRmat fasp_dbsr_diaginv4 (dBSRmat *A, REAL *diaginv)
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
 * \note A is preordered that the first block of each row is the diagonal block!
 *
 * \author Xiaozhe Hu
 * \date 03/12/2011
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/26/2012
 */
dBSRmat fasp_dbsr_diaginv4 (dBSRmat *A, 
                            REAL *diaginv)
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
    INT nthreads = 1, use_openmp = FALSE;
    INT myid, mybegin, myend;

#ifdef _OPENMP
    if (ROW > OPENMP_HOLDS) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
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
                FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    ibegin = IA[i]; iend = IA[i+1];
                    // get the diagonal sub-blocks (It is the first block of each row)
                    m = ibegin*4;
                    memcpy(diaginv+i*4, val+m, 4*sizeof(REAL));
                    fasp_smat_identity_nc2(valb+m);
    
                    // compute the inverses of the diagonal sub-blocks 
                    //fasp_blas_smat_inv_nc2(diaginv+i*9);
                    fasp_blas_smat_inv_nc2(diaginv+i*4); // Modified by Zheng Li

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
                //fasp_blas_smat_inv_nc2(diaginv+i*9);
                fasp_blas_smat_inv_nc2(diaginv+i*4); // Modified by Zheng Li

                // compute D^{-1}*A
                for (k = ibegin+1; k < iend; ++k) {
                    m = k*4;
                    fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
                }
            }// end of main loop
        }
    
        break;    
    
    case 3:
        if (use_openmp) {
#ifdef _openmp
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, m, k)
#endif
            for (myid = 0; myid < nthreads; myid++) {
                FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    ibegin = IA[i]; iend = IA[i+1];
                    // get the diagonal sub-blocks (It is the first block of each row)
                    m = ibegin*9;
                    memcpy(diaginv+i*9, val+m, 9*sizeof(REAL));
                    fasp_smat_identity_nc3(valb+m);
                    // compute the inverses of the diagonal sub-blocks 
                    fasp_blas_smat_inv_nc3(diaginv+i*9);
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
                fasp_blas_smat_inv_nc3(diaginv+i*9);
    
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
                FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    // get the diagonal sub-blocks
                    ibegin = IA[i]; iend = IA[i+1];
                    m = ibegin*25;
                    memcpy(diaginv+i*25, val+m, 25*sizeof(REAL));
                    fasp_smat_identity_nc5(valb+m);
    
                    // compute the inverses of the diagonal sub-blocks 
                    fasp_blas_smat_inv_nc5(diaginv+i*25);
    
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
                fasp_blas_smat_inv_nc5(diaginv+i*25);
    
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
                FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    // get the diagonal sub-blocks
                    ibegin = IA[i]; iend = IA[i+1];    
                    m = ibegin*49;
                    memcpy(diaginv+i*49, val+m, 49*sizeof(REAL));
                    fasp_smat_identity_nc7(valb+m);
    
                    // compute the inverses of the diagonal sub-blocks 
                    fasp_blas_smat_inv_nc7(diaginv+i*49);
    
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
                fasp_blas_smat_inv_nc7(diaginv+i*49);
    
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
                FASP_GET_START_END(myid, nthreads, ROW, &mybegin, &myend);
                for (i = mybegin; i < myend; ++i) {
                    // get the diagonal sub-blocks
                    ibegin = IA[i]; iend = IA[i+1];
                    m = ibegin*nb2;
                    memcpy(diaginv+i*nb2, val+m, nb2*sizeof(REAL));
                    fasp_smat_identity(valb+m, nb, nb2);
    
                    // compute the inverses of the diagonal sub-blocks 
                    fasp_blas_smat_inv(diaginv+i*nb2, nb);
    
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
                fasp_blas_smat_inv(diaginv+i*nb2, nb);
    
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
 * \fn fasp_dbsr_getdiag ( INT n, dBSRmat *A, REAL *diag )
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
void fasp_dbsr_getdiag (INT n, 
                        dBSRmat *A, 
                        REAL *diag )
{
    const INT nb2 = A->nb*A->nb;
    
    INT i,k;

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
 * \fn dBSRmat fasp_dbsr_diagLU(dBSRmat *A, REAL *DL, REAL *DU)
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
dBSRmat fasp_dbsr_diagLU(dBSRmat *A,
                         REAL *DL,
                         REAL *DU)
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
    //REAL s;
    
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
                        
                        // compute s
                        //s = temp[3] - (temp[1]*temp[2])/temp[0];
                        
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
            printf("only works for nb = 2 or 3 now!!");
            break;
            
    } // end of switch
    
    // compute B = DL*A*DU
    switch (nb){
            
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
                        
                    } // end if (JA[j] == i)
                    
                } // end for (j=IA[i]; j<IA[i+1]; j++)
                
            } // end of for (i=0; i<ROW; i++)
            
            break;
            
        default:
            printf("only works for nb = 2 or 3 now!!");
            break;
            
    }
    
    // return
    return B;
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
