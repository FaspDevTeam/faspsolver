/*! \file sparse_bsr.c
 *  \brief Simple operations for BSR format
 *
 */

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
 * \fn void fasp_dbsr_alloc (INT ROW, INT COL, INT NNZ, INT nb, INT storage_manner, dBSRmat *A)
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
        A->val = (REAL*)fasp_mem_calloc(NNZ*nb*nb, sizeof(REAL));
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
    for (i=0;i<m;++i) AT->IA[i] = 0;
    
    for (j=0;j<nnz;++j) {
        i=N2C(A->JA[j]); // column number of A = row number of A'
        if (i<m-1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val) {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
    
            for(p=ibegin;p<iend1;p++) {
                j=A->JA[N2C(p)]+1;
                k=AT->IA[N2C(j)];
                AT->JA[N2C(k)]=C2N(i);
                
                for(inb=0;inb<nb;inb++)
                    for (jnb=0;jnb<nb;jnb++)
                        AT->val[ nb2*N2C(k) + inb*nb + jnb ] =A->val[nb2*N2C(p) + jnb*nb + inb ];
                
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    
    }
    else {
        for (i=0;i<n;++i) {
            INT ibegin=A->IA[i], iend1=A->IA[i+1];
    
            for(p=ibegin;p<iend1;p++) {
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
 * \note Reordering is done in place. 
 *
 */
SHORT fasp_dbsr_diagpref (dBSRmat *A)
{    
    const INT     num_rowsA = A -> ROW; 
    const INT     num_colsA = A -> COL; 
    const INT     nb = A->nb;
    const INT      nb2 = nb*nb;
    const INT    *A_i    = A -> IA;
    INT          *A_j    = A -> JA;
    REAL         *A_data = A -> val;
    
    INT   i, j, tempi, row_size;

    /* the matrix should be square */
    if (num_rowsA != num_colsA) return ERROR_INPUT_PAR;
    
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
    
    // free tempd
    fasp_mem_free(tempd);
    
    return SUCCESS;
}

/**
 * \fn dBSRmat fasp_dbsr_diaginv ( dBSRmat *A )
 *
 * \brief Compute B := D^{-1}*A, where 'D' is the block diagonal part of A.
 *
 * \param A   Pointer to the dBSRmat matrix
 *
 * \author Zhou Zhiyang
 * \date   2010/10/26
 *
 * \note Works for general nb (Xiaozhe)
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
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;        
    
    memcpy(IAb, IA, (ROW+1)*sizeof(INT));
    memcpy(JAb, JA, NNZ*sizeof(INT));
    
    // allocate memory   
    diaginv = (REAL *)fasp_mem_calloc(size, sizeof(REAL));         
    
    // get all the diagonal sub-blocks   
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i+1]; ++k) {
            if (JA[k] == i)
                memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(REAL));
        }
    }
    
    // compute the inverses of all the diagonal sub-blocks   
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
    
    // compute D^{-1}*A
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
 * \author Zhou Zhiyang
 * \date  2010/11/07
 *
 * \note Works for general nb (Xiaozhe)
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
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;        
    
    memcpy(IAb, IA, (ROW+1)*sizeof(INT));
    memcpy(JAb, JA, NNZ*sizeof(INT));
    
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
 * \note Works for general nb (Xiaozhe)
 *
 * \author Xiaozhe Hu
 * \date   12/25/2010
 * \date   05/26/2012
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
    INT myid;
    INT mybegin;
    INT stride_i;
    INT myend;

	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || ROW <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
                
	    #pragma omp parallel
		nthreads = omp_get_num_threads();
	}
    
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
            stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                mybegin = myid*stride_i;
                if(myid < nthreads-1)  myend = mybegin+stride_i;
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
            stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                mybegin = myid*stride_i;
                if(myid < nthreads-1)  myend = mybegin+stride_i;
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
            stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                mybegin = myid*stride_i;
                if(myid < nthreads-1)  myend = mybegin+stride_i;
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
            stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                mybegin = myid*stride_i;
                if(myid < nthreads-1)  myend = mybegin+stride_i;
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
            stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend,i,k,m,j) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                mybegin = myid*stride_i;
                if(myid < nthreads-1)  myend = mybegin+stride_i;
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
 * \author Xiaozhe Hu
 * \date 03/12/2011
 *
 * \note Works for general nb (Xiaozhe)
 * \note A is preordered that the first block of each row is the diagonal block!!
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
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;        
    
    memcpy(IAb, IA, (ROW+1)*sizeof(INT));
    memcpy(JAb, JA, NNZ*sizeof(INT));
    
    switch (nb) {

    case 2:
        // main loop 
        for (i = 0; i < ROW; ++i) {
            ibegin = IA[i]; iend = IA[i+1];
            // get the diagonal sub-blocks (It is the first block of each row)
            m = ibegin*4;
            memcpy(diaginv+i*4, val+m, 4*sizeof(REAL));
            fasp_smat_identity_nc2(valb+m);
    
            // compute the inverses of the diagonal sub-blocks 
            fasp_blas_smat_inv_nc2(diaginv+i*9);
    
            // compute D^{-1}*A
            for (k = ibegin+1; k < iend; ++k) {
                m = k*4;
                fasp_blas_smat_mul_nc2(diaginv+i*4, val+m, valb+m);
            }
        }// end of main loop
    
        break;    
    
    case 3:
        // main loop 
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
    
        break;
    
    case 5: 
        // main loop 
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
    
        break;
    
    case 7:
        // main loop 
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
    
        break;
    
    default:
        // main loop 
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
 * \author Zhou Zhiyang
 * \date   2010/10/26 
 *
 * \note Works for general nb (Xiaozhe)
 */
void fasp_dbsr_getdiag (INT n, 
                        dBSRmat *A, 
                        REAL *diag )
{
    const INT nb2 = A->nb*A->nb;
    
    INT i,k;
    
    for (i = 0; i < n; ++i) {
        for (k = A->IA[i]; k < A->IA[i+1]; ++k) {
            if (A->JA[k] == i) {
                memcpy(diag+i*nb2, A->val+k*nb2, nb2*sizeof(REAL));
                break;
            }
        }
    }
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
