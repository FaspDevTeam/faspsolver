/*! \file sparse_csr.c
 *  \brief Functions for CSR sparse matrices. 
 */

#include <math.h>
#include <omp.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dCSRmat fasp_dcsr_create (const INT m, const INT n, const INT nnz)
 *
 * \brief Create CSR sparse matrix data memory space
 *
 * \param m    Number of rows
 * \param n    Number of columns
 * \param nnz  Number of nonzeros
 *
 * \return A   the new dCSRmat matrix
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
 */
dCSRmat fasp_dcsr_create (const INT m, 
                          const INT n, 
                          const INT nnz)
{    
    dCSRmat A;
    
    if ( m > 0 ) {
        A.IA = (INT*)fasp_mem_calloc(m+1,sizeof(INT));
    }
    else {
        A.IA = NULL;
    }
    
    if ( n > 0 ) {
        A.JA = (INT*)fasp_mem_calloc(nnz,sizeof(INT));
    }
    else {
        A.JA = NULL;
    }    
    
    if ( nnz > 0 ) {
        A.val=(REAL*)fasp_mem_calloc(nnz,sizeof(REAL));
    }
    else {
        A.val = NULL;
    }
    
    A.row=m; A.col=n; A.nnz=nnz;
    
    return A;
}

/**
 * \fn void fasp_dcsr_alloc (const INT m, const INT n, const INT nnz, dCSRmat *A)
 *
 * \brief Allocate CSR sparse matrix memory space
 *
 * \param m      Number of rows
 * \param n      Number of columns
 * \param nnz    Number of nonzeros
 * \param A      Pointer to the dCSRmat matrix
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
 */
void fasp_dcsr_alloc (const INT m, 
                      const INT n, 
                      const INT nnz, 
                      dCSRmat *A)
{    
    if ( m > 0 ) {
        A->IA=(INT*)fasp_mem_calloc(m+1,sizeof(INT));
    }
    else {
        A->IA = NULL;
    }
    
    if ( n > 0 ) {
        A->JA=(INT*)fasp_mem_calloc(nnz,sizeof(INT));
    }
    else {
        A->JA = NULL;
    }    
    
    if ( nnz > 0 ) {
        A->val=(REAL*)fasp_mem_calloc(nnz,sizeof(REAL));
    }
    else {
        A->val = NULL;
    }
    
    A->row=m; A->col=n; A->nnz=nnz;
    
    return;    
}

/**
 * \fn void fasp_dcsr_free (dCSRmat *A)
 *
 * \brief Free CSR sparse matrix data memeory space
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \author Chensong Zhang
 * \date   2010/04/06 
 */
void fasp_dcsr_free (dCSRmat *A)
{    
    if (A==NULL) return;
    
    fasp_mem_free(A->IA);  A->IA  = NULL;
    fasp_mem_free(A->JA);  A->JA  = NULL;
    fasp_mem_free(A->val); A->val = NULL;
}

/**
 * \fn void fasp_icsr_free (iCSRmat *A)
 *
 * \brief Free CSR sparse matrix data memeory space
 *
 * \param A   Pointer to the iCSRmat matrix
 *
 * \author Chensong Zhang
 * \date   2010/04/06 
 */
void fasp_icsr_free (iCSRmat *A)
{    
    if (A==NULL) return;
    
    fasp_mem_free(A->IA);  A->IA  = NULL;
    fasp_mem_free(A->JA);  A->JA  = NULL;
    fasp_mem_free(A->val); A->val = NULL;
}

/**
 * \fn void fasp_dcsr_null (dCSRmat *A)
 *
 * \brief Initialize CSR sparse matrix
 * 
 * \param A   Pointer to the dCSRmat matrix
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
 */
void fasp_dcsr_null (dCSRmat *A)
{    
    A->row = A->col = A->nnz = 0;
    A->IA  = A->JA  = NULL;
    A->val = NULL;
}

/**
 * \fn dCSRmat fasp_dcsr_perm (dCSRmat *A, INT *P)
 *
 * \brief Apply permutation of A, i.e. Aperm=PAP' by the orders given in P
 *
 * \param A  Pointer to the oringinal dCSRmat matrix
 * \param P  Pointer to orders
 *
 * \return   The new ordered dCSRmat matrix if succeed, NULL if fail
 *
 * \author Shiquan Zhang
 * \date   03/10/2010
 * 
 * \note   P[i] = k means k-th row and column become i-th row and column!
 */
dCSRmat fasp_dcsr_perm (dCSRmat *A, 
                        INT *P)
{
    const INT n=A->row,nnz=A->nnz;
    const INT *ia=A->IA, *ja=A->JA;
    const REAL *Aval=A->val;
    INT i,j,k,jaj,i1,i2,start;
    
    dCSRmat Aperm = fasp_dcsr_create(n,n,nnz);
    
    // form the tranpose of P
    INT *Pt = (INT*)fasp_mem_calloc(n,sizeof(INT)); 
    for (i=0; i<n; ++i) Pt[P[i]] = i;
    
    // compute IA of P*A (row permutation)
    Aperm.IA[0] = 0;
    for (i=0; i<n; ++i) {
        k = P[i];
        Aperm.IA[i+1] = Aperm.IA[i]+(ia[k+1]-ia[k]);
    }
    
    // perform actual P*A
    for (i=0; i<n; ++i) {
        i1 = Aperm.IA[i]; i2 = Aperm.IA[i+1]-1; 
        k = P[i];
        start = ia[k];
        for (j=i1; j<=i2; ++j) {
            jaj = start+j-i1;
            Aperm.JA[j] = ja[jaj];
            Aperm.val[j] = Aval[jaj];
        }    
    }
    
    // perform P*A*P' (column permutation)
    for (k=0; k<nnz; ++k) {
        j = Aperm.JA[k];
        Aperm.JA[k] = Pt[j];
    }
    
    fasp_mem_free(Pt);
    
    return(Aperm);    
}

/**
 * \fn void fasp_dcsr_sort (dCSRmat *A)
 *
 * \brief Sort each row of A in ascending order w.r.t. column indices.
 *
 * \param A   Pointer to the dCSRmat matrix
 *
 * \author Shiquan Zhang
 * \date   06/10/2010
 */
void fasp_dcsr_sort (dCSRmat *A)
{
    const INT n=A->col;
    INT i,j,start,row_length;
    
    // temp memory for sorting rows of A
    INT *index, *ja; 
    REAL *a;
    
    index = (INT*)fasp_mem_calloc(n,sizeof(INT));     
    ja    = (INT*)fasp_mem_calloc(n,sizeof(INT));     
    a     = (REAL*)fasp_mem_calloc(n,sizeof(REAL));     
    
    for (i=0;i<n;++i) {
        start=A->IA[i];
        row_length=A->IA[i+1]-start;
    
        for (j=0;j<row_length;++j) index[j]=j;
    
        fasp_aux_iQuickSortIndex(&(A->JA[start]), 0, row_length-1, index);
    
        for (j=0;j<row_length;++j) {
            ja[j]=A->JA[start+index[j]];
            a[j]=A->val[start+index[j]];
        }
    
        for (j=0;j<row_length;++j) {
            A->JA[start+j]=ja[j];
            A->val[start+j]=a[j];
        }
    }
    
    // clean up memory
    fasp_mem_free(index);
    fasp_mem_free(ja);
    fasp_mem_free(a);
}

/**
 * \fn void fasp_dcsr_getdiag (INT n, dCSRmat *A, dvector *diag)
 *
 * \brief Get first n diagonal entries of a CSR matrix A
 *
 * \param n     Number of diag entries to get (if n=0, then get all diagonal entries)
 * \param A     Pointer to dCSRmat CSR matrix
 * \param diag  Pointer to the diagonal as a dvector
 *
 * \author Chensong Zhang
 * \date   05/20/2009
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */

void fasp_dcsr_getdiag (INT n, 
                        dCSRmat *A, 
                        dvector *diag) 
{
    INT i,k,j,ibegin,iend;    
    
    if (n==0) n=MIN(A->row,A->col);
	INT nthreads = 1, use_openmp = FALSE;

	if(FASP_USE_OPENMP && n > OPENMP_HOLDS){
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    fasp_dvec_alloc(n, diag);
    
    if (use_openmp) {
        INT mybegin,myend,myid;
#if FASP_USE_OPENMP
#pragma omp parallel for private(myid, mybegin, myend, i, ibegin, iend, k, j) 
#endif
        for (myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
            for (i=mybegin; i<myend; i++) {
                ibegin=A->IA[i]; iend=A->IA[i+1];
                for (k=ibegin;k<iend;++k) {
                    j=N2C(A->JA[N2C(k)]);
                    if ((j-i)==0) {
                        diag->val[i] = A->val[N2C(k)]; break;
                    } // end if
                } // end for k
            } // end for i
        }
    }
    else {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i]; iend=A->IA[i+1];
            for (k=ibegin;k<iend;++k) {
                j=N2C(A->JA[N2C(k)]);
                if ((j-i)==0) {
                    diag->val[i] = A->val[N2C(k)]; break;
                } // end if
            } // end for k
        } // end for i
    }
}

/**
 * \fn void fasp_dcsr_getcol (const INT n, dCSRmat *A, REAL *col)
 *
 * \brief Get the n-th column of a CSR matrix A
 *
 * \param n    Index of a column of A (0 <= n <= A.col-1) 
 * \param A    Pointer to dCSRmat CSR matrix
 * \param col  Pointer to the column
 *
 * \author Xiaozhe Hu
 * \date   11/07/2009
 */
void fasp_dcsr_getcol (const INT n, 
                       dCSRmat *A, 
                       REAL *col) 
{
    INT i,j, row_begin, row_end;
    INT nrow = A->row, ncol = A->col;
    INT status = SUCCESS;
    
    // check the column index n
    if (n < 0 || n >= ncol) {
        printf("### ERROR: the column index n is illegal!!\n");
        status = ERROR_DUMMY_VAR;
        goto FINISHED;
    }
    
    // get the column
    for (i=0; i<nrow; ++i) {
        // set the entry to zero
        col[i] = 0.0;
    
        row_begin = A->IA[i]; row_end = A->IA[i+1];
    
        for (j=row_begin; j<row_end; ++j) {
            if (A->JA[j] == n) {
                col[i] = A->val[j];
            }
        } // end for (j=row_begin; j<row_end; ++j)
    
    } // end for (i=0; i<nrow; ++i)
    
 FINISHED:
    fasp_chkerr(status,"fasp_dcsr_getcol");
}

/*!
 * \fn void fasp_dcsr_diagpref ( dCSRmat *A )
 *
 * \brief Re-order the column and data arrays of a CSR matrix, 
 *        so that the first entry in each row is the diagonal
 *
 * \param A   Pointer to the matrix to be re-ordered
 *
 * \author Zhiyang Zhou
 * \date   09/09/2010
 *
 * \note Reordering is done in place. 
 *
 * Modified by Chensong Zhang on Dec/21/2012
 */
void fasp_dcsr_diagpref (dCSRmat *A)
{
    const INT   num_rowsA = A -> row; 
    REAL      * A_data    = A -> val;
    INT       * A_i       = A -> IA;
    INT       * A_j       = A -> JA;
    
    // Local variable
    INT    i, j;
    INT    tempi, row_size;
    REAL   tempd;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_dcsr_diagpref ...... [Start]\n");
#endif
    
    for (i = 0; i < num_rowsA; i ++) {
        row_size = A_i[i+1] - A_i[i];
        
        // check whether the first entry is already diagonal        
        if (A_j[0] != i) { 
            
            for (j = 1; j < row_size; j ++) {
                if (A_j[j] == i) {
#if DEBUG_MODE
                    printf("### DEBUG: Switch entry_%d with entry_0\n", j);
#endif
                    
                    tempi  = A_j[0];
                    A_j[0] = A_j[j];
                    A_j[j] = tempi;
                    
                    tempd     = A_data[0];
                    A_data[0] = A_data[j];
                    A_data[j] = tempd;
                    
                    break;
                }
            }
            
            if (j == row_size) {
                printf("### ERROR: Diagonal entry in row %d is either missing or zero!\n", i);
                exit(ERROR_MISC); // Diagonal is zero    
            }
        }
        
        A_j    += row_size;
        A_data += row_size;
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_dcsr_diagpref ...... [Finish]\n");
#endif
}

/**
 * \fn SHORT fasp_dcsr_regdiag (dCSRmat *A, REAL value)
 *
 * \brief Regularize diagonal entries of a CSR sparse matrix
 *
 * \param A      Pointer to the dCSRmat matrix
 * \param value  Set a value on diag(A) which is too close to zero to "value"
 * 
 * \return       SUCCESS if no diagonal entry is close to zero, else ERROR
 *
 * \author Shiquan Zhang
 * \date   11/07/2009
 */
SHORT fasp_dcsr_regdiag (dCSRmat *A, 
                         REAL value)
{
    const INT    m = A->row;
    const INT * ia = A->IA, * ja = A->JA;
    REAL      * aj = A->val;
    
    // Local variables
    INT i,j,k,begin_row,end_row;
    SHORT status=RUN_FAIL;
    
    for (i=0;i<m;++i) {
        begin_row=ia[i],end_row=ia[i+1];
        for (k=begin_row;k<end_row;++k) {
            j=ja[k];
            if (i==j) {
                if (aj[k] < 0.0) goto FINISHED;
                else if (aj[k] < SMALLREAL) aj[k]=value;    
            }
        } // end for k
    } // end for i
    
    status = SUCCESS;
    
 FINISHED:    
    return status;
}

/**
 * \fn void fasp_dcsr_cp (dCSRmat *A, dCSRmat *B)
 *
 * \brief copy a dCSRmat to a new one B=A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param B   Pointer to the dCSRmat matrix
 *
 * \author Chensong Zhang
 * \date   04/06/2010  
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue  
 * \date   05/23/2012    
 */

void fasp_dcsr_cp (dCSRmat *A, 
                   dCSRmat *B) 
{
    B->row=A->row;
    B->col=A->col;
    B->nnz=A->nnz;
    
    fasp_iarray_cp (A->row+1, A->IA, B->IA);
    fasp_iarray_cp (A->nnz, A->JA, B->JA);
    fasp_array_cp  (A->nnz, A->val, B->val);
}

/**
 * \fn void fasp_icsr_trans (iCSRmat *A, iCSRmat *AT)
 *
 * \brief Find transpose of iCSRmat matrix A
 *
 * \param A   Pointer to the iCSRmat matrix A
 * \param AT  Pointer to the iCSRmat matrix A'
 *
 * \return    The transpose of iCSRmat matrix A
 *
 * \author Chensong Zhang
 * \date   04/06/2010  
 *
 *  Modified by Chunsheng Feng, Zheng Li
 * \date   06/20/2012   
 */
void fasp_icsr_trans (iCSRmat *A, 
                      iCSRmat *AT)
{
    const INT n=A->row, m=A->col, nnz=A->nnz, m1=m-1;
    
    // Local variables
    INT i,j,k,p;
    INT ibegin, iend;
    
#if DEBUG_MODE
    printf("### DEBUG: m=%d, n=%d, nnz=%d\n",m,n,nnz);
#endif
    
    AT->row=m; AT->col=n; AT->nnz=nnz;
    
    AT->IA=(INT*)fasp_mem_calloc(m+1,sizeof(INT));
    
    AT->JA=(INT*)fasp_mem_calloc(nnz,sizeof(INT));
    
    if (A->val) { 
        AT->val=(INT*)fasp_mem_calloc(nnz,sizeof(INT)); 
    }
    else { 
        AT->val=NULL; 
    }    
    
    // first pass: find the Number of nonzeros in the first m-1 columns of A 
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1
	
    //for (i=0;i<m;++i) AT->IA[i] = 0;   //Here is a Bug.
    for (i=0; i<=m; ++i) AT->IA[i] = 0;  //Chunsheng Feng ,Zheng Li, June/20/2012
    
    for (j=0;j<nnz;++j) {
        i=N2C(A->JA[j]); // column Number of A = row Number of A'
        if (i<m1) AT->IA[i+2]++;
    }
    
    for (i=2;i<=m;++i) AT->IA[i]+=AT->IA[i-1];
    
    // second pass: form A'
    if (A->val != NULL) {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i], iend=A->IA[i+1];    
            for(p=ibegin;p<iend;p++) {
                j=A->JA[N2C(p)]+1;
                k=AT->IA[N2C(j)];
                AT->JA[N2C(k)]=C2N(i);
                AT->val[N2C(k)]=A->val[N2C(p)];
                AT->IA[j]=k+1;
            } // end for p
        } // end for i
    }
    else {
        for (i=0;i<n;++i) {
            ibegin=A->IA[i], iend=A->IA[i+1];    
            for(p=ibegin;p<iend;p++) {
                j=A->JA[N2C(p)]+1;
                k=AT->IA[N2C(j)];
                AT->JA[N2C(k)]=C2N(i);
                AT->IA[j]=k+1;
            } // end for p
        } // end for i    
    } // end if
}    

/**
 * \fn void fasp_dcsr_trans (dCSRmat *A, dCSRmat *AT)
 *
 * \brief Find tranpose of dCSRmat matrix A
 *
 * \param A   Pointer to the dCSRmat matrix
 * \param AT  Pointer to the transpose of dCSRmat matrix A (output)
 *
 * \author Chensong Zhang
 * \date   04/06/2010   
 *
 *  Modified by Chunsheng Feng, Zheng Li
 * \date   06/20/2012   
 */
INT fasp_dcsr_trans (dCSRmat *A, 
                     dCSRmat *AT)
{
    const INT n=A->row, m=A->col, nnz=A->nnz;    
    
    // Local variables
    INT i,j,k,p;
    
    AT->row=m;
    AT->col=n;
    AT->nnz=nnz;
    
    AT->IA=(INT*)fasp_mem_calloc(m+1,sizeof(INT));
    
    AT->JA=(INT*)fasp_mem_calloc(nnz,sizeof(INT));
    
    if (A->val) { 
        AT->val=(REAL*)fasp_mem_calloc(nnz,sizeof(REAL)); 
    
    }
    else { AT->val=NULL; }
    
    // first pass: find the Number of nonzeros in the first m-1 columns of A 
    // Note: these Numbers are stored in the array AT.IA from 1 to m-1
	
    //for (i=0;i<m;++i) AT->IA[i] = 0;   //Here is a Bug.
    for (i=0; i<=m; ++i) AT->IA[i] = 0;  //Chunsheng Feng ,Zheng Li, June/20/2012
    
    for (j=0;j<nnz;++j) {
        i=N2C(A->JA[j]); // column Number of A = row Number of A'
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
                AT->val[N2C(k)]=A->val[N2C(p)];
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
    
    return SUCCESS;
}    

/* 
 * \fn void fasp_dcsr_transpose (INT *row[2], INT *col[2], REAL *val[2], INT *nn, INT *tniz)
 *
 * \brief Transpose of an CSR matrix.
 *
 * \param row     Pointers of the rows of the matrix and its transpose
 * \param col     Pointers of the columns of the matrix and its transpose
 * \param val     Pointers to the values of the matrix and its transpose
 * \param nn      Number of rows and columns of the matrix
 * \param tniz    Number of the nonzeros in the matrices A and A'
 *
 * \author Shuo Zhang
 * \date   07/06/2009   
 *
 * \note This subroutine transpose in CSR format IN ORDER. 
 */
void fasp_dcsr_transpose (INT *row[2], 
                          INT *col[2], 
                          REAL *val[2], 
                          INT *nn, 
                          INT *tniz)
{
    const INT nca=nn[1]; // Number of columns
    
    INT *izc    = (INT *)fasp_mem_calloc(nn[1],sizeof(INT));    
    INT *izcaux = (INT *)fasp_mem_calloc(nn[1],sizeof(INT));
    
    // Local variables
    INT i,m,itmp;
    
    // first pass: to set order right    
    for (i=0;i<tniz[0];++i) izc[N2C(col[0][i])]++;
    
    izcaux[0]=0;
    for (i=1;i<nca;++i) izcaux[i]=izcaux[i-1]+izc[i-1];
    
    // second pass: form transpose
    memset(izc,0,nca*sizeof(INT));
    
    for (i=0;i<tniz[0];++i) {
        m=col[0][i]; 
        itmp=izcaux[m]+izc[m];
        row[1][itmp]=m;
        col[1][itmp]=row[0][i];
        val[1][itmp]=val[0][i];
        izc[m]++;
    }
    
    fasp_mem_free(izc);
    fasp_mem_free(izcaux);
}

/**
 * \fn void fasp_dcsr_compress (dCSRmat *A, dCSRmat *B, REAL dtol)
 *
 * \brief Compress a CSR matrix A and store in CSR matrix B by
 *        dropping small entries abs(aij)<=dtol
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param B     Pointer to dCSRmat CSR matrix
 * \param dtol  Drop tolerance 
 *
 * \author Shiquan Zhang
 * \date   03/10/2010
 */
void fasp_dcsr_compress (dCSRmat *A, 
                         dCSRmat *B, 
                         REAL dtol)
{
    INT i, j, k;
    INT ibegin,iend1;    
    
    INT *index=(INT*)fasp_mem_calloc(A->nnz,sizeof(INT));
    
    B->row=A->row; B->col=A->col;
    
    B->IA=(INT*)fasp_mem_calloc(A->row+1,sizeof(INT));
    
    B->IA[0]=A->IA[0];
    
    // first pass: determine the size of B
    k=0;
    for (i=0;i<A->row;++i) {
        ibegin=A->IA[i]; iend1=A->IA[i+1];    
        for (j=ibegin;j<iend1;++j)
            if (ABS(A->val[N2C(j)])>dtol) {
                index[k]=N2C(j);
                ++k;
            } /* end of j */
        B->IA[i+1]=C2N(k);
    } /* end of i */
    B->nnz=k;
    
    B->JA=(INT*)fasp_mem_calloc(B->nnz,sizeof(INT));    
    B->val=(REAL*)fasp_mem_calloc(B->nnz,sizeof(REAL));
    
    // second pass: generate the index and element to B
    for (i=0;i<B->nnz;++i) {
        B->JA[i]=A->JA[index[i]];
        B->val[i]=A->val[index[i]]; 
    }  
    
    fasp_mem_free(index);
}

/**
 * \fn SHORT fasp_dcsr_compress_inplace (dCSRmat *A, REAL dtol)
 *
 * \brief Compress a CSR matrix A IN PLACE by
 *        dropping small entries abs(aij)<=dtol
 *
 * \param A     Pointer to dCSRmat CSR matrix
 * \param dtol  Drop tolerance
 *
 * \author Xiaozhe Hu
 * \date   12/25/2010
 *
 * \note This routine can be modified for filtering.
 */
SHORT fasp_dcsr_compress_inplace (dCSRmat *A, 
                                  REAL dtol)
{
    const INT row=A->row;
    const INT nnz=A->nnz;
    
    INT i, j, k;
    INT ibegin,iend=A->IA[0];    
    SHORT status = SUCCESS;
    
    k=0;
    for (i=0;i<row;++i) {
        ibegin=iend; iend=A->IA[i+1];    
        for (j=ibegin;j<iend;++j)
            if (ABS(A->val[N2C(j)])>dtol) {
                A->JA[N2C(k)] = A->JA[N2C(j)];
                A->val[N2C(k)] = A->val[N2C(j)];
                ++k;
            } /* end of j */
        A->IA[i+1]=C2N(k);
    } /* end of i */
    
    if (k<=nnz) 
        {
            A->nnz=k;
            A->JA = (INT*)fasp_mem_realloc(A->JA, k*sizeof(INT));
            A->val = (REAL *)fasp_mem_realloc(A->val, k*sizeof(REAL));
        }
    else 
        {
            status = RUN_FAIL;
            printf("### ERROR: Size of compressed matrix is larger than the original!!\n");
        }
    
    return (status);
}

/**
 * \fn void fasp_dcsr_shift (dCSRmat *A, INT offset)
 *
 * \brief Reindex a REAL matrix in CSR format to make the index starting from 0 or 1
 *
 * \param A         Pointer to CSR matrix
 * \param  offset   Size of offset (1 or -1)
 *
 * \author Chensong Zhang
 * \date   04/06/2010  
 */
void fasp_dcsr_shift (dCSRmat *A, 
                      INT offset)
{
    const INT nnz=A->nnz;
    const INT n=A->row;
    INT i, *ai=A->IA, *aj=A->JA;
    
    if (offset == 0) offset = ISTART;
    
    for (i=0; i<=n; ++i) ai[i]+=offset;
    
    for (i=0;i<nnz;++i) aj[i]+=offset;
}

/**
 * \fn void fasp_dcsr_symdiagscale (dCSRmat *A, dvector *diag)
 *
 * \brief Symmetric diagonal scaling D^{-1/2}AD^{-1/2}
 *
 * \param A     Pointer to the dCSRmat matrix
 * \param diag  Pointer to the diagonal entries
 *
 * \author Xiaozhe Hu
 * \date   01/31/2011
 */
void fasp_dcsr_symdiagscale (dCSRmat *A, 
                             dvector *diag)
{
    // information about matrix A
    const INT n = A->row;
    INT *IA = A->IA;
    INT *JA = A->JA;
    REAL *val = A->val;
    
    // local variables
    unsigned INT i, j, k, row_start, row_end;
    
    if (diag->row != n) {
        printf("### ERROR: Size of diag = %d and size of matrix = %d mismatch!!", 
               diag->row, n);
        exit(ERROR_MISC);
    }
    
    // work space
    REAL *work = (REAL *)fasp_mem_calloc(n, sizeof(REAL));
    
    // square root of diagonal entries
    for (i=0; i<n; i++) work[i] = sqrt(diag->val[i]);
    
    // main loop
    for (i=0; i<n; i++) {
        row_start = IA[i]; row_end = IA[i+1];
        for (j=row_start; j<row_end; j++) {
            k = JA[j];
            val[j] = val[j]/(work[i]*work[k]);
        }
    }
    
    // free work space
    if (work) fasp_mem_free(work);
}

/**
 * \fn dCSRmat fasp_dcsr_sympat(dCSRmat *A)
 * \brief Symmetrize the parttarn of a dCSRmat matrix
 *
 * \param *A      pointer to the dCSRmat matrix
 *
 * \return symmetrized the dCSRmat matrix
 *
 * \author Xiaozhe Hu
 * \date 03/21/2011
 */
dCSRmat fasp_dcsr_sympat(dCSRmat *A)
{
	//local variable
	dCSRmat AT;
	
	//return variable
	dCSRmat SA;
	
	// get the transpose of A
	fasp_dcsr_trans (A,  &AT);
	
	// get symmetrized A
	fasp_blas_dcsr_add (A, 1.0, &AT, 0.0, &SA);
	
	// clean 
	fasp_dcsr_free(&AT);
	
	// return 
	return SA;
	
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
