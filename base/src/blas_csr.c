/*! \file blas_csr.c
 *  \brief BLAS operations for sparse matrices in CSR format.
 *
 *  \note Sparse functions usually contain three runs. The three runs are all the 
 *        same but thy serve different purpose. 
 * 
 *  Example: If you do c=a+b: 
 *    - first do a dry run to find the number of non-zeroes in the result and form ic; 
 *    - allocate space (memory) for jc and form this one; 
 *    - if you only care about a "boolean" result of the addition, you stop here;
 *    - you call another routine, which uses ic and jc to perform the addition. 
 *
 */
#include <omp.h>
#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_dcsr_add (dCSRmat *A, const REAL alpha, dCSRmat *B,  
 *                              const REAL beta, dCSRmat *C)
 *
 * \brief compute C = alpha*A + beta*B in CSR format 
 *
 * \param A      Pointer to CSR matrix
 * \param alpha  REAL factor alpha 
 * \param B      Pointer to CSR matrix
 * \param beta   REAL factor beta 
 * \param C      Pointer to CSR matrix
 *
 * \return       SUCCESS if succees, RUN_FAIL if not
 *
 * \author Xiaozhe Hu
 * \date   11/07/2009
 *
 * Modified by Chunsheng Feng, Zheng Li on 06/29/2012
 */
INT fasp_blas_dcsr_add (dCSRmat *A, 
                        const REAL alpha, 
                        dCSRmat *B, 
                        const REAL beta, 
                        dCSRmat *C)
{
    INT i,j,k,l;
    INT count=0, added, countrow;
    INT status = SUCCESS;

    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP  
	if ( A->nnz > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

    if (A->row != B->row || A->col != B->col) {
#if DEBUG_MODE
        printf("### DEBUG: The dim of two matrices does not match --fasp_blas_dcsr_add!\n");
#endif
        status = ERROR_DATA_STRUCTURE;
        goto FINISHED;
    }
    
    if (A == NULL && B == NULL) {
        C->row=0; C->col=0; C->nnz=0;
        status=SUCCESS; goto FINISHED;
    }
    
    if (A->nnz == 0 && B->nnz == 0) {
        C->row=A->row; C->col=A->col; C->nnz=A->nnz;
        status=SUCCESS; goto FINISHED;
    }
    
    // empty matrix A
    if (A->nnz == 0 || A == NULL) {
        fasp_dcsr_alloc(B->row,B->col,B->nnz,C);
        memcpy(C->IA,B->IA,(B->row+1)*sizeof(INT));
        memcpy(C->JA,B->JA,(B->nnz)*sizeof(INT));

        if (use_openmp) {
#ifdef _OPENMP 
            INT mybegin, myend, myid;
#pragma omp parallel private(myid, mybegin, myend, i)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, A->nnz, &mybegin, &myend);
                for (i=mybegin;i<myend;++i) C->val[i]=B->val[i]*beta;
            }
#endif
        }
        else {
            for (i=0;i<A->nnz;++i) C->val[i]=B->val[i]*beta;
        }

        status=SUCCESS; 
        goto FINISHED;
    }
    
    // empty matrix B
    if (B->nnz == 0 || B == NULL) {
        fasp_dcsr_alloc(A->row,A->col,A->nnz,C);
        memcpy(C->IA,A->IA,(A->row+1)*sizeof(INT));
        memcpy(C->JA,A->JA,(A->nnz)*sizeof(INT));

        if (use_openmp) {
#ifdef _OPENMP 
            INT mybegin, myend, myid;
#pragma omp parallel private(myid, mybegin, myend, i)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, A->nnz, &mybegin, &myend);
                for (i=mybegin;i<myend;++i) C->val[i]=A->val[i]*alpha;
            }
#endif
        }
        else {
            for (i=0;i<A->nnz;++i) C->val[i]=A->val[i]*alpha;
        }

        status=SUCCESS; goto FINISHED;
    }
    
    C->row=A->row; C->col=A->col;
    
    C->IA=(INT*)fasp_mem_calloc(C->row+1,sizeof(INT));
    
    // allocate work space for C->JA and C->val
    C->JA=(INT *)fasp_mem_calloc(A->nnz+B->nnz,sizeof(INT));
    
    C->val=(REAL *)fasp_mem_calloc(A->nnz+B->nnz,sizeof(REAL));

    // initial C->IA 
    memset(C->IA, 0, sizeof(INT)*(C->row+1));
    
    for (i=0; i<A->row; ++i) {
        countrow = 0;
        for (j=A->IA[i]; j<A->IA[i+1]; ++j) {
            C->val[count] = alpha * A->val[N2C(j)];
            C->JA[count] = A->JA[N2C(j)];
            C->IA[i+1]++;
            count++;
            countrow++;
        } // end for j
    
        for (k=B->IA[i]; k<B->IA[i+1]; ++k) {
            added = 0;
    
            for (l=C->IA[i]; l<C->IA[i]+countrow+1; l++) {
                if (B->JA[N2C(k)] == C->JA[N2C(l)]) {
                    C->val[N2C(l)] = C->val[N2C(l)] + beta * B->val[N2C(k)];
                    added = 1;
                    break;
                }
            } // end for l
    
            if (added == 0) {
                C->val[count] = beta * B->val[N2C(k)];
                C->JA[count] = B->JA[N2C(k)];
                C->IA[i+1]++;
                count++;
            }    
    
        } // end for k
    
        C->IA[i+1] += C->IA[i];
    
    }
    
    C->nnz = count;
    C->JA  = (INT *)fasp_mem_realloc(C->JA, (count)*sizeof(INT));
    C->val = (REAL *)fasp_mem_realloc(C->val, (count)*sizeof(REAL));
    
 FINISHED:
    return status;
}

/**
 * \fn void fasp_blas_dcsr_axm (dCSRmat *A, const REAL alpha)
 *
 * \brief Multiply a sparse matrix A in CSR format by a scalar alpha.
 *
 * \param A      Pointer to CSR matrix A
 * \param alpha  REAL factor alpha 
 *
 * \author Chensong Zhang
 * \date   07/01/2009
 *
 * Modified by Chunsheng Feng, Zheng Li on 06/29/2012
 */
void fasp_blas_dcsr_axm (dCSRmat *A, 
                         const REAL alpha)
{
    const INT nnz=A->nnz;
    
    // A direct calculation can be written as:
    fasp_blas_array_ax(nnz, alpha, A->val);
}

/**
 * \fn void fasp_blas_dcsr_mxv (dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \author Chensong Zhang
 * \date   07/01/2009
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/26/2012    
 *
 */
void fasp_blas_dcsr_mxv (dCSRmat *A, 
                         REAL *x, 
                         REAL *y) 
{
    const INT m=A->row;
    const INT *ia=A->IA, *ja=A->JA;
    const REAL *aj=A->val;
    INT i, k, begin_row, end_row, nnz_num_row;
    register REAL temp;
    
    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP  
	if ( m > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

   if (use_openmp) {
        INT myid, mybegin, myend;

#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, nnz_num_row, k) 
#endif
        for (myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, m, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) {
                temp=0.0; 
                begin_row = ia[i];
                end_row = ia[i+1];
                nnz_num_row = end_row - begin_row;
                switch(nnz_num_row) {
                case 3:
                    k = begin_row;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    break;
                case 4:
                    k = begin_row;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    break;
                case 5:
                    k = begin_row;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    break;
                case 6:
                    k = begin_row;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    break;
                case 7:
                    k = begin_row;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    k ++;
                    temp += aj[k]*x[ja[k]];
                    break;
                default:
                    for (k=begin_row; k<end_row; ++k) {
                        temp += aj[k]*x[ja[k]];
                    }
                    break;
                }
                y[i]=temp;
            }
        }
    }

    else {
        for (i=0;i<m;++i) {
            temp=0.0;
            begin_row=ia[i];
            end_row=ia[i+1];
            nnz_num_row = end_row - begin_row;
            switch(nnz_num_row) {
            case 3:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 4:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 5:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 6:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            case 7:
                k=begin_row;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                k ++;
                temp+=aj[k]*x[ja[k]];
                break;
            default:
                for (k=begin_row; k<end_row; ++k) {
                    temp+=aj[k]*x[ja[k]];
                }
                break;
            }
            y[i]=temp;
        }
    }
}


/**
 * \fn void fasp_blas_dcsr_mxv_agg (dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = A*x, where the entries of A are all ones.
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \author Xiaozhe Hu
 * \date   02/22/2011
 *
 * Modified by Chunsheng Feng, Zheng Li
 * \date   08/28/2012
 */

void fasp_blas_dcsr_mxv_agg (dCSRmat *A, 
                             REAL *x, 
                             REAL *y)
{
    const INT  m  = A->row;
    const INT *ia = A->IA, *ja = A->JA;
    INT i, k, begin_row, end_row;    
    register REAL temp;
    
#ifdef _OPENMP
    if (m > OPENMP_HOLDS) {
        INT myid, mybegin, myend;
	INT nthreads = FASP_GET_NUM_THREADS();
#pragma omp parallel for private(myid, i, mybegin, myend, temp, begin_row, end_row, k)
	for (myid=0; myid<nthreads; myid++) {
	    FASP_GET_START_END(myid, nthreads, m, &mybegin, &myend);
	    for (i=mybegin; i<myend; i++) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                y[i]=temp;
            }
        }
    }
    else {
#endif
        for (i=0;i<m;++i) {
            temp=0.0; 
            begin_row=ia[i]; end_row=ia[i+1]; 
            for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
            y[i]=temp;
        }
#ifdef _OPENMP
    }
#endif

}

/**
 * \fn void fasp_blas_dcsr_aAxpy (const REAL alpha, dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 *
 * \param alpha  REAL factor alpha 
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Chensong Zhang
 * \date   07/01/2009
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/26/2012    
 */

void fasp_blas_dcsr_aAxpy (const REAL alpha, 
                           dCSRmat *A, 
                           REAL *x,
                           REAL *y)
{
    const INT  m  = A->row;
    const INT *ia = A->IA, *ja = A->JA;
    const REAL *aj = A->val;
    INT i, k, begin_row, end_row;
    register REAL temp;

    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP  
    if ( m > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif

    if ( alpha == 1.0 ) {
        if (use_openmp) {
            INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k)
#endif
            for (myid = 0; myid < nthreads; myid++ ) {
                FASP_GET_START_END(myid, nthreads, m, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) {
                    temp=0.0; 
                    begin_row=ia[i]; end_row=ia[i+1]; 
                    for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
                        y[i]+=temp;
                }
            }
        }
        else {
            for (i=0;i<m;++i) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
                y[i]+=temp;    
            }
        }
    }
    
    else if ( alpha == -1.0 ) {
        if (use_openmp) {
            INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k) 
#endif
            for (myid = 0; myid < nthreads; myid++ ) {
                FASP_GET_START_END(myid, nthreads, m, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) {
                    temp=0.0; 
                    begin_row=ia[i]; end_row=ia[i+1]; 
                    for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
                    y[i]-=temp;
                }
            }
        }
        else {
            for (i=0;i<m;++i) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
                y[i]-=temp;
            }
        }
    }
    
    else {
        if (use_openmp) {
            INT myid, mybegin, myend;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k) 
#endif
            for (myid = 0; myid < nthreads; myid++ ) {
                FASP_GET_START_END(myid, nthreads, m, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) {
                    temp=0.0;
                    begin_row=ia[i]; end_row=ia[i+1];
                    for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
                    y[i]+=temp*alpha;
                }
            }
        }
        else {
            for (i=0;i<m;++i) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
                y[i]+=temp*alpha;
            }
        }
    }
}

/**
 * \fn void fasp_blas_dcsr_aAxpy_agg (const REAL alpha, dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief Matrix-vector multiplication y = alpha*A*x + y, where the entries of A are all ones
 *
 * \param alpha  REAL factor alpha 
 * \param A      Pointer to dCSRmat matrix A
 * \param x      Pointer to array x
 * \param y      Pointer to array y
 *
 * \author Xiaozhe Hu
 * \date   02/22/2011
 *
 * Modified by Chunsheng Feng, Zheng Li
 * \date   08/28/2012
 */
void fasp_blas_dcsr_aAxpy_agg (const REAL alpha, 
                               dCSRmat *A, 
                               REAL *x, 
                               REAL *y)
{
    const INT  m  = A->row;
    const INT *ia = A->IA, *ja = A->JA;
    
    INT i, k, begin_row, end_row;    
    register REAL temp;
    
    if ( alpha == 1.0 ) {
#ifdef _OPENMP
        if (m > OPENMP_HOLDS) {
            INT myid, mybegin, myend;
	    INT nthreads = FASP_GET_NUM_THREADS();
#pragma omp parallel for private(myid, i, mybegin, myend, temp, begin_row, end_row, k)
            for (myid=0; myid<nthreads; myid++) {
                FASP_GET_START_END(myid, nthreads, m, &mybegin， &myend);
                for (i=mybegin;i<myend;++i) {
                    temp=0.0; 
                    begin_row=ia[i]; end_row=ia[i+1]; 
                    for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                    y[i]+=temp;    
                }
            }
        }
	else {
#endif
            for (i=0;i<m;++i) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                y[i]+=temp;    
            }
#ifdef _OPENMP
	}
#endif
    }
    else if ( alpha == -1.0 ) {
#ifdef _OPENMP
        if (m > OPENMP_HOLDS) {
            INT myid, mybegin, myend;
	    INT nthreads = FASP_GET_NUM_THREADS();
#pragma omp parallel for private(myid, i, mybegin, myend, temp, begin_row, end_row, k)
            for (myid=0; myid<nthreads; myid++) {
                FASP_GET_START_END(myid, nthreads, m, &mybegin， &myend);
                for (i=mybegin;i<myend;++i) {
                    temp=0.0; 
                    begin_row=ia[i]; end_row=ia[i+1]; 
                    for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                    y[i]-=temp;    
                }
            }
        }
	else {
#endif
            for (i=0;i<m;++i) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                y[i]-=temp;    
            }
#ifdef _OPENMP
	}
#endif
    }
    
    else {
#ifdef _OPENMP
        if (m > OPENMP_HOLDS) {
            INT myid, mybegin, myend;
	    INT nthreads = FASP_GET_NUM_THREADS();
#pragma omp parallel for private(myid, i, mybegin, myend, temp, begin_row, end_row, k)
            for (myid=0; myid<nthreads; myid++) {
                FASP_GET_START_END(myid, nthreads, m, &mybegin， &myend);
                for (i=mybegin;i<myend;++i) {
                    temp=0.0; 
                    begin_row=ia[i]; end_row=ia[i+1]; 
                    for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                    y[i]+=temp*alpha;    
                }
            }
        }
	else {
#endif
            for (i=0;i<m;++i) {
                temp=0.0; 
                begin_row=ia[i]; end_row=ia[i+1]; 
                for (k=begin_row; k<end_row; ++k) temp+=x[ja[k]];
                y[i]+=temp*alpha;    
            }
#ifdef _OPENMP
	}
#endif
    }
    
}

/**
 * \fn REAL fasp_blas_dcsr_vmv (dCSRmat *A, REAL *x, REAL *y) 
 *
 * \brief vector-Matrix-vector multiplication alpha = y'*A*x
 *
 * \param A   Pointer to dCSRmat matrix A
 * \param x   Pointer to array x
 * \param y   Pointer to array y
 *
 * \author Chensong Zhang
 * \date   07/01/2009
 */
REAL fasp_blas_dcsr_vmv (dCSRmat *A, 
                         REAL *x, 
                         REAL *y)
{
    register REAL value=0.0;
    const INT m=A->row;
    const INT *ia=A->IA, *ja=A->JA;
    const REAL *aj=A->val;
    INT i, k, begin_row, end_row;
    register REAL temp;

    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP  
    if ( m > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel for reduction(+:value) private(i,temp,begin_row,end_row,k) 
#endif
        for (i=0;i<m;++i) {
            temp=0.0;
            begin_row=ia[i]; end_row=ia[i+1];
            for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
            value+=y[i]*temp;
        }
    }
    else {
        for (i=0;i<m;++i) {
            temp=0.0;
            begin_row=ia[i]; end_row=ia[i+1];
            for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
            value+=y[i]*temp;
        }
    }
    return value;
}

/**
 * \fn void fasp_blas_dcsr_mxm (dCSRmat *A, dCSRmat *B, dCSRmat *C)
 *
 * \brief Sparse matrix multiplication C=A*B
 *
 * \param A   Pointer to the dCSRmat matrix A
 * \param B   Pointer to the dCSRmat matrix B
 * \param C   Pointer to dCSRmat matrix equal to A*B
 *
 * \author Xiaozhe Hu
 * \date   11/07/2009
 *
 * \note This fct will be replaced! --Chensong
 */
void fasp_blas_dcsr_mxm (dCSRmat *A, 
                         dCSRmat *B, 
                         dCSRmat *C)
{    
    INT i,j,k,l,count;
    INT *JD = (INT *)fasp_mem_calloc(B->col,sizeof(INT));
    
    C->row=A->row;
    C->col=B->col;
    
    C->val = NULL;
    C->JA  = NULL;    
    C->IA  = (INT*)fasp_mem_calloc(C->row+1,sizeof(INT));
    
    for (i=0;i<B->col;++i) JD[i]=-1;
    
    // step 1: Find first the structure IA of C
    for(i=0;i<C->row;++i) {
        count=0;
    
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<count;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
    
                if (l==count) {
                    JD[count]=B->JA[j];
                    count++;
                }
            }
        }
        C->IA[i+1]=count;
        for (j=0;j<count;++j) {
            JD[j]=-1;
        }    
    }    
    
    for (i=0;i<C->row;++i) C->IA[i+1]+=C->IA[i];
    
    // step 2: Find the structure JA of C
    INT countJD;
    
    C->JA=(INT*)fasp_mem_calloc(C->IA[C->row],sizeof(INT));
    
    for (i=0;i<C->row;++i) {
        countJD=0;
        count=C->IA[i];
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            for (j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];++j) {
                for (l=0;l<countJD;l++) {
                    if (JD[l]==B->JA[j]) break;
                }
    
                if (l==countJD) {
                    C->JA[count]=B->JA[j];
                    JD[countJD]=B->JA[j];
                    count++;
                    countJD++;
                }
            }
        }
    
        for (j=0;j<countJD;++j) JD[j]=-1;
    }    
    
    fasp_mem_free(JD);
    
    // step 3: Find the structure A of C
    C->val=(REAL*)fasp_mem_calloc(C->IA[C->row],sizeof(REAL));
    
    for (i=0;i<C->row;++i) {
        for (j=C->IA[i];j<C->IA[i+1];++j) {
            C->val[j]=0;
            for (k=A->IA[i];k<A->IA[i+1];++k) {
                for (l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++) {
                    if (B->JA[l]==C->JA[j]) {
                        C->val[j]+=A->val[k]*B->val[l];
                    } // end if
                } // end for l
            } // end for k
        } // end for j
    }    // end for i
    
    C->nnz = C->IA[C->row]-C->IA[0];
    
}

/**
 * \fn void fasp_blas_dcsr_rap (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param R   Pointer to the dCSRmat matrix R
 * \param A   Pointer to the dCSRmat matrix A
 * \param P   Pointer to the dCSRmat matrix P
 * \param RAP Pointer to dCSRmat matrix equal to R*A*P
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   05/10/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/26/2012 
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *       Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void fasp_blas_dcsr_rap( dCSRmat  *R,
                         dCSRmat  *A,
                         dCSRmat  *P,
                         dCSRmat  *RAP)
{
    INT n_coarse = R->row;
    INT *R_i = R->IA;
    INT *R_j = R->JA;
    REAL *R_data = R->val;
    
    INT n_fine = A->row;
    INT *A_i = A->IA;
    INT *A_j = A->JA;
    REAL *A_data = A->val;
    
    INT *P_i = P->IA;
    INT *P_j = P->JA;
    REAL *P_data = P->val;
    
    INT RAP_size;
    INT *RAP_i = NULL;
    INT *RAP_j = NULL;
    REAL *RAP_data = NULL;

#ifdef _OPENMP 
    INT *P_marker = NULL;
    INT *A_marker = NULL;
#endif

    INT *Ps_marker = NULL;
    INT *As_marker = NULL;
    INT *RAP_temp = NULL;
    INT *part_end = NULL;
    
    INT ic, i1, i2, i3, jj1, jj2, jj3;
    INT jj_counter, jj_row_begining;
    REAL r_entry, r_a_product, r_a_p_product;
    
    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP  
    if ( n_coarse > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;
    
    Ps_marker = (INT *)fasp_mem_calloc(total_calloc, sizeof(INT));
    As_marker = Ps_marker + coarse_mul_nthreads;
    RAP_temp = As_marker + fine_mul_nthreads;
    part_end = RAP_temp + coarse_add_nthreads;
    
    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
    RAP_i = (INT *)fasp_mem_calloc(n_coarse+1, sizeof(INT));
    
    fasp_iarray_set(minus_one_length, Ps_marker, -1);
    
    if (use_openmp) {
#ifdef _OPENMP 
        INT myid, mybegin, myend, Ctemp;
#pragma omp parallel private(myid, mybegin, myend, P_marker, A_marker, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            jj_counter = 0;
            for (ic = mybegin; ic < myend; ic ++) {
                P_marker[ic] = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter ++;
                
                for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic) {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3] = jj_counter;
                                    jj_counter ++;
                                }
                            }
                        }
                    }
                }
                
                RAP_temp[ic+myid] = jj_row_begining;
            }
            RAP_temp[myend+myid] = jj_counter;
            
            part_end[myid] = myend + myid + 1;
        }
        fasp_iarray_cp(part_end[0], RAP_temp, RAP_i);
        jj_counter = part_end[0];
        Ctemp = 0;
        for (i1 = 1; i1 < nthreads; i1 ++) {
            Ctemp += RAP_temp[part_end[i1-1]-1];
            for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++) {
                RAP_i[jj_counter] = RAP_temp[jj1] + Ctemp;
                jj_counter ++;
            }
        }
        RAP_size = RAP_i[n_coarse];
#endif
    }

    else {
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic ++) {
            Ps_marker[ic] = jj_counter;
            jj_row_begining = jj_counter;
            jj_counter ++;
            
            for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
                i1 = R_j[jj1];
                
                for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic) {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3] = jj_counter;
                                jj_counter ++;
                            }
                        }
                    }
                }
            }
            
            RAP_i[ic] = jj_row_begining;
        }
        
        RAP_i[n_coarse] = jj_counter;
        RAP_size = jj_counter;
    }
    //    printf("A_H NNZ = %d\n", RAP_size);
    RAP_j = (INT *)fasp_mem_calloc(RAP_size, sizeof(INT));
    RAP_data = (REAL *)fasp_mem_calloc(RAP_size, sizeof(REAL));
    
    fasp_iarray_set(minus_one_length, Ps_marker, -1);
    
    if (use_openmp) {
#ifdef _OPENMP 
        INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, P_marker, A_marker, jj_counter, ic, jj_row_begining, \
                             jj1, r_entry, i1, jj2, r_a_product, i2, jj3, r_a_p_product, i3)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n_coarse, &mybegin, &myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            jj_counter = RAP_i[mybegin];
            for (ic = mybegin; ic < myend; ic ++) {
                P_marker[ic] = jj_counter;
                jj_row_begining = jj_counter;
                RAP_j[jj_counter] = ic;
                RAP_data[jj_counter] = 0.0;
                jj_counter ++;
                for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
                    r_entry = R_data[jj1];
                    
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
                        r_a_product = r_entry * A_data[jj2];
                        
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic) {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
                                r_a_p_product = r_a_product * P_data[jj3];
                                
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining) {
                                    P_marker[i3] = jj_counter;
                                    RAP_data[jj_counter] = r_a_p_product;
                                    RAP_j[jj_counter] = i3;
                                    jj_counter ++;
                                }
                                else {
                                    RAP_data[P_marker[i3]] += r_a_p_product;
                                }
                            }
                        }
                        else {
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
                                i3 = P_j[jj3];
                                r_a_p_product = r_a_product * P_data[jj3];
                                RAP_data[P_marker[i3]] += r_a_p_product;
                            }
                        }
                    }
                }
            }
        }
#endif
    }

    else {
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic ++) {
            Ps_marker[ic] = jj_counter;
            jj_row_begining = jj_counter;
            RAP_j[jj_counter] = ic;
            RAP_data[jj_counter] = 0.0;
            jj_counter ++;
            
            for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++) {
                r_entry = R_data[jj1];
                
                i1 = R_j[jj1];
                for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++) {
                    r_a_product = r_entry * A_data[jj2];
                    
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic) {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
                            r_a_p_product = r_a_product * P_data[jj3];
                            
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining) {
                                Ps_marker[i3] = jj_counter;
                                RAP_data[jj_counter] = r_a_p_product;
                                RAP_j[jj_counter] = i3;
                                jj_counter ++;
                            }
                            else {
                                RAP_data[Ps_marker[i3]] += r_a_p_product;
                            }
                        }
                    }
                    else {
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++) {
                            i3 = P_j[jj3];
                            r_a_p_product = r_a_product * P_data[jj3];
                            RAP_data[Ps_marker[i3]] += r_a_p_product;
                        }
                    }
                }
            }
        }
    }
    
    RAP->row = n_coarse;
    RAP->col = n_coarse;
    RAP->nnz = RAP_size;
    RAP->IA = RAP_i;
    RAP->JA = RAP_j;
    RAP->val = RAP_data;
    
    fasp_mem_free(Ps_marker);
}

/**
 * \fn void fasp_blas_dcsr_rap_agg (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P, where the entries of R and P are all ones.
 *
 * \param R   Pointer to the dCSRmat matrix R 
 * \param A   Pointer to the dCSRmat matrix A
 * \param P   Pointer to the dCSRmat matrix P
 * \param B   Pointer to dCSRmat matrix equal to R*A*P
 *
 * \author Xiaozhe Hu
 * \date   02/21/2011
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *       Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void fasp_blas_dcsr_rap_agg (dCSRmat *R, 
                             dCSRmat *A, 
                             dCSRmat *P, 
                             dCSRmat *B)
{
    const INT row=R->row, col=P->col;
    INT nB=A->nnz;
    INT i,i1,j,jj,k,length;    
    INT begin_row,end_row,begin_rowA,end_rowA,begin_rowR,end_rowR;
    INT istart,iistart,count;
    
    REAL *aj=A->val, *acj;
    INT *ir=R->IA, *ia=A->IA, *ip=P->IA, *iac;
    INT *jr=R->JA, *ja=A->JA, *jp=P->JA, *jac;
    
    INT *index  = (INT *)fasp_mem_calloc(A->col,sizeof(INT));
    
    INT *iindex = (INT *)fasp_mem_calloc(col,sizeof(INT));    
    
    for (i=0; i<A->col; ++i) index[i] = -2;
    
    memcpy(iindex,index,col*sizeof(INT));
    
    jac = (INT*)fasp_mem_calloc(nB,sizeof(INT));    
    
    iac = (INT*)fasp_mem_calloc(row+1,sizeof(INT));    
    
    REAL *temp = (REAL*)fasp_mem_calloc(A->col,sizeof(REAL));
    
    iac[0] = 0;
    
    // First loop: form sparsity partern of R*A*P
    for (i=0; i < row; ++i) {    
        // reset istart and length at the begining of each loop
        istart = -1; length = 0; i1 = i+1;
    
        // go across the rows in R
        begin_rowR=ir[i]; end_rowR=ir[i1];
        for (jj=begin_rowR; jj<end_rowR; ++jj) {
            j = jr[N2C(jj)];
            // for each column in A
            begin_rowA=ia[j]; end_rowA=ia[j+1];    
            for (k=begin_rowA; k<end_rowA; ++k) {
                if (index[N2C(ja[N2C(k)])] == -2) {
                    index[N2C(ja[N2C(k)])] = istart;
                    istart = ja[N2C(k)];
                    ++length;
                }
            }
        }    
    
        // book-keeping [reseting length and setting iistart]
        count = length; iistart = -1; length = 0;
    
        // use each column that would have resulted from R*A
        for (j=0; j < count; ++j) {
            jj = istart;
            istart = index[istart];
            index[N2C(jj)] = -2;
    
            // go across the row of P
            begin_row=ip[jj]; end_row=ip[jj+1];
            for (k=begin_row; k<end_row; ++k) {
                // pull out the appropriate columns of P
                if (iindex[N2C(jp[N2C(k)])] == -2) {
                    iindex[N2C(jp[N2C(k)])] = iistart;
                    iistart = jp[N2C(k)];
                    ++length;
                }
            } // end for k
        } // end for j
    
        // set B->IA
        iac[i1]=iac[i]+length;
    
        if (iac[i1]>nB) {
            nB=nB*2;
            jac=(INT*)fasp_mem_realloc(jac, nB*sizeof(INT));
        }
    
        // put the correct columns of p into the column list of the products
        begin_row=iac[i]; end_row=iac[i1];
        for (j=begin_row; j<end_row; ++j) {
            // put the value in B->JA
            jac[N2C(j)] = iistart;
            // set istart to the next value
            iistart = iindex[N2C(iistart)];
            // set the iindex spot to 0
            iindex[N2C(jac[j])] = -2;
        } // end j
    
    } // end i: First loop
    
    jac=(INT*)fasp_mem_realloc(jac,(iac[row])*sizeof(INT));
    
    acj=(REAL*)fasp_mem_calloc(iac[row],sizeof(REAL));

    INT *BTindex=(INT*)fasp_mem_calloc(col,sizeof(INT));
    
    // Second loop: compute entries of R*A*P
    for (i=0; i<row; ++i) {
        i1 = i+1;
    
        // each col of B
        begin_row=iac[i]; end_row=iac[i1];    
        for (j=begin_row; j<end_row; ++j) {
            BTindex[N2C(jac[N2C(j)])]=j;
        }
    
        // reset istart and length at the begining of each loop
        istart = -1; length = 0;
    
        // go across the rows in R
        begin_rowR=ir[i]; end_rowR=ir[i1];    
        for ( jj=begin_rowR; jj<end_rowR; ++jj ) {
            j = jr[N2C(jj)];
    
            // for each column in A
            begin_rowA=ia[j]; end_rowA=ia[j+1];    
            for (k=begin_rowA; k<end_rowA; ++k) {
                if (index[N2C(ja[N2C(k)])] == -2) {
                    index[N2C(ja[N2C(k)])] = istart;
                    istart = ja[N2C(k)];
                    ++length;
                }
                temp[N2C(ja[N2C(k)])]+=aj[N2C(k)];
            }
        } 
    
        // book-keeping [reseting length and setting iistart]
        // use each column that would have resulted from R*A
        for (j=0; j<length; ++j) {
            jj = N2C(istart);
            istart = index[N2C(istart)];
            index[N2C(jj)] = -2;
    
            // go across the row of P
            begin_row=ip[jj]; end_row=ip[jj+1];    
            for (k=begin_row; k<end_row; ++k) {
                // pull out the appropriate columns of P
                acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj];
            }
            temp[jj]=0.0;
        }
    
    } // end for i: Second loop
    
    // setup coarse matrix B
    B->row=row; B->col=col;
    B->IA=iac; B->JA=jac; B->val=acj;    
    B->nnz=B->IA[B->row]-B->IA[0];
    
    fasp_mem_free(temp);
    fasp_mem_free(index);
    fasp_mem_free(iindex);
    fasp_mem_free(BTindex);
}

/**
 * \fn void fasp_blas_dcsr_ptap (dCSRmat *Pt, dCSRmat *A, dCSRmat *P, dCSRmat *Ac)
 *
 * \brief Triple sparse matrix multiplication B=P'*A*P
 *
 * \param Pt  Pointer to the restriction matrix
 * \param A   Pointer to the fine coefficient matrix
 * \param P   Pointer to the prolongation matrix
 * \param Ac  Pointer to the coarse coefficient matrix (output)
 *
 * \author Ludmil Zikatanov, Chensong Zhang
 * \date   05/10/2010
 *
 * \note Driver to compute triple matrix product P'*A*P using ltz CSR format. 
 *       In ltx format: ia[0]=1, ja[0] and a[0] are used as usual. When called 
 *       from Fortran, ia[0], ja[0] and a[0] will be just ia(1),ja(1),a(1).
 *       For the indices, 
 *           ia_ltz[k] = ia_usual[k]+1, 
 *           ja_ltz[k] = ja_usual[k]+1,
 *            a_ltz[k] =  a_usual[k].
 */
void fasp_blas_dcsr_ptap (dCSRmat *Pt,
                          dCSRmat *A, 
                          dCSRmat *P, 
                          dCSRmat *Ac)
{
    const INT nc=Pt->row, n=Pt->col, nnzP=P->nnz, nnzA=A->nnz;
    INT i, maxrpout;
    
    // shift A from usual to ltz format
    for (i=0;i<=n;++i)   { A->IA[i]++; P->IA[i]++; }
    for (i=0;i<nnzA;++i) { A->JA[i]++; }
    for (i=0;i<=nc;++i)  { Pt->IA[i]++; }
    for (i=0;i<nnzP;++i) { P->JA[i]++; Pt->JA[i]++; }
    
    // compute P' A P
    dCSRmat PtAP = fasp_blas_dcsr_rap2(Pt->IA,Pt->JA,Pt->val,A->IA,A->JA,A->val,
                                       Pt->IA,Pt->JA,Pt->val,n,nc,&maxrpout,
                                       P->IA,P->JA);
    
    Ac->row = PtAP.row;
    Ac->col = PtAP.col;
    Ac->nnz = PtAP.nnz;
    Ac->IA  = PtAP.IA;
    Ac->JA  = PtAP.JA;
    Ac->val = PtAP.val;
    
    // shift A back from ltz format
    for (i=0;i<=Ac->row;++i) Ac->IA[i]--;
    for (i=0;i< Ac->nnz;++i) Ac->JA[i]--;
    
    for (i=0;i<=n;++i)   A->IA[i]--;
    for (i=0;i<nnzA;++i) A->JA[i]--;

    for (i=0;i<=n;++i)   P->IA[i]--;
    for (i=0;i<nnzP;++i) P->JA[i]--;
    
    for (i=0;i<=nc;++i)  Pt->IA[i]--;
    for (i=0;i<nnzP;++i) Pt->JA[i]--;
    
    return;
}

/**
 * \fn void fasp_blas_dcsr_rap4 (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, INT *icor_ysk)
 *
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param R   pointer to the dCSRmat matrix
 * \param A   pointer to the dCSRmat matrix
 * \param P   pointer to the dCSRmat matrix
 * \param B   pointer to dCSRmat matrix equal to R*A*P
 * \param icor_ysk pointer to the array
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date   08/02/2011
 *
 * \note Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *          Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 */
void fasp_blas_dcsr_rap4 (dCSRmat *R, 
                          dCSRmat *A, 
                          dCSRmat *P, 
                          dCSRmat *B,
                          INT *icor_ysk) 
{
    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP  
	if ( R->row > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif
    
    if (use_openmp) {
        const INT row=R->row, col=P->col;
        INT *ir=R->IA, *ia=A->IA, *ip=P->IA;
        INT *jr=R->JA, *ja=A->JA, *jp=P->JA;
        REAL *rj=R->val, *aj=A->val, *pj=P->val;
        INT istart, iistart;
        INT end_row, end_rowA, end_rowR;
        INT i, j, jj, k, length, myid, mybegin, myend;
        INT jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
        INT *index = NULL;
        INT *iindex = NULL;
        INT *BTindex = NULL;
        REAL *temp = NULL;
        INT FiveMyid, min_A, min_P, A_pos, P_pos, FiveIc;
        INT minus_one_length_A = icor_ysk[5*nthreads];
        INT minus_one_length_P = icor_ysk[5*nthreads+1];
        INT minus_one_length = minus_one_length_A + minus_one_length_P;
    
        INT *iindexs = (INT *)fasp_mem_calloc(minus_one_length+minus_one_length_P, sizeof(INT));
#if CHMEM_MODE
        total_alloc_mem += minus_one_length*sizeof(INT);
#endif
        INT *indexs = iindexs + minus_one_length_P;
        INT *BTindexs = indexs + minus_one_length_A;
        
        INT *iac=(INT*)fasp_mem_calloc(row+1,sizeof(INT));
#if CHMEM_MODE
        total_alloc_mem += (row+1)*sizeof(INT);
#endif
        INT *part_end=(INT*)fasp_mem_calloc(2*nthreads+row,sizeof(INT));
#if CHMEM_MODE
        total_alloc_mem += (2*nthreads+row)*sizeof(INT);
#endif
        INT *iac_temp=part_end + nthreads;
        INT **iindex_array = (INT **)fasp_mem_calloc(nthreads, sizeof(INT *));
        INT **index_array = (INT **)fasp_mem_calloc(nthreads, sizeof(INT *));
    
        fasp_iarray_set(minus_one_length, iindexs, -2);
#ifdef _OPENMP 
#pragma omp parallel for private(myid, FiveMyid, mybegin, myend, min_A, min_P, index, iindex, A_pos, P_pos, ic, FiveIc, jj_counter, jj_row_begining, end_rowR, jj1, i1, end_rowA, jj2, i2, end_row, jj3, i3)
#endif
        for (myid = 0; myid < nthreads; myid ++) {
            FiveMyid = myid * 5;
            mybegin = icor_ysk[FiveMyid];
            if (myid == nthreads-1) {
                myend = row;
            }
            else {
                myend = icor_ysk[FiveMyid+5];
            }
            min_A = icor_ysk[FiveMyid+2];
            min_P = icor_ysk[FiveMyid+4];
            A_pos = 0;
            P_pos = 0;
            for (ic = myid-1; ic >= 0; ic --) {
                FiveIc = ic * 5;
                A_pos += icor_ysk[FiveIc+1];
                P_pos += icor_ysk[FiveIc+3];
            }
            iindex_array[myid] = iindex= iindexs + P_pos - min_P;
            index_array[myid] = index = indexs + A_pos - min_A;
            jj_counter = 0;
            for (ic = mybegin; ic < myend; ic ++) {
                iindex[ic] = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter ++;
                end_rowR = ir[ic+1];
                for (jj1 = ir[ic]; jj1 < end_rowR; jj1 ++) {
                    i1 = jr[jj1];
                    end_rowA = ia[i1+1];
                    for (jj2 = ia[i1]; jj2 < end_rowA; jj2 ++) {
                        i2 = ja[jj2];
                        if (index[i2] != ic) {
                            index[i2] = ic;
                            end_row = ip[i2+1];
                            for (jj3 = ip[i2]; jj3 < end_row; jj3 ++) {
                                i3 = jp[jj3];
                                if (iindex[i3] < jj_row_begining) {
                                    iindex[i3] = jj_counter;
                                    jj_counter ++;
                                }
                            }
                        }
                    }
                }
                iac_temp[ic+myid] = jj_row_begining;
            }
            iac_temp[myend+myid] = jj_counter;
            part_end[myid] = myend + myid + 1;
        }
        fasp_iarray_cp(part_end[0], iac_temp, iac);
        jj_counter = part_end[0];
        INT Ctemp = 0;
        for (i1 = 1; i1 < nthreads; i1 ++) {
            Ctemp += iac_temp[part_end[i1-1]-1];
            for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++) {
                iac[jj_counter] = iac_temp[jj1] + Ctemp;
                jj_counter ++;
            }
        }
        INT *jac=(INT*)fasp_mem_calloc(iac[row],sizeof(INT));
#if CHMEM_MODE
        total_alloc_mem += iac[row]*sizeof(INT);
#endif
        fasp_iarray_set(minus_one_length, iindexs, -2);
#ifdef _OPENMP 
#pragma omp parallel for private(myid, index, iindex, FiveMyid, mybegin, myend, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, iistart, end_row)
#endif
        for (myid = 0; myid < nthreads; myid ++) {
            iindex = iindex_array[myid];
            index = index_array[myid];
            FiveMyid = myid * 5;
            mybegin = icor_ysk[FiveMyid];
            if (myid == nthreads-1) {
                myend = row;
            }
            else {
                myend = icor_ysk[FiveMyid+5];
            }
            for (i = mybegin; i < myend; ++ i) {
                istart = -1;
                length = 0;
                i1 = i+1;
                // go across the rows in R
                end_rowR = ir[i1];
                for (jj = ir[i]; jj < end_rowR; ++ jj) {
                    j = jr[N2C(jj)];
                    // for each column in A
                    end_rowA = ia[j+1];
                    for (k = ia[j]; k < end_rowA; ++ k) {
                        if (index[N2C(ja[N2C(k)])] == -2) {
                            index[N2C(ja[N2C(k)])] = istart;
                            istart = ja[N2C(k)];
                            ++ length;
                        }
                    }
                }
                // book-keeping [reseting length and setting iistart]
                //count = length;
                iistart = -1;
                //length = 0;
                // use each column that would have resulted from R*A
                //for (j = 0; j < count; ++ j) {
                for (j = 0; j < length; ++ j) {
                    jj = istart;
                    istart = index[istart];
                    index[N2C(jj)] = -2;
                    // go across the row of P
                    end_row = ip[jj+1];
                    for (k = ip[jj]; k < end_row; ++ k) {
                        // pull out the appropriate columns of P
                        if (iindex[N2C(jp[N2C(k)])] == -2) {
                            iindex[N2C(jp[N2C(k)])] = iistart;
                            iistart = jp[N2C(k)];
                            //++length;
                        }
                    } // end for k
                } // end for j
                // put the correct columns of p into the column list of the products
                end_row = iac[i1];
                for (j = iac[i]; j < end_row; ++ j) {
                    // put the value in B->JA
                    jac[N2C(j)] = iistart;
                    // set istart to the next value
                    iistart = iindex[N2C(iistart)];
                    // set the iindex spot to 0
                    iindex[N2C(jac[j])] = -2;
                } // end j
            }
        }
        // Third loop: compute entries of R*A*P
        REAL *acj=(REAL*)fasp_mem_calloc(iac[row],sizeof(REAL));
#if CHMEM_MODE
        total_alloc_mem += iac[row]*sizeof(REAL);
#endif
        REAL *temps=(REAL*)fasp_mem_calloc(minus_one_length_A,sizeof(REAL));
#if CHMEM_MODE
        total_alloc_mem += minus_one_length_A*sizeof(REAL);
#endif

#ifdef _OPENMP 
#pragma omp parallel for private(myid, index, FiveMyid, mybegin, myend, min_A, min_P, A_pos, P_pos, ic, FiveIc, BTindex, temp, i, i1, end_row, j, istart, length, end_rowR, jj, end_rowA, k)
#endif
        for (myid = 0; myid < nthreads; myid ++) {
            index = index_array[myid];
            FiveMyid = myid * 5;
            mybegin = icor_ysk[FiveMyid];
            if (myid == nthreads-1) {
                myend = row;
            }
            else {
                myend = icor_ysk[FiveMyid+5];
            }
            min_A = icor_ysk[FiveMyid+2];
            min_P = icor_ysk[FiveMyid+4];
            A_pos = 0;
            P_pos = 0;
            for (ic = myid-1; ic >= 0; ic --) {
                FiveIc = ic * 5;
                A_pos += icor_ysk[FiveIc+1];
                P_pos += icor_ysk[FiveIc+3];
            }
            BTindex = BTindexs + P_pos - min_P;
            temp = temps + A_pos - min_A;
            for (i = mybegin; i < myend; ++ i) {
                i1 = i+1;
                // each col of B
                end_row = iac[i1];
                for (j = iac[i]; j < end_row; ++ j) {
                    BTindex[N2C(jac[N2C(j)])] = j;
                }
                // reset istart and length at the begining of each loop
                istart = -1;
                length = 0;
                // go across the rows in R
                end_rowR = ir[i1];
                for (jj = ir[i]; jj < end_rowR; ++ jj) {
                    j = jr[N2C(jj)];
                    // for each column in A
                    end_rowA = ia[j+1];
                    for (k = ia[j]; k < end_rowA; ++ k) {
                        if (index[N2C(ja[N2C(k)])] == -2) {
                            index[N2C(ja[N2C(k)])] = istart;
                            istart = ja[N2C(k)];
                            ++ length;
                        }
                        temp[N2C(ja[N2C(k)])] += rj[N2C(jj)]*aj[N2C(k)];
                    }
                }
                // book-keeping [reseting length and setting iistart]
                // use each column that would have resulted from R*A
                for (j = 0; j < length; ++ j) {
                    jj = N2C(istart);
                    istart = index[N2C(istart)];
                    index[N2C(jj)] = -2;
                    // go across the row of P
                    end_row = ip[jj+1];
                    for (k = ip[jj]; k < end_row; ++ k) {
                        // pull out the appropriate columns of P
                        acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
                    }
                    temp[jj]=0.0;
                }
            }
        }
        // setup coarse matrix B
        B->row = row;
        B->col = col;
        B->IA = iac;
        B->JA = jac;
        B->val = acj;
        B->nnz = B->IA[B->row] - B->IA[0];
        fasp_mem_free(temps);
        fasp_mem_free(iindexs);
        fasp_mem_free(part_end);
        fasp_mem_free(iindex_array);
        fasp_mem_free(index_array);
    }
    else {
        fasp_blas_dcsr_rap (R, A, P, B);
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
