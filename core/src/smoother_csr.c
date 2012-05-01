/*! \file smoother_csr.c
 *  \brief Smoothers for sparse matrix in CSR format
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dcsr_jacobi (dvector *u, const INT i_1, const INT i_n, const INT s, 
 *                                     dCSRmat *A, dvector *b, INT L)
 *
 * \brief Jacobi method as a smoother
 *
 * \param u    Initial guess (in) and the new approximation after smoothing
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to stiffness matrix
 * \param b    Pointer to right hand side
 * \param L    Number of iterations
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/26/2009
 */
void fasp_smoother_dcsr_jacobi (dvector *u, 
                                const INT i_1, 
                                const INT i_n, 
                                const INT s, 
                                dCSRmat *A, 
                                dvector *b, 
                                INT L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT i,j,k,begin_row,end_row;
    
    REAL *t = (REAL *)fasp_mem_calloc(N,sizeof(REAL));    
#if CHMEM_MODE
    total_alloc_mem += (N)*sizeof(REAL);
#endif    
    
    REAL *d = (REAL *)fasp_mem_calloc(N,sizeof(REAL));
#if CHMEM_MODE
    total_alloc_mem += (N)*sizeof(REAL);
#endif    
    
    while (L--) {
    
    if (s>0) {
    for (i=i_1;i<=i_n;i+=s) {
    t[i]=bval[i];
    begin_row=ia[i],end_row=ia[i+1];
    for (k=begin_row;k<end_row;++k) {
    j=ja[k];
    if (i!=j) t[i]-=aj[k]*uval[j];
    else d[i]=aj[k];
    }    
    }
    
    for (i=i_1;i<=i_n;i+=s) {    
    if (ABS(d[i])>SMALLREAL) uval[i]=t[i]/d[i];
    }    
    } 
    
    else {
    
    for (i=i_1;i>=i_n;i+=s) {
    t[i]=bval[i];
    begin_row=ia[i],end_row=ia[i+1];
    for (k=begin_row;k<end_row;++k) {
    j=ja[k];
    if (i!=j) t[i]-=aj[k]*uval[j];
    else d[i]=aj[k];
    }
    
    }
    
    for (i=i_1;i>=i_n;i+=s) {
    if (ABS(d[i])>SMALLREAL) uval[i]=t[i]/d[i];
    }
    } 
    
    } // end while
    
    fasp_mem_free(t);
    fasp_mem_free(d);
    
    return;    
}

/**
 * \fn void fasp_smoother_dcsr_gs (dvector *u, const INT i_1, const INT i_n, const INT s, 
 *                                 dCSRmat *A, dvector *b, INT L)
 *
 * \brief Gauss-Seidel method as a smoother
 *
 * \param u    Initial guess (in) and the new approximation after smoothing
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to stiffness matrix
 * \param b    Pointer to right hand side
 * \param L    Number of iterations
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/26/2009
 */
void fasp_smoother_dcsr_gs (dvector *u, 
                            const INT i_1, 
                            const INT i_n, 
                            const INT s, 
                            dCSRmat *A, 
                            dvector *b, 
                            INT L)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0.0;
    
    if (s > 0) {
    
    while (L--) {
            
    for (i=i_1;i<=i_n;i+=s) {
    
                t = bval[i];
    begin_row=ia[i],end_row=ia[i+1];
                
#if DIAGONAL_PREF // diagonal first
                d=aj[begin_row];
                if (ABS(d)>SMALLREAL) {
                    for (k=begin_row+1;k<end_row;++k) {
                        j=ja[k];
                        t-=aj[k]*uval[j]; 
                    }
                    uval[i]=t/d;
                }                    
#else // general order
                for (k=begin_row;k<end_row;++k) {
    j=ja[k];
    if (i!=j) 
    t-=aj[k]*uval[j]; 
    else if (ABS(aj[k])>SMALLREAL) d=1.e+0/aj[k];    
    }
                uval[i]=t*d;
#endif
                
            } // end for i
        } // end while    
        
    } // if s
    else {    
        
        while (L--) {

            for (i=i_1;i>=i_n;i+=s) {
            
                t=bval[i];
                begin_row=ia[i],end_row=ia[i+1];

#if DIAGONAL_PREF // diagonal first
                d=aj[begin_row];
                if (ABS(d)>SMALLREAL) {
                    for (k=begin_row+1;k<end_row;++k) {
                        j=ja[k];
                        t-=aj[k]*uval[j]; 
                    }
                    uval[i]=t/d;
                }                    
#else // general order
                for (k=begin_row;k<end_row;++k) {
    j=ja[k];
    if (i!=j) 
    t-=aj[k]*uval[j]; 
    else if (ABS(aj[k])>SMALLREAL) d=1.e+0/aj[k];    
    }
                uval[i]=t*d;
#endif

            } // end for i
        } // end while    
        
    } // end if
    return;
}

/**
 * \fn void fasp_smoother_dcsr_gs_cf (dvector *u, dCSRmat *A, dvector *b, INT L, INT *mark, const INT order)
 *
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \param u      Initial guess (in) and the new approximation after smoothing
 * \param A      Pointer to stiffness matrix
 * \param b      Pointer to right hand side
 * \param L      Number of iterations
 * \param mark   C/F marker array
 * \param order  C/F ordering: -1: F-first; 1: C-first
 *
 * \author Zhiyang Zhou
 * \date   2010/11/12 
 */
void fasp_smoother_dcsr_gs_cf (dvector *u, 
                               dCSRmat *A, 
                               dvector *b, 
                               INT L, 
                               INT *mark, 
                               const INT order )
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT    i,j,k,begin_row,end_row;
    INT    size = b->row;
    REAL   t,d=0.0;
    
    if (order == -1) // F-point first
    {    
        while (L--) 
        {
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] != 1)
                {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
            
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] == 1)
                {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
        } // end while    
    }
    else
    {
        while (L--) 
        {
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] == 1)
                {
                    t = bval[i];
                    begin_row = ia[i],end_row = ia[i+1]-1;
                    for (k = begin_row; k <= end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
            
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] != 1)
                {
                    t = bval[i];
                    begin_row = ia[i],end_row = ia[i+1]-1;
                    for (k = begin_row; k <= end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d) > SMALLREAL) uval[i] = t/d;
                }
            } // end for i
        } // end while    
    }    
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_sgs (dvector *u, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Symmetric Gauss-Seidel method as a smoother
 *
 * \param u    Initial guess (in) and the new approximation after smoothing
 * \param A    Pointer to stiffness matrix
 * \param b    Pointer to right hand side
 * \param L    Number of iterations
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
void fasp_smoother_dcsr_sgs (dvector *u, 
                             dCSRmat *A, 
                             dvector *b, 
                             INT L)
{
    const INT    nm1=b->row-1;
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d;
    
    while (L--) {
        
        // forward sweep
        for (i=0;i<=nm1;++i) {
            t=bval[i];
            begin_row=ia[i], end_row=ia[i+1];
            for (k=begin_row;k<end_row;++k) {
                j=ja[k];
                if (i!=j) t-=aj[k]*uval[j];
                else d=aj[k];
            } // end for k
            if (ABS(d)>SMALLREAL) uval[i]=t/d;
        } // end for i
        
        // backward sweep
        for (i=nm1-1;i>=0;--i) {
            t=bval[i];
            begin_row=ia[i], end_row=ia[i+1];
            for (k=begin_row;k<end_row;++k) {
                j=ja[k];
                if (i!=j) t-=aj[k]*uval[j];
                else d=aj[k];
            } // end for k
            if (ABS(d)>SMALLREAL) uval[i]=t/d;
        } // end for i
        
    } // end while    
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_sor (dvector *u, const INT i_1, const INT i_n, const INT s, 
 *                                  dCSRmat *A, dvector *b, INT L, const REAL w)
 *
 * \brief SOR method as a smoother
 *
 * \param u    Initial guess (in) and the new approximation after smoothing
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to stiffness matrix
 * \param b    Pointer to right hand side
 * \param L    Number of iterations
 * \param w    Over-relaxation weight
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 */
void fasp_smoother_dcsr_sor (dvector *u, 
                             const INT i_1, 
                             const INT i_n, 
                             const INT s, 
                             dCSRmat *A, 
                             dvector *b, 
                             INT L, 
                             const REAL w)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;

    // local variables
    INT    i,j,k,begin_row,end_row;
    REAL   t, d;
    
    while (L--) {
        if (s>0) {
            for (i=i_1;i<=i_n;i+=s)
            {
                t=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k)
                {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else
                        d=aj[k];
                }
                if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
            }
        } 
        else {
            for (i=i_1;i>=i_n;i+=s)
            {
                t=bval[i];
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k)
                {
                    j=ja[k];
                    if (i!=j)
                        t-=aj[k]*uval[j];
                    else
                        d=aj[k];
                }
                if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
            }
        } 
    }  // end while
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_sor_cf (dvector *u, dCSRmat *A, dvector *b, INT L, 
 *                                     const REAL w, INT *mark, const INT order)
 *
 * \brief SOR smoother with C/F ordering for Au=b
 *
 * \param u      Initial guess (in) and the new approximation after smoothing
 * \param A      Pointer to stiffness matrix
 * \param b      Pointer to right hand side
 * \param L      Number of iterations
 * \param w      Over-relaxation weight
 * \param mark   C/F marker array
 * \param order  C/F ordering: -1: F-first; 1: C-first
 *
 * \author Zhiyang Zhou
 * \date   2010/11/12 
 */
void fasp_smoother_dcsr_sor_cf (dvector *u, 
                                dCSRmat *A, 
                                dvector *b, 
                                INT L, 
                                const REAL w, 
                                INT *mark, 
                                const INT order )
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT    i,j,k,begin_row,end_row;
    INT    size = b->row;
    REAL   t,d=0.0;
    
    if (order == -1) // F-point first
    {    
        while (L--) 
        {
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] == 0 || mark[i] == 2)
                {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
            
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] == 1)
                {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
        } // end while    
    }
    else
    {
        while (L--) 
        {
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] == 1)
                {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
            
            for (i = 0; i < size; i ++) 
            {
                if (mark[i] != 1)
                {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i+1];
                    for (k = begin_row; k < end_row; k ++) 
                    {
                        j = ja[k];
                        if (i!=j) t -= aj[k]*uval[j]; 
                        else d = aj[k];    
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
            } // end for i
        } // end while    
    }    
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_ilu (dCSRmat *A, dvector *b, dvector *x, void *data)
 *
 * \brief ILU method as a smoother
 *
 * \param A     Pointer to stiffness matrix
 * \param b     Pointer to right hand side
 * \param x     Pointer to current solution
 * \param data  Pointer to user defined data
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   2010/11/12 
 */
void fasp_smoother_dcsr_ilu (dCSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             void *data)
{
    const unsigned INT m=A->row, m2=2*m, memneed=3*m;
    const ILU_data *iludata=(ILU_data *)data;    
    
    REAL *zz = iludata->work; 
    REAL *zr = iludata->work+m;
    REAL *z  = iludata->work+m2;
    
    if (iludata->nwork<memneed) goto MEMERR; 

    {
        INT i, j, jj, begin_row, end_row;
        REAL *lu = iludata->luval;
        INT *ijlu = iludata->ijlu;
        REAL *xval = x->val, *bval = b->val;
        
        /** form residual zr = b - A x */
        fasp_array_cp(m,bval,zr); fasp_blas_dcsr_aAxpy(-1.0,A,xval,zr);
        
        // forward sweep: solve unit lower matrix equation L*zz=zr
        zz[0]=zr[0];
        for (i=1;i<m;++i) {
            begin_row=ijlu[i]; end_row=ijlu[i+1];
            for (j=begin_row;j<end_row;++j) {
                jj=ijlu[j];
                if (jj<i) zr[i]-=lu[j]*zz[jj];
                else break;
            }
            zz[i]=zr[i];
        }
        
        // backward sweep: solve upper matrix equation U*z=zz
        z[m-1]=zz[m-1]*lu[m-1];
        for (i=m-2;i>=0;--i) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=end_row;j>=begin_row;--j) {
                jj=ijlu[j];
                if (jj>i) zz[i]-=lu[j]*z[jj];
                else break;
            } 
            z[i]=zz[i]*lu[i];
        }
        
        fasp_blas_array_axpy(m,1,z,xval);    
    }
    
    return;
    
MEMERR:
    printf("### ERROR: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
    exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_smoother_dcsr_kaczmarz (dvector *u, const INT i_1, const INT i_n, const INT s, 
 *                                       dCSRmat *A, dvector *b, INT L, const REAL w)
 *
 * \brief Kaczmarz method as a smoother
 *
 * \param u    Initial guess (in) and the new approximation after smoothing
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to stiffness matrix
 * \param b    Pointer to right hand side
 * \param L    Number of iterations
 * \param w    Relaxation weight
 *
 * \author Xiaozhe Hu
 * \date   2010/11/12 
 */
void fasp_smoother_dcsr_kaczmarz (dvector *u, 
                                  const INT i_1, 
                                  const INT i_n, 
                                  const INT s, 
                                  dCSRmat *A, 
                                  dvector *b, 
                                  INT L, 
                                  const REAL w)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  temp1,temp2,alpha;
    
    if (s > 0) {
        
        while (L--) {
            for (i=i_1;i<=i_n;i+=s) {
                temp1 = 0; temp2 = 0;
                
                begin_row=ia[i], end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    temp1 += aj[k]*aj[k];
                    temp2 += aj[k]*uval[j];
                } // end for k
                
                alpha = (bval[i] - temp2)/temp1;
                
                for (k=begin_row;k<end_row;++k){
                    j = ja[k];
                    uval[j] += w*alpha*aj[k];
                }// end for k
            } // end for i
        } // end while    
        
    } // if s
    
    else {    
        
        while (L--) {
            for (i=i_1;i>=i_n;i+=s) {
                temp1 = 0; temp2 = 0;
                
                begin_row=ia[i], end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    temp1 += aj[k]*aj[k];
                    temp2 += aj[k]*uval[j];
                } // end for k
                
                alpha = (bval[i] - temp2)/temp1;
                
                for (k=begin_row;k<end_row;++k){
                    j = ja[k];
                    uval[j] += w*alpha*aj[k];
                }// end for k
            } // end for i
        } // end while    
        
    } // end if
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_L1diag (dvector *u, const INT i_1, const INT i_n, const INT s, 
 *                                     dCSRmat *A, dvector *b, INT L)
 *
 * \brief Diagonal scaling (using L1 norm) as a smoother
 *
 * \param u    Initial guess (in) and the new approximation after smoothing
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to stiffness matrix
 * \param b    Pointer to right hand side
 * \param L    Number of iterations
 *
 * \author Xiaozhe Hu, James Brannick
 * \date   01/26/2011
 */
void fasp_smoother_dcsr_L1diag (dvector *u, 
                                const INT i_1, 
                                const INT i_n, 
                                const INT s, 
                                dCSRmat *A, 
                                dvector *b, 
                                INT L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    
    // Checks should be outside of for; t,d can be allocated before calling!!! --Chensong
    REAL *t = (REAL *)fasp_mem_calloc(N,sizeof(REAL));    
#if CHMEM_MODE
    total_alloc_mem += (N)*sizeof(REAL);
#endif    
    REAL *d = (REAL *)fasp_mem_calloc(N,sizeof(REAL));
#if CHMEM_MODE
    total_alloc_mem += (N)*sizeof(REAL);
#endif    
    
    while (L--) {
        
        if (s>0) {
            for (i=i_1;i<=i_n;i+=s) {
                t[i]=bval[i]; d[i]=0.0;
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    t[i]-=aj[k]*uval[j];
                    d[i]+=ABS(aj[k]);
                }    
            }
            
            for (i=i_1;i<=i_n;i+=s) {    
                if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
            }    
        } 
        
        else {
            
            for (i=i_1;i>=i_n;i+=s) {
                t[i]=bval[i];d[i]=0.0;
                begin_row=ia[i],end_row=ia[i+1];
                for (k=begin_row;k<end_row;++k) {
                    j=ja[k];
                    t[i]-=aj[k]*uval[j];
                    d[i]+=ABS(aj[k]);
                }
                
            }
            
            for (i=i_1;i>=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
            }
        } 
        
    } // end while
    
    fasp_mem_free(t);
    fasp_mem_free(d);
    
    return;    
}

#if 0
/**
 * \fn dCSRmat static form_contractor(dCSRmat *A, const INT smoother, const INT steps, 
 *                                    const INT ndeg, const REAL relax, const REAL dtol)
 *
 * \brief Form contractor I-BA
 *
 * \param A          Pointer to the dCSRmat
 * \param smoother   Smoother type
 * \param steps      Smoothing steps
 * \param ndeg       Degree of the polynomial smoother
 * \param relax      Relaxation parameter for SOR smoother
 * \param dtol       Drop tplerance for droping small entries in matrix
 *
 * \return The contractor in dCSRmat format
 *
 * \author Xiaozhe Hu, James Brannick
 * \date   01/26/2011
 *
 * \note Not an O(N) algorithm, need to be modified!!!!
 */
static dCSRmat form_contractor (dCSRmat *A, 
                                const INT smoother, 
                                const INT steps, 
                                const INT ndeg, 
                                const REAL relax, 
                                const REAL dtol)
{
    const INT    n=A->row;
    unsigned INT i;
    
    REAL *work = (REAL *)fasp_mem_calloc(2*n,sizeof(REAL));
    
#if CHMEM_MODE
    total_alloc_mem += (2*n)*sizeof(REAL);
#endif    
    
    dvector b, x;
    b.row=x.row=n;
    b.val=work; x.val=work+n;
    
    INT *index = (INT *)fasp_mem_calloc(n,sizeof(INT));    
    
#if CHMEM_MODE
    total_alloc_mem += (n)*sizeof(INT);
#endif
    
    for (i=0; i<n; ++i) index[i]=i;
    
    dCSRmat B = fasp_dcsr_create(n, n, n*n); // too much memory required, need to change!!
    dCSRmat C, D;
    
    for (i=0; i<n; ++i){
        
        // get i-th column
        fasp_dcsr_getcol(i, A, b.val);
        
        // set x =0.0 
        fasp_dvec_set(n, &x, 0.0);
        
        // smooth
        switch (smoother) {
            case GS:
                fasp_smoother_dcsr_gs(&x, 0, n-1, 1, A, &b, steps);
                break;
            case POLY:
                fasp_smoother_dcsr_poly(A, &b, &x, n, ndeg, steps); 
                break;
            case JACOBI:
                fasp_smoother_dcsr_jacobi(&x, 0, n-1, 1, A, &b, steps);
                break;
            case SGS:
                fasp_smoother_dcsr_sgs(&x, A, &b, steps);
                break;
            case SOR:
                fasp_smoother_dcsr_sor(&x, 0, n-1, 1, A, &b, steps, relax);
                break;
            case SSOR:
                fasp_smoother_dcsr_sor(&x, 0, n-1, 1, A, &b, steps, relax);
                fasp_smoother_dcsr_sor(&x, n-1, 0,-1, A, &b, steps, relax);
                break;
            case GSOR:
                fasp_smoother_dcsr_gs(&x, 0, n-1, 1, A, &b, steps);
                fasp_smoother_dcsr_sor(&x, n-1, 0, -1, A, &b, steps, relax);
                break;
            case SGSOR:
                fasp_smoother_dcsr_gs(&x, 0, n-1, 1, A, &b, steps);
                fasp_smoother_dcsr_gs(&x, n-1, 0,-1, A, &b, steps);
                fasp_smoother_dcsr_sor(&x, 0, n-1, 1, A, &b, steps, relax);
                fasp_smoother_dcsr_sor(&x, n-1, 0,-1, A, &b, steps, relax);
                break;
            default:
                printf("### ERROR: Wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
        } 
        
        // store to B
        B.IA[i] = i*n;
        memcpy(&(B.JA[i*n]), index, n*sizeof(INT));
        memcpy(&(B.val[i*n]), x.val, x.row*sizeof(REAL));
        
    }
    
    B.IA[n] = n*n;
    
    // drop small entries
    compress_dCSRmat(&B, &D, dtol);
    
    // get contractor
    fasp_dcsr_trans(&D, &C);
    
    // clean up
    fasp_mem_free(work);
    fasp_dcsr_free(&B);
    fasp_dcsr_free(&D);
    
    return C;
}
#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
