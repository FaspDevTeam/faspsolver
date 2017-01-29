/*! \file ItrSmootherCSR.c
 *
 *  \brief Smoothers for dCSRmat matrices
 *
 *  \note This file contains Level-2 (Itr) functions. It requires
 *        AuxArray.c, AuxMemory.c, AuxMessage.c, BlaArray.c, and BlaSpmvCSR.c
 */

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/

static void swep2db (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
                     INT nbegy, INT *mark, INT nx, INT ny);
static void swep3db (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
                     INT nbegy, INT nbegz, INT *mark, INT nx, INT ny, INT nz);
static void rb0b2d (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f,
                    INT *mark, INT nx, INT ny, INT nsweeps);
static void rb0b3d (INT *ia, INT *ja, REAL *aa,REAL *u, REAL *f,
                    INT *mark, INT nx, INT ny, INT nz, INT nsweeps);
static void swep2df (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
                     INT nbegy, INT *mark, INT nx, INT ny);
static void swep3df (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
                     INT nbegy, INT nbegz, INT *mark, INT nx, INT ny, INT nz);
static void rb0f2d (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f,
                    INT *mark, INT nx, INT ny, INT nsweeps);
static void rb0f3d (INT *ia, INT *ja, REAL *aa,REAL *u, REAL *f, INT *mark,
                    INT nx, INT ny, INT nz, INT nsweeps);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dcsr_jacobi (dvector *u, const INT i_1, const INT i_n,
 *                                     const INT s, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Jacobi method as a smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/26/2009
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/29/2012
 */
void fasp_smoother_dcsr_jacobi (dvector    *u,
                                const INT   i_1,
                                const INT   i_n,
                                const INT   s,
                                dCSRmat    *A,
                                dvector    *b,
                                INT         L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    REAL         w = 0.8;
    
    // local variables
    INT i,j,k,begin_row,end_row;
    
    // OpenMP variables
#ifdef _OPENMP
    INT myid, mybegin, myend;
    INT nthreads = fasp_get_num_threads();
#endif
    
    REAL *t = (REAL *)fasp_mem_calloc(N,sizeof(REAL));
    REAL *d = (REAL *)fasp_mem_calloc(N,sizeof(REAL));
    
    while (L--) {
        if (s>0) {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, begin_row, end_row, i, k, j)
                for (myid=0; myid<nthreads; ++myid) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin += i_1; myend += i_1;
                    for (i=mybegin; i<myend; i+=s) {
                        t[i]=bval[i];
                        begin_row=ia[i],end_row=ia[i+1];
                        for (k=begin_row; k<end_row; ++k) {
                            j=ja[k];
                            if (i!=j) t[i]-=aj[k]*uval[j];
                            else d[i]=aj[k];
                        }
                    }
                }
            }
            else {
#endif
                for (i=i_1;i<=i_n;i+=s) {
                    t[i]=bval[i];
                    begin_row=ia[i],end_row=ia[i+1];
                    for (k=begin_row;k<end_row;++k) {
                        j=ja[k];
                        if (i!=j) t[i]-=aj[k]*uval[j];
                        else d[i]=aj[k];
                    }
                }
#ifdef _OPENMP
            }
#endif
            
#ifdef _OPENMP
#pragma omp parallel for private (i)
#endif
            for (i=i_1;i<=i_n;i+=s) {
                //if (ABS(d[i])>SMALLREAL) uval[i]= w*t[i]/d[i];
                if (ABS(d[i])>SMALLREAL) uval[i]=(1-w)*uval[i]+ w*t[i]/d[i];
            }
        }
        else {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, begin_row, end_row, k, j)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin = i_1-mybegin; myend = i_1-myend;
                    for (i=mybegin; i>myend; i+=s) {
                        t[i]=bval[i];
                        begin_row=ia[i],end_row=ia[i+1];
                        for (k=begin_row; k<end_row; ++k) {
                            j=ja[k];
                            if (i!=j) t[i]-=aj[k]*uval[j];
                            else d[i]=aj[k];
                        }
                    }
                }
            }
            else {
#endif
                for (i=i_1;i>=i_n;i+=s) {
                    t[i]=bval[i];
                    begin_row=ia[i],end_row=ia[i+1];
                    for (k=begin_row;k<end_row;++k) {
                        j=ja[k];
                        if (i!=j) t[i]-=aj[k]*uval[j];
                        else d[i]=aj[k];
                    }
                }
#ifdef _OPENMP
            }
#endif
            
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
            for (i=i_1;i>=i_n;i+=s) {
                if (ABS(d[i])>SMALLREAL) uval[i]=(1-w)*uval[i]+ w*t[i]/d[i];
            }
        }
        
    } // end while
    
    fasp_mem_free(t);
    fasp_mem_free(d);
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_gs (dvector *u, const INT i_1, const INT i_n,
 *                                 const INT s, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Gauss-Seidel method as a smoother
 *
 * \param u    Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1  Starting index
 * \param i_n  Ending index
 * \param s    Increasing step
 * \param A    Pointer to dBSRmat: the coefficient matrix
 * \param b    Pointer to dvector: the right hand side
 * \param L    Number of iterations
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/26/2009
 *
 * Modified by Chunsheng Feng, Zheng Li on 09/01/2012
 */
void fasp_smoother_dcsr_gs (dvector    *u,
                            const INT   i_1,
                            const INT   i_n,
                            const INT   s,
                            dCSRmat    *A,
                            dvector    *b,
                            INT         L)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0.0;
    
#ifdef _OPENMP
    // variables for OpenMP
    const INT    N = ABS(i_n - i_1)+1;
    INT   myid, mybegin, myend;
    INT   nthreads = fasp_get_num_threads();
#endif
    
    if (s > 0) {
        while (L--) {
#ifdef _OPENMP
            if (N >OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, t, begin_row, end_row, d, k, j)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin += i_1, myend += i_1;
                    for (i=mybegin; i<myend; i+=s) {
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
                }
            }
            else {
#endif
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
#ifdef _OPENMP
            }
#endif
        } // end while
        
    } // if s
    else {
        
        while (L--) {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, begin_row, end_row, d, k, j, t)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin = i_1 - mybegin; myend = i_1 - myend;
                    for (i=mybegin; i>myend; i+=s) {
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
                            else if (ABS(aj[k])>SMALLREAL) d=1.0/aj[k];
                        }
                        uval[i]=t*d;
#endif
                    } // end for i
                }
            }
            else {
#endif
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
                        else if (ABS(aj[k])>SMALLREAL) d=1.0/aj[k];
                    }
                    uval[i]=t*d;
#endif
                } // end for i
#ifdef _OPENMP
            }
#endif
        } // end while
    } // end if
    return;
}

/**
 * \fn void fasp_smoother_dcsr_gs_cf (dvector *u, dCSRmat *A, dvector *b, INT L,
 *                                    INT *mark, const INT order)
 *
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 * \param mark   C/F marker array
 * \param order  C/F ordering: -1: F-first; 1: C-first
 *
 * \author Zhiyang Zhou
 * \date   11/12/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/24/2012
 */
void fasp_smoother_dcsr_gs_cf (dvector   *u,
                               dCSRmat   *A,
                               dvector   *b,
                               INT        L,
                               INT       *mark,
                               const INT  order)
{
    const INT    nrow = b->row; // number of rows
    const INT   *ia = A->IA, *ja = A->JA;
    const REAL  *aj = A->val, *bval = b->val;
    REAL        *uval = u->val;
    
    INT i,j,k,begin_row,end_row;
    REAL t,d=0.0;
    
#ifdef _OPENMP
    INT myid,mybegin,myend;
    INT nthreads = fasp_get_num_threads();
#endif
    
    // F-point first, C-point second
    if (order == FPFIRST) {
        
        while (L--) {
            
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i=mybegin; i<myend; i++) {
                        if (mark[i] != 1) {
                            t = bval[i];
                            begin_row = ia[i], end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                            d = aj[begin_row];
                            for (k = begin_row+1; k < end_row; k ++) {
                                j = ja[k];
                                t -= aj[k]*uval[j];
                            } // end for k
#else
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
#endif // end if DIAG_PREF
                            if (ABS(d) > SMALLREAL) uval[i] = t/d;
                        }
                    } // end for i
                }
            }
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] != 1) {
                        t = bval[i];
                        begin_row = ia[i], end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                        d = aj[begin_row];
                        for (k = begin_row+1; k < end_row; k ++) {
                            j = ja[k];
                            t -= aj[k]*uval[j];
                        } // end for k
#else
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
#endif // end if DIAG_PREF
                        if (ABS(d) > SMALLREAL) uval[i] = t/d;
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
            
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private(myid,mybegin,myend,i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i=mybegin; i<myend; i++) {
                        if (mark[i] == 1) {
                            t = bval[i];
                            begin_row = ia[i], end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                            d = aj[begin_row];
                            for (k = begin_row+1; k < end_row; k ++) {
                                j = ja[k];
                                t -= aj[k]*uval[j];
                            } // end for k
#else
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
#endif // end if DIAG_PREF
                            if (ABS(d) > SMALLREAL) uval[i] = t/d;
                        }
                    } // end for i
                }
            }
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] == 1) {
                        t = bval[i];
                        begin_row = ia[i], end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                        d = aj[begin_row];
                        for (k = begin_row+1; k < end_row; k ++) {
                            j = ja[k];
                            t -= aj[k]*uval[j];
                        } // end for k
#else
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
#endif // end if DIAG_PREF
                        if (ABS(d) > SMALLREAL) uval[i] = t/d;
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
        } // end while
        
    }
    
    // C-point first, F-point second
    else {
        
        while (L--) {
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private(myid,mybegin,myend,t,i,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i=mybegin; i<myend; i++) {
                        if (mark[i] == 1) {
                            t = bval[i];
                            begin_row = ia[i],end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                            d = aj[begin_row];
                            for (k = begin_row+1; k < end_row; k ++) {
                                j = ja[k];
                                t -= aj[k]*uval[j];
                            } // end for k
#else
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
#endif // end if DIAG_PREF
                            if (ABS(d) > SMALLREAL) uval[i] = t/d;
                        }
                    } // end for i
                }
            }
            else {
#endif
                for (i = 0; i < nrow; i ++)  {
                    if (mark[i] == 1) {
                        t = bval[i];
                        begin_row = ia[i],end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 09/22/2012
                        d = aj[begin_row];
                        for (k = begin_row+1; k < end_row; k ++) {
                            j = ja[k];
                            t -= aj[k]*uval[j];
                        } // end for k
#else
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
#endif // end if DIAG_PREF
                        if (ABS(d) > SMALLREAL) uval[i] = t/d;
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
            
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i=mybegin; i<myend; i++) {
                        if (mark[i] != 1) {
                            t = bval[i];
                            begin_row = ia[i],end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 01/17/2013
                            d = aj[begin_row];
                            for (k = begin_row+1; k < end_row; k ++) {
                                j = ja[k];
                                t -= aj[k]*uval[j];
                            } // end for k
#else
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
#endif // end if DIAG_PREF
                            if (ABS(d) > SMALLREAL) uval[i] = t/d;
                        }
                    } // end for i
                }
            }
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] != 1) {
                        t = bval[i];
                        begin_row = ia[i],end_row = ia[i+1];
#if DIAGONAL_PREF // Added by Chensong on 09/22/2012
                        d = aj[begin_row];
                        for (k = begin_row+1; k < end_row; k ++) {
                            j = ja[k];
                            t -= aj[k]*uval[j];
                        } // end for k
#else
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
#endif // end if DIAG_PREF
                        if (ABS(d) > SMALLREAL) uval[i] = t/d;
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
        } // end while
        
    } // end if order
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_sgs (dvector *u, dCSRmat *A, dvector *b, INT L)
 *
 * \brief Symmetric Gauss-Seidel method as a smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 *
 * Modified by Chunsheng Feng, Zheng Li on 09/01/2012
 */
void fasp_smoother_dcsr_sgs (dvector *u,
                             dCSRmat *A,
                             dvector *b,
                             INT      L)
{
    const INT    nm1=b->row-1;
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  t,d=0;
    
#ifdef _OPENMP
    // variables for OpenMP
    INT  myid, mybegin, myend, up;
    INT  nthreads = fasp_get_num_threads();
#endif
    
    while (L--) {
        // forward sweep
#ifdef _OPENMP
        up = nm1 + 1;
        if (up > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, t, begin_row, end_row, j, k, d)
            for (myid=0; myid<nthreads; myid++) {
                fasp_get_start_end(myid, nthreads, up, &mybegin, &myend);
                for (i=mybegin; i<myend; i++) {
                    t=bval[i];
                    begin_row=ia[i], end_row=ia[i+1];
                    for (k=begin_row;k<end_row;++k) {
                        j=ja[k];
                        if (i!=j) t-=aj[k]*uval[j];
                        else d=aj[k];
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=t/d;
                } // end for i
            }
        }
        else {
#endif
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
#ifdef _OPENMP
        }
#endif
        
        // backward sweep
#ifdef _OPENMP
        up = nm1;
        if (up > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, t, begin_row, end_row, k, j, d)
            for (myid=0; myid<nthreads; myid++) {
                fasp_get_start_end(myid, nthreads, up, &mybegin, &myend);
                mybegin = nm1-1-mybegin; myend = nm1-1-myend;
                for (i=mybegin; i>myend; i--) {
                    t=bval[i];
                    begin_row=ia[i], end_row=ia[i+1];
                    for (k=begin_row; k<end_row; k++) {
                        j=ja[k];
                        if (i!=j) t-=aj[k]*uval[j];
                        else d=aj[k];
                    } // end for k
                    if (ABS(d)>SMALLREAL) uval[i]=t/d;
                } // end for i
            }
        }
        else {
#endif
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
#ifdef _OPENMP
        }
#endif
    } // end while
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_sor (dvector *u, const INT i_1, const INT i_n, const INT s,
 *                                  dCSRmat *A, dvector *b, INT L, const REAL w)
 *
 * \brief SOR method as a smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 * \param w      Over-relaxation weight
 *
 * \author Xiaozhe Hu
 * \date   10/26/2010
 *
 * Modified by Chunsheng Feng, Zheng Li on 09/01/2012
 */
void fasp_smoother_dcsr_sor (dvector    *u,
                             const INT   i_1,
                             const INT   i_n,
                             const INT   s,
                             dCSRmat    *A,
                             dvector    *b,
                             INT         L,
                             const REAL  w)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT    i,j,k,begin_row,end_row;
    REAL   t, d=0;
    
#ifdef _OPENMP
    // variables for OpenMP
    const INT    N = ABS(i_n - i_1)+1;
    INT myid, mybegin, myend;
    INT nthreads = fasp_get_num_threads();
#endif
    
    while (L--) {
        if (s>0) {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, t, begin_row, end_row, k, j, d)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin += i_1, myend += i_1;
                    for (i=mybegin; i<myend; i+=s) {
                        t=bval[i];
                        begin_row=ia[i], end_row=ia[i+1];
                        for (k=begin_row; k<end_row; k++) {
                            j=ja[k];
                            if (i!=j)
                                t-=aj[k]*uval[j];
                            else
                                d=aj[k];
                        }
                        if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                    }
                }
                
            }
            else {
#endif
                for (i=i_1; i<=i_n; i+=s) {
                    t=bval[i];
                    begin_row=ia[i], end_row=ia[i+1];
                    for (k=begin_row; k<end_row; ++k) {
                        j=ja[k];
                        if (i!=j)
                            t-=aj[k]*uval[j];
                        else
                            d=aj[k];
                    }
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
#ifdef _OPENMP
            }
#endif
        }
        else {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, t, begin_row, end_row, k, j, d)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin = i_1 - mybegin, myend = i_1 - myend;
                    for (i=mybegin; i>myend; i+=s) {
                        t=bval[i];
                        begin_row=ia[i],end_row=ia[i+1];
                        for (k=begin_row;k<end_row;++k) {
                            j=ja[k];
                            if (i!=j)
                                t-=aj[k]*uval[j];
                            else
                                d=aj[k];
                        }
                        if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                    }
                }
            }
            else {
#endif
                for (i=i_1;i>=i_n;i+=s) {
                    t=bval[i];
                    begin_row=ia[i],end_row=ia[i+1];
                    for (k=begin_row;k<end_row;++k) {
                        j=ja[k];
                        if (i!=j)
                            t-=aj[k]*uval[j];
                        else
                            d=aj[k];
                    }
                    if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                }
#ifdef _OPENMP
            }
#endif
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
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 * \param w      Over-relaxation weight
 * \param mark   C/F marker array
 * \param order  C/F ordering: -1: F-first; 1: C-first
 *
 * \author Zhiyang Zhou
 * \date   2010/11/12
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/29/2012
 */
void fasp_smoother_dcsr_sor_cf (dvector    *u,
                                dCSRmat    *A,
                                dvector    *b,
                                INT         L,
                                const REAL  w,
                                INT        *mark,
                                const INT   order )
{
    const INT    nrow = b->row; // number of rows
    const INT   *ia = A->IA, *ja=A->JA;
    const REAL  *aj = A->val,*bval=b->val;
    REAL        *uval = u->val;
    
    // local variables
    INT    i,j,k,begin_row,end_row;
    REAL   t,d=0.0;
    
#ifdef _OPENMP
    INT    myid, mybegin, myend;
    INT    nthreads = fasp_get_num_threads();
#endif
    
    // F-point first
    if (order == -1) {
        while (L--) {
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private (myid, mybegin, myend, i, t, begin_row, end_row, k, j, d)
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i = mybegin; i < myend; i ++) {
                        if (mark[i] == 0 || mark[i] == 2) {
                            t = bval[i];
                            begin_row = ia[i], end_row = ia[i+1];
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
                            if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                        }
                    }
                }
            } // end for i
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] == 0 || mark[i] == 2) {
                        t = bval[i];
                        begin_row = ia[i], end_row = ia[i+1];
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
                        if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
            
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, i, mybegin, myend, t, begin_row, end_row, k, j, d)
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i = mybegin; i < myend; i++) {
                        if (mark[i] == 1) {
                            t = bval[i];
                            begin_row = ia[i], end_row = ia[i+1];
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
                            if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                        }
                    } // end for i
                }
            }
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] == 1) {
                        t = bval[i];
                        begin_row = ia[i], end_row = ia[i+1];
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
                        if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
        } // end while
    }
    else {
        while (L--) {
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, t, k, j, d, begin_row, end_row)
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i = mybegin; i < myend; i++) {
                        if (mark[i] == 1) {
                            t = bval[i];
                            begin_row = ia[i], end_row = ia[i+1];
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
                            if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                        }
                    } // end for i
                }
            }
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] == 1) {
                        t = bval[i];
                        begin_row = ia[i], end_row = ia[i+1];
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
                        if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
            
#ifdef _OPENMP
            if (nrow > OPENMP_HOLDS) {
#pragma omp parallel for private (myid, mybegin, myend, i, t, begin_row, end_row, k, j, d)
                for (myid = 0; myid < nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, nrow, &mybegin, &myend);
                    for (i = mybegin; i < myend; i++) {
                        if (mark[i] != 1) {
                            t = bval[i];
                            begin_row = ia[i], end_row = ia[i+1];
                            for (k = begin_row; k < end_row; k ++) {
                                j = ja[k];
                                if (i!=j) t -= aj[k]*uval[j];
                                else d = aj[k];
                            } // end for k
                            if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                        }
                    }
                }
            } // end for i
            else {
#endif
                for (i = 0; i < nrow; i ++) {
                    if (mark[i] != 1) {
                        t = bval[i];
                        begin_row = ia[i], end_row = ia[i+1];
                        for (k = begin_row; k < end_row; k ++) {
                            j = ja[k];
                            if (i!=j) t -= aj[k]*uval[j];
                            else d = aj[k];
                        } // end for k
                        if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
                    }
                } // end for i
#ifdef _OPENMP
            }
#endif
        } // end while
    }
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_ilu (dCSRmat *A, dvector *b, dvector *x, void *data)
 *
 * \brief ILU method as a smoother
 *
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param x      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param data   Pointer to user defined data
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date   2010/11/12
 */
void fasp_smoother_dcsr_ilu (dCSRmat *A,
                             dvector *b,
                             dvector *x,
                             void    *data)
{
    const INT m=A->row, m2=2*m, memneed=3*m;
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
        fasp_darray_cp(m,bval,zr); fasp_blas_dcsr_aAxpy(-1.0,A,xval,zr);
        
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
        
        fasp_blas_darray_axpy(m,1,z,xval);
    }
    
    return;
    
MEMERR:
    printf("### ERROR: ILU needs %d memory, only %d available! %s : %d\n",
           memneed, iludata->nwork, __FILE__, __LINE__);
    fasp_chkerr(ERROR_ALLOC_MEM, __FUNCTION__);
}

/**
 * \fn void fasp_smoother_dcsr_kaczmarz (dvector *u, const INT i_1, const INT i_n,
 *                                       const INT s, dCSRmat *A, dvector *b,
 *                                       INT L, const REAL w)
 *
 * \brief Kaczmarz method as a smoother
 *
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 * \param w      Over-relaxation weight
 *
 * \author Xiaozhe Hu
 * \date   2010/11/12
 *
 * Modified by Chunsheng Feng, Zheng Li on 2012/09/01
 */
void fasp_smoother_dcsr_kaczmarz (dvector    *u,
                                  const INT   i_1,
                                  const INT   i_n,
                                  const INT   s,
                                  dCSRmat    *A,
                                  dvector    *b,
                                  INT         L,
                                  const REAL  w)
{
    const INT   *ia=A->IA,*ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    REAL  temp1,temp2,alpha;
    
#ifdef _OPENMP
    // variables for OpenMP
    const INT    N = ABS(i_n - i_1)+1;
    INT   myid, mybegin, myend;
    INT   nthreads = fasp_get_num_threads();
#endif
    
    if (s > 0) {
        
        while (L--) {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, temp1, temp2, begin_row, end_row, k, alpha, j)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin += i_1, myend += i_1;
                    for (i=mybegin; i<myend; i+=s) {
                        temp1 = 0; temp2 = 0;
                        begin_row=ia[i], end_row=ia[i+1];
                        for (k=begin_row; k<end_row; k++) {
                            j=ja[k];
                            temp1 += aj[k]*aj[k];
                            temp2 += aj[k]*uval[j];
                        } // end for k
                    }
                    alpha = (bval[i] - temp2)/temp1;
                    for (k=begin_row; k<end_row; ++k){
                        j = ja[k];
                        uval[j] += w*alpha*aj[k];
                    }// end for k
                } // end for i
            }
            else {
#endif
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
#ifdef _OPENMP
            }
#endif
        } // end while
        
    } // if s
    
    else {
        while (L--) {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, temp1, temp2, begin_row, end_row, k, alpha, j)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin = i_1 - mybegin, myend = i_1 - myend;
                    for (i=mybegin; i>myend; i+=s) {
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
                }
            }
            else {
#endif
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
#ifdef _OPENMP
            }
#endif
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
 * \param u      Pointer to dvector: the unknowns (IN: initial, OUT: approximation)
 * \param i_1    Starting index
 * \param i_n    Ending index
 * \param s      Increasing step
 * \param A      Pointer to dBSRmat: the coefficient matrix
 * \param b      Pointer to dvector: the right hand side
 * \param L      Number of iterations
 *
 * \author Xiaozhe Hu, James Brannick
 * \date   01/26/2011
 *
 * Modified by Chunsheng Feng, Zheng Li on 09/01/2012
 */
void fasp_smoother_dcsr_L1diag (dvector    *u,
                                const INT   i_1,
                                const INT   i_n,
                                const INT   s,
                                dCSRmat    *A,
                                dvector    *b,
                                INT         L)
{
    const INT    N = ABS(i_n - i_1)+1;
    const INT   *ia=A->IA, *ja=A->JA;
    const REAL  *aj=A->val,*bval=b->val;
    REAL        *uval=u->val;
    
    // local variables
    INT   i,j,k,begin_row,end_row;
    
#ifdef _OPENMP
    // variables for OpenMP
    INT   myid, mybegin, myend;
    INT   nthreads = fasp_get_num_threads();
#endif
    
    // Checks should be outside of for; t,d can be allocated before calling!!! --Chensong
    REAL *t = (REAL *)fasp_mem_calloc(N,sizeof(REAL));
    REAL *d = (REAL *)fasp_mem_calloc(N,sizeof(REAL));
    
    while (L--) {
        if (s>0) {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, begin_row, end_row, k, j)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin += i_1, myend += i_1;
                    for (i=mybegin; i<myend; i+=s) {
                        t[i]=bval[i]; d[i]=0.0;
                        begin_row=ia[i], end_row=ia[i+1];
                        for (k=begin_row; k<end_row; k++) {
                            j=ja[k];
                            t[i]-=aj[k]*uval[j];
                            d[i]+=ABS(aj[k]);
                        }
                    }
                }
#pragma omp parallel for private(i)
                for (i=i_1;i<=i_n;i+=s) {
                    if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
                }
            }
            else {
#endif
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
#ifdef _OPENMP
            }
#endif
        }
        else {
#ifdef _OPENMP
            if (N > OPENMP_HOLDS) {
#pragma omp parallel for private(myid, mybegin, myend, i, k, j, begin_row, end_row)
                for (myid=0; myid<nthreads; myid++) {
                    fasp_get_start_end(myid, nthreads, N, &mybegin, &myend);
                    mybegin = i_1 - mybegin, myend = i_1 - myend;
                    for (i=mybegin; i>myend; i+=s) {
                        t[i]=bval[i]; d[i]=0.0;
                        begin_row=ia[i], end_row=ia[i+1];
                        for (k=begin_row; k<end_row; k++) {
                            j=ja[k];
                            t[i]-=aj[k]*uval[j];
                            d[i]+=ABS(aj[k]);
                        }
                    }
                }
#pragma omp parallel for private(i)
                for (i=i_1;i>=i_n;i+=s) {
                    if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
                }
            }
            else {
#endif
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
#ifdef _OPENMP
            }
#endif
        }
        
    } // end while
    
    fasp_mem_free(t);
    fasp_mem_free(d);
    
    return;
}

/**
 * \fn void fasp_smoother_dcsr_gs_rb3d (dvector *u, dCSRmat *A, dvector *b, INT L, 
 *                                      const INT order, INT *mark, const INT maximap,
 *                                      const INT nx, const INT ny, const INT nz)
 *
 * \brief       Colored Gauss-Seidel smoother for Au=b
 *
 * \param u        Initial guess and the new approximation to the solution
 * \param A        Pointer to stiffness matrix
 * \param b        Pointer to right hand side
 * \param L        Number of iterations
 * \param order    Ordering: -1: Forward; 1: Backward
 * \param mark     Marker for C/F points
 * \param maximap  Size of IMAP
 * \param nx       Number vertex of X direction
 * \param ny       Number vertex of Y direction
 * \param nz       Number vertex of Z direction
 *
 * \author Chunsheng Feng
 * \date   02/08/2012
 */
void fasp_smoother_dcsr_gs_rb3d (dvector    *u,
                                 dCSRmat    *A,
                                 dvector    *b,
                                 INT         L,
                                 const INT   order,
                                 INT        *mark,
                                 const INT   maximap,
                                 const INT   nx,
                                 const INT   ny,
                                 const INT   nz)
{
    const INT   nrow = b->row; // number of rows
    INT        *ia=A->IA,*ja=A->JA;
    REAL       *aa=A->val, *bval=b->val;
    REAL       *uval=u->val;
    
    INT         i,j,k,begin_row,end_row;
    REAL        t,d=0.0;
    
    // forward
    if (order == 1) {
        while (L--) {
            rb0f3d( ia, ja, aa, uval, bval, mark, nx,  ny,  nz, 1);
#if 1 // TODO: Check! Why? --Chensong
            for ( i = maximap; i < nrow; i++ ) {
                t = bval[i];
                begin_row = ia[i], end_row = ia[i+1];
                for ( k = begin_row; k < end_row; k++ ) {
                    j = ja[k];
                    if (i!=j) t -= aa[k]*uval[j];
                    else d = aa[k];
                } // end for k
                if (ABS(d) > SMALLREAL) uval[i] = t/d;
            } // end for i
#endif
        } // end while
    }
    
    // backward
    else {
        while (L--) {
#if 1 // TODO: Check! Why? --Chensong
            for ( i = nrow-1; i >= maximap; i-- ) {
                t = bval[i];
                begin_row = ia[i],end_row = ia[i+1];
                for ( k = begin_row; k < end_row; k ++ ) {
                    j = ja[k];
                    if (i!=j) t -= aa[k]*uval[j];
                    else d = aa[k];
                } // end for k
                if (ABS(d) > SMALLREAL) uval[i] = t/d;
            } // end for i
#endif
            rb0b3d( ia, ja, aa, uval, bval, mark, nx,  ny,  nz, 1);
        } // end while
    }
    return;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void swep2db (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
 *                          INT nbegy, INT *mark, INT nx, INT ny)
 *
 * \brief      Gauss-Seidel backward smoother for certain color
 *
 * \param ia     Pointer to start location of each row
 * \param ja     Pointer to column index of nonzero elements
 * \param aa     Pointer to nonzero elements of
 * \param u      Pointer to initial guess
 * \param f      Pointer to right hand
 * \param nbegx  The stride between the same color nodes in x direction
 * \param nbegy  The stride between the same color nodes in y direction
 * \param mark   Pointer to order of nodes
 * \param nx     Number of nodes in x direction
 * \param ny     Number of nodes in y direction
 *
 * \author  Chunsheng Feng, Zheng Li
 * \date    02/06/2012
 *
 * \note The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void swep2db (INT   *ia,
                     INT   *ja,
                     REAL  *aa,
                     REAL  *u,
                     REAL  *f,
                     INT    nbegx,
                     INT    nbegy,
                     INT   *mark,
                     INT    nx,
                     INT    ny)
{
    INT j, j0, i, i0;
    INT begin_row, end_row, ii, jj;
    REAL t, d;
    
    nbegx = nx + nbegx;
    nbegy = ny + nbegy;
    
#ifdef _OPENMP
#pragma omp parallel for private(j,j0,i,i0,t,begin_row,end_row,ii,jj,d)
#endif
    for (j = nbegy; j >=0; j-=2) {
        j0= j*nx;
        for (i = nbegx; i >=0; i-=2) {
            i0 = i + j0;
            i0 = mark[i0]-1; // Fortran to C
            if (i0>=0 ) {
                t = f[i0];
                begin_row = ia[i0], end_row = ia[i0+1];
                for (ii = begin_row; ii < end_row; ii ++) {
                    jj = ja[ii];
                    if (i0!=jj) t -= aa[ii]*u[jj];
                    else d = aa[ii];
                } // end for ii
                if (ABS(d) > SMALLREAL) u[i0] = t/d;
            } //if (i0>=0 )
        }
    }
}


/**
 * \fn static void swep3db (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
 *                          INT nbegy, INT nbegz, INT *mark, INT nx, INT ny, INT nz)
 *
 * \brief      Gauss-Seidel backward smoother for certain color
 *
 * \param ia     Pointer to start location of each row
 * \param ja     Pointer to column index of nonzero elements
 * \param aa     Pointer to nonzero elements of
 * \param u      Pointer to initial guess
 * \param f      Pointer to right hand
 * \param nbegx  The stride between the same color nodes in x direction
 * \param nbegy  The stride between the same color nodes in y direction
 * \param nbegz  The stride between the same color nodes in z direction
 * \param mark   Pointer to order of nodes
 * \param nx     Number of nodes in x direction
 * \param ny     Number of nodes in y direction
 * \param nz     Number of nodes in z direction
 *
 * \author  Chunsheng Feng, Zheng Li
 * \date    02/06/2012
 *
 * \note The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void swep3db (INT   *ia,
                     INT   *ja,
                     REAL  *aa,
                     REAL  *u,
                     REAL  *f,
                     INT    nbegx,
                     INT    nbegy,
                     INT    nbegz,
                     INT   *mark,
                     INT    nx,
                     INT    ny,
                     INT    nz)
{
    INT nxy, k, k0, j, j0, i, i0;
    INT begin_row, end_row, ii, jj;
    REAL t, d;
    
    nxy=nx*ny;
    nbegx = nx + nbegx;
    nbegy = ny + nbegy;
    nbegz = nz + nbegz;
    
#ifdef _OPENMP
#pragma omp parallel for private(k,k0,j,j0,i,i0,t,begin_row,end_row,ii,jj,d)
#endif
    for (k=nbegz; k >=0; k-=2) {
        k0= k*nxy;
        for (j = nbegy; j >=0; j-=2) {
            j0= j*nx;
            for (i = nbegx; i >=0; i-=2) {
                i0 = i   +  j0    + k0;
                i0 = mark[i0]-1;  // Fortran to C
                if (i0>=0 ) {
                    t = f[i0];
                    begin_row = ia[i0], end_row = ia[i0+1];
                    for (ii = begin_row; ii < end_row; ii ++) {
                        jj = ja[ii];
                        if (i0!=jj) t -= aa[ii]*u[jj];
                        else d = aa[ii];
                    } // end for ii
                    
                    if (ABS(d) > SMALLREAL) u[i0] = t/d;
                } //if (i0>=0 )
            }
        }
    }
}

/*
 * \fn static void rb0b2d (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f,
 *                         INT *mark, INT nx, INT ny, INT nsweeps)
 *
 * \brief Colored Gauss-Seidel backward smoother for Au=b
 *
 * \param ia       Pointer to start location of each row
 * \param ja       Pointer to column index of nonzero elements
 * \param aa       Pointer to nonzero elements of
 * \param u        Pointer to initial guess
 * \param f        Pointer to right hand
 * \param mark     Pointer to order of nodes
 * \param nx       Number of nodes in x direction
 * \param ny       Number of nodes in y direction
 * \param nsweeps  Number of relaxation sweeps
 *
 * \author  Chunsheng Feng, Zheng Li
 * \date    02/06/2012
 *
 * \note This subroutine is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void rb0b2d (INT   *ia,
                    INT   *ja,
                    REAL  *aa,
                    REAL  *u,
                    REAL  *f,
                    INT   *mark,
                    INT    nx,
                    INT    ny,
                    INT    nsweeps)
{
#if 1
    INT n0e,n0o,isweep;
    
    n0e=0;
    n0o=1;
    
    for (isweep = 1; isweep <= nsweeps; isweep++) {
        /*...  e-e */
        swep2df(ia,ja,aa,u,f,n0e,n0e,mark,nx,ny);
        /*...  e-o */
        swep2df(ia,ja,aa,u,f,n0e,n0o,mark,nx,ny);
        /*...  o-e */
        swep2df(ia,ja,aa,u,f,n0o,n0e,mark,nx,ny);
        /*...  o-o */
        swep2df(ia,ja,aa,u,f,n0o,n0o,mark,nx,ny);
    }
#else
    INT n0e = -1, n0o = -2, isweep;
    //INT ex, ey, ez;
    //INT ox, oy, oz;
    
    for (isweep = 1; isweep <= nsweeps; isweep++) {
		if ((nx%2==0) &&(ny%2 ==0)) {
			/*...  e-e (and going backwards) */
			swep2db(ia,ja,aa,u,f,n0e,n0e,mark,nx,ny);
			/*...  e-o */
			swep2db(ia,ja,aa,u,f,n0e,n0o,mark,nx,ny);
			/*...  o-e */
			swep2db(ia,ja,aa,u,f,n0o,n0e,mark,nx,ny);
			/*...  o-o */
			swep2db(ia,ja,aa,u,f,n0o,n0o,mark,nx,ny);
		}
		else if ((nx%2==0)&&(ny%2 ==1)) {
			/*...  e-o (and going backwards) */
			swep2db(ia,ja,aa,u,f,n0e,n0o,mark,nx,ny);
			/*...  e-e */
			swep2db(ia,ja,aa,u,f,n0e,n0e,mark,nx,ny);
			/*...  o-o */
			swep2db(ia,ja,aa,u,f,n0o,n0o,mark,nx,ny);
			/*...  o-e */
			swep2db(ia,ja,aa,u,f,n0o,n0e,mark,nx,ny);
        }
		else if ((nx%2==1)&&(ny%2 ==0)) {
			/*...  o-e (and going backwards) */
			swep2db(ia,ja,aa,u,f,n0o,n0e,mark,nx,ny);
			/*...  o-o */
			swep2db(ia,ja,aa,u,f,n0o,n0o,mark,nx,ny);
			/*...  e-e */
			swep2db(ia,ja,aa,u,f,n0e,n0e,mark,nx,ny);
			/*...  e-o */
			swep2db(ia,ja,aa,u,f,n0e,n0o,mark,nx,ny);
		}
		else if ((nx%2==1)&&(ny%2 ==1)) {
			/*...  o-o (and going backwards) */
			swep2db(ia,ja,aa,u,f,n0o,n0o,mark,nx,ny);
			/*...  o-e */
			swep2db(ia,ja,aa,u,f,n0o,n0e,mark,nx,ny);
			/*...  e-o */
			swep2db(ia,ja,aa,u,f,n0e,n0o,mark,nx,ny);
			/*...  e-e */
			swep2db(ia,ja,aa,u,f,n0e,n0e,mark,nx,ny);
		}
	}
#endif
}

/*
 * \fn static void rb0b3d (INT *ia, INT *ja, REAL *aa,REAL *u, REAL *f,
 *                         INT *mark, INT nx, INT ny, INT nz, INT nsweeps)
 *
 * \brief  Colores Gauss-Seidel backward smoother for Au=b
 *
 * \param ia       Pointer to start location of each row
 * \param ja       Pointer to column index of nonzero elements
 * \param aa       Pointer to nonzero elements of
 * \param u        Pointer to initial guess
 * \param f        Pointer to right hand
 * \param mark     Pointer to order of nodes
 * \param nx       Number of nodes in x direction
 * \param ny       Number of nodes in y direction
 * \param nz       Number of nodes in z direction
 * \param nsweeps  Number of relaxation sweeps
 *
 * \author  Chunsheng Feng
 * \date    02/06/2012
 *
 * \note This subroutine is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void rb0b3d (INT   *ia,
                    INT   *ja,
                    REAL  *aa,
                    REAL  *u,
                    REAL  *f,
                    INT   *mark,
                    INT    nx,
                    INT    ny,
                    INT    nz,
                    INT    nsweeps)
{
#if 1
    INT n0e,n0o,isweep;
    
    n0e=0;
    n0o=1;
    
    for (isweep = 1; isweep <= nsweeps; isweep++) {
        /*...  e-e-e */
        swep3df(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
        /*...  e-e-o */
        swep3df(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
        /*...  e-o-e */
        swep3df(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
        /*...  e-o-o */
        swep3df(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
        /*...  o-e-e */
        swep3df(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
        /*...  o-e-o */
        swep3df(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
        /*...  o-o-e */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
        /*...  o-o-o */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
    }
#else
    INT n0e = -1, n0o = -2, isweep;
    //INT ex, ey, ez;
    //INT ox, oy, oz;
    
    for (isweep = 1; isweep <= nsweeps; isweep++) {
		if ((nx%2==0) &&(ny%2 ==0)  &&(nz%2==0)) {
			/*...  e-e-e (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  e-e-e (and going backwards) */
			//swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
		}
		else if ((nx%2==0) &&(ny%2 ==0)  &&(nz%2==1)) {
			/*...  e-e-o (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  e-e-o (and going backwards) */
			//swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
		}
		else if ((nx%2==0)&&(ny%2 ==1)&&(nz%2==0)) {
			/*...  e-o-e (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  e-o-e (and going backwards) */
            //	swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
		}
		else if ((nx%2==0)&&(ny%2 ==1)&&(nz%2==1)) {
			/*...  e-o-o (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  e-o-o (and going backwards) */
            //	swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
		}
		else if ((nx%2==1)&&(ny%2 ==0)&&(nz%2==0)) {
			/*...  o-e-e (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  o-e-e (and going backwards) */
            //	swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
		}
		else if ((nx%2==1)&&(ny%2 ==0)&&(nz%2==1)) {
			/*...  o-e-o (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  o-e-o (and going backwards) */
            //	swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
		}
		else if ((nx%2==1)&&(ny%2 ==1)&&(nz%2==0)) {
			/*...  o-o-e (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  o-o-o */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  o-o-e (and going backwards) */
			//swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
		}
		else if ((nx%2==1)&&(ny%2 ==1)&&(nz%2==1)) {
			/*...  o-o-o (and going backwards) */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
			/*...  o-o-e */
			swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
			/*...  o-e-o */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
			/*...  o-e-e */
			swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
			/*...  e-o-o */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
			/*...  e-o-e */
			swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
			/*...  e-e-o */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
			/*...  e-e-e */
			swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
			/*...  o-o-o (and going backwards) */
			//swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            
		}
	}
#endif
}

/**
 * \fn static void swep2df (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
 *                          INT nbegy, INT *mark, INT nx, INT ny)
 *
 * \brief      Gauss-Seidel forward smoother for certain color
 *
 * \param ia     Pointer to start location of each row
 * \param ja     Pointer to column index of nonzero elements
 * \param aa     Pointer to nonzero elements of
 * \param u      Pointer to initial guess
 * \param f      Pointer to right hand
 * \param nbegx  The stride between the same color nodes in x direction
 * \param nbegy  The stride between the same color nodes in y direction
 * \param mark   Pointer to order of nodes
 * \param nx     Number of nodes in x direction
 * \param ny     Number of nodes in y direction
 *
 * \author  Chunsheng Feng, Zheng Li
 * \date    02/06/2012
 *
 * \note The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void swep2df (INT   *ia,
                     INT   *ja,
                     REAL  *aa,
                     REAL  *u,
                     REAL  *f,
                     INT    nbegx,
                     INT    nbegy,
                     INT   *mark,
                     INT    nx,
                     INT    ny)
{
    INT j,j0,i,i0;
    INT begin_row,end_row,ii,jj;
    REAL t,d=0;
    
#ifdef _OPENMP
#pragma omp parallel for private(j,j0,i,i0,t,begin_row,end_row,ii,jj,d)
#endif
    for (j = nbegy; j < ny; j+=2) {
        j0= j*nx;
        for (i = nbegx; i < nx; i+=2)    /*!*/
        {
            i0 = i + j0;
            i0 = mark[i0]-1; //Fortran to C
            if (i0>=0 ) {
                t = f[i0];
                begin_row = ia[i0], end_row = ia[i0+1];
                for (ii = begin_row; ii < end_row; ii ++) {
                    jj = ja[ii];
                    if (i0!=jj) t -= aa[ii]*u[jj];
                    else d = aa[ii];
                } // end for ii
                if (ABS(d) > SMALLREAL) u[i0] = t/d;
            } //    if (i0>=0 )
        }
    }
    
}


/**
 * \fn static void swep3df (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f, INT nbegx,
 *                          INT nbegy, INT nbegz, INT *mark, INT nx, INT ny, INT nz)
 *
 * \brief      Gauss-Seidel forward smoother for certain color
 *
 * \param ia     Pointer to start location of each row
 * \param ja     Pointer to column index of nonzero elements
 * \param aa     Pointer to nonzero elements of
 * \param u      Pointer to initial guess
 * \param f      Pointer to right hand
 * \param nbegx  The stride between the same color nodes in x direction
 * \param nbegy  The stride between the same color nodes in y direction
 * \param nbegz  The stride between the same color nodes in z direction
 * \param mark   Pointer to order of nodes
 * \param nx     Number of nodes in x direction
 * \param ny     Number of nodes in y direction
 * \param nz     Number of nodes in z direction
 *
 * \author  Chunsheng Feng, Zheng Li
 * \date    02/06/2012
 *
 * Note: The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void swep3df (INT   *ia,
                     INT   *ja,
                     REAL  *aa,
                     REAL  *u,
                     REAL  *f,
                     INT    nbegx,
                     INT    nbegy,
                     INT    nbegz,
                     INT   *mark,
                     INT    nx,
                     INT    ny,
                     INT    nz)
{
    INT nxy=nx*ny,k,k0,j,j0,i,i0;
    INT begin_row,end_row,ii,jj;
    REAL t,d=0;
        
#ifdef _OPENMP
#pragma omp parallel for private(k,k0,j,j0,i,i0,t,begin_row,end_row,ii,jj,d)
#endif
    for (k=nbegz; k < nz; k+=2) {
        k0= k*nxy;
        for (j = nbegy; j < ny; j+=2) {
            j0= j*nx;
            for (i = nbegx; i < nx; i+=2)    /*!*/
            {
                i0 = i   +  j0    + k0;
                i0 = mark[i0]-1; //Fortran to C

                if (i0>=0 ) {
                    t = f[i0];
                    begin_row = ia[i0], end_row = ia[i0+1];
                    for (ii = begin_row; ii < end_row; ii ++) {
                        jj = ja[ii];
                        if (i0!=jj) t -= aa[ii]*u[jj];
                        else d = aa[ii];
                    } // end for ii
                    
                    if (ABS(d) > SMALLREAL) u[i0] = t/d;
                } //    if (i0>=0 )
            }
        }
    }
    
}

/*
 * \fn static void rb0f2d (INT *ia, INT *ja, REAL *aa, REAL *u, REAL *f,
 *                         INT *mark, INT nx, INT ny, INT nsweeps)
 *
 * \brief  Colores Gauss-Seidel forward smoother for Au = b
 *
 * \param ia       Pointer to the start location to of row
 * \param ja       Pointer to the column index of nonzero elements
 * \param aa       Pointer to the values of the nonzero elements
 * \param u        Pointer to initial value
 * \param f        Pointer to right hand
 * \param mark     Pointer to the order index of nodes
 * \param nx       Number of nodes in x direction
 * \param ny       Number of nodes in y direction
 * \param nsweeps  Number of relaxation sweeps
 *
 * \author  Chunsheng Feng, Zheng Li
 * \data    02/06/2012
 *
 * NOTE: The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void rb0f2d (INT   *ia,
                    INT   *ja,
                    REAL  *aa,
                    REAL  *u,
                    REAL  *f,
                    INT   *mark,
                    INT    nx,
                    INT    ny,
                    INT    nsweeps)
{
    INT n0e,n0o,isweep;
    
    n0e=0;
    n0o=1;
    
    for (isweep = 1; isweep <= nsweeps; isweep++) {
        /*...  o-o */
        swep2df(ia,ja,aa,u,f,n0o,n0o,mark,nx,ny);
        /*...  o-e */
        swep2df(ia,ja,aa,u,f,n0o,n0e,mark,nx,ny);
        /*...  e-o */
        swep2df(ia,ja,aa,u,f,n0e,n0o,mark,nx,ny);
        /*...  e-e */
        swep2df(ia,ja,aa,u,f,n0e,n0e,mark,nx,ny);
    }
}

/*
 * \fn static void rb0f3d (INT *ia, INT *ja, REAL *aa,REAL *u, REAL *f, INT *mark,
 *                         INT nx, INT ny, INT nz, INT nsweeps)
 *
 * \brief  Colores Gauss-Seidel forward smoother for Au=b
 *
 * \param ia       Pointer to the start location to of row
 * \param ja       Pointer to the column index of nonzero elements
 * \param aa       Pointer to the values of the nonzero elements
 * \param u        Pointer to initial value
 * \param f        Pointer to right hand
 * \param mark     Pointer to the order index of nodes
 * \param nx       Number of nodes in x direction
 * \param ny       Number of nodes in y direction
 * \param nz       Number of nodes in z direction
 * \param nsweeps  Number of relaxation sweeps
 *
 * \author      Chunsheng Feng, Zheng Li
 * \data        02/06/2012
 *
 * NOTE: The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 * (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */
static void rb0f3d (INT   *ia,
                    INT   *ja,
                    REAL  *aa,
                    REAL  *u,
                    REAL  *f,
                    INT   *mark,
                    INT    nx,
                    INT    ny,
                    INT    nz,
                    INT    nsweeps)
{
    INT n0e,n0o,isweep;
    
    n0e=0;
    n0o=1;
    
    for (isweep = 1; isweep <= nsweeps; isweep++) {
        /*...  o-o-o */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
        /*...  o-o-e */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
        /*...  o-e-o */
        swep3df(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
        /*...  o-e-e */
        swep3df(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
        /*...  e-o-o */
        swep3df(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
        /*...  e-o-e */
        swep3df(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
        /*...  e-e-o */
        swep3df(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
        /*...  e-e-e */
        swep3df(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
        /*...  o-o-o */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
        
    }
}

#if 0
/**
 * \fn static dCSRmat form_contractor (dCSRmat *A, const INT smoother, const INT steps,
 *                                     const INT ndeg, const REAL relax, const REAL dtol)
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
 * \note This is NOT an O(N) algorithm, need to be modified!!!!
 */
static dCSRmat form_contractor (dCSRmat    *A,
                                const INT   smoother,
                                const INT   steps,
                                const INT   ndeg,
                                const REAL  relax,
                                const REAL  dtol)
{
    const INT   n=A->row;
    INT         i;
    
    REAL *work = (REAL *)fasp_mem_calloc(2*n,sizeof(REAL));
    
    dvector b, x;
    b.row=x.row=n;
    b.val=work; x.val=work+n;
    
    INT *index = (INT *)fasp_mem_calloc(n,sizeof(INT));
    
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
                printf("### ERROR: Wrong smoother type!\n");
                fasp_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
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
