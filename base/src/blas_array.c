/*! \file blas_array.c
 *  \brief BLAS operations for arrays.
 * 
 *  Some simple array operations ...
 *
 */

#include <math.h>
#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_array_ax (const INT n, const REAL a, REAL *x)
 *
 * \brief x = a*x
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012
 *
 * \note x is reused to store the resulting array.
 */
void fasp_blas_array_ax(const INT n, 
                        const REAL a, 
                        REAL *x)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;
   
#ifdef _OPENMP 
	if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

    if (a == 1.0) {
    
    }
    else {
        if (use_openmp) {
#ifdef _OPENMP 
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) 
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) x[i] *= a;
            }
#endif
        }
        else {
            for (i=0; i<n; ++i) x[i] *= a;
        }
    }
}

/**
 * \fn void fasp_blas_array_axpy (const INT n, const REAL a, REAL *x, REAL *y)
 *
 * \brief y = a*x + y
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 * \param y    Pointer to y
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 *
 * \note y is reused to store the resulting array.
 */
void fasp_blas_array_axpy (const INT n, 
                           const REAL a, 
                           REAL *x, 
                           REAL *y)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP 
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    if (a==1.0) {
        if (use_openmp) {
#ifdef _OPENMP 
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) y[i] += x[i];
            }
#endif
        }
        else {
            for (i=0; i<n; ++i) y[i] += x[i];
        }
    }
    else if (a==-1.0) {
        if (use_openmp) {
#ifdef _OPENMP 
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) y[i] -= x[i];
            }
#endif
        }
        else {
            for (i=0; i<n; ++i) y[i] -= x[i];
        }
    }
    else {
        if (use_openmp) {
#ifdef _OPENMP 
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) y[i] += a*x[i];
            }
#endif
        }
        else {
            for (i=0; i<n; ++i) y[i] += a*x[i];
        }
    }
}

/**
 * \fn void fasp_blas_array_axpyz(const INT n, const REAL a, REAL *x,
 *                                REAL *y, REAL *z)
 *
 * \brief z = a*x + y
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 * \param y    Pointer to y
 * \param z    Pointer to z
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */

void fasp_blas_array_axpyz (const INT n, 
                            const REAL a, 
                            REAL *x, 
                            REAL *y, 
                            REAL *z)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP 
	if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

    if (use_openmp) {
#ifdef _OPENMP 
        INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) z[i] = a*x[i]+y[i];
        }
#endif
    }
    else {
        for (i=0; i<n; ++i) z[i] = a*x[i]+y[i];
    }
}

/**
 * \fn void fasp_blas_array_axpby (const INT n, const REAL a, REAL *x, 
 *                                 const REAL b, REAL *y)
 *
 * \brief y = a*x + b*y
 *
 * \param n    Number of variables
 * \param a    Factor a
 * \param x    Pointer to x
 * \param b    Factor b
 * \param y    Pointer to y
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 * 
 * \note y is reused to store the resulting array.
 */
void fasp_blas_array_axpby(const INT n, 
                           const REAL a, 
                           REAL *x, 
                           const REAL b, 
                           REAL *y) 
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP 
	if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

    if (use_openmp) {
#ifdef _OPENMP 
        INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) y[i] = a*x[i]+b*y[i];
        }
#endif
    }
    else {
        for (i=0; i<n; ++i) y[i] = a*x[i]+b*y[i];
    }

}

/**
 * \fn REAL fasp_blas_array_dotprod (const INT n, REAL *x, REAL *y) 
 *
 * \brief Inner product of two arraies (x,y)
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 * \param y    Pointer to y
 *
 * \return     Inner product (x,y)
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */
REAL fasp_blas_array_dotprod(const INT n, 
                             REAL *x, 
                             REAL *y)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;
    REAL value=0.0;

#ifdef _OPENMP 
	if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif
    
    if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel for reduction(+:value) private(i)
#endif
        for (i=0;i<n;++i) value+=x[i]*y[i];
    }
    else {
        for (i=0;i<n;++i) value+=x[i]*y[i];
    }
    return value;
}

/**
 * \fn REAL fasp_blas_array_norm1 (const INT n, REAL *x)
 *
 * \brief L1 norm of array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     L1 norm of x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */

REAL fasp_blas_array_norm1 (const INT n, 
                            REAL *x)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;
    REAL onenorm = 0.;
 
#ifdef _OPENMP 
	if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

    if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel for reduction(+:onenorm) private(i)
#endif
        for (i=0;i<n;++i) onenorm+=ABS(x[i]);
    }
    else {
        for (i=0;i<n;++i) onenorm+=ABS(x[i]);
    }
    return onenorm;
}

/**
 * \fn REAL fasp_blas_array_norm2 (const INT n, REAL *x) 
 *
 * \brief L2 norm of array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     L2 norm of x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */

REAL fasp_blas_array_norm2 (const INT n, 
                            REAL *x) 
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;
    REAL twonorm = 0.;

#ifdef _OPENMP 
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif

    if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel for reduction(+:twonorm) private(i)
#endif
        for (i=0;i<n;++i) twonorm+=x[i]*x[i];
    }
    else {
        for (i=0;i<n;++i) twonorm+=x[i]*x[i];
    }
    return sqrt(twonorm);
}

/**
 * \fn REAL fasp_blas_array_norminf (const INT n, REAL *x) 
 *
 * \brief Linf norm of array x
 *
 * \param n    Number of variables
 * \param x    Pointer to x
 *
 * \return     L_inf norm of x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Zheng Li
 * \date   06/28/2012
 */

REAL fasp_blas_array_norminf (const INT n, 
                              REAL *x)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;
    REAL infnorm = 0.0;

#ifdef _OPENMP 
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
	
    if (use_openmp) {
#ifdef _OPENMP 
        INT myid, mybegin, myend;
        REAL infnorm_loc = 0.0;
#pragma omp parallel firstprivate(infnorm_loc) private(myid, mybegin, myend, i)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i) infnorm_loc = MAX(infnorm_loc, ABS(x[i]));
                 if (infnorm_loc > infnorm) {
                     #pragma omp critical 
                     infnorm = MAX(infnorm_loc, infnorm);
                 }
        }
#endif
    }
    else {
        for (i=0;i<n;++i) infnorm=MAX(infnorm,ABS(x[i]));
    }

    return infnorm;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
