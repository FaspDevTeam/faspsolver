/*! \file BlaArray.c
 *
 *  \brief BLAS1 operations for arrays
 *
 *  \note This file contains Level-1 (Bla) functions.
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
 * \fn inline void fasp_blas_array_ax (const INT n, const REAL a, REAL *x)
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 *
 * \note x is reused to store the resulting array.
 */
void fasp_blas_array_ax (const INT    n,
                         const REAL   a,
                         REAL        *x)
{
    SHORT use_openmp = FALSE;
    INT   i;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    if ( a == 1.0 ) {
        // do nothing
    }
    else {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i)
            {
                myid = omp_get_thread_num();
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for ( i = mybegin; i < myend; ++i ) x[i] *= a;
            }
#endif
        }
        else {
            for ( i = 0; i < n; ++i ) x[i] *= a;
        }
    }
}

/**
 * \fn void fasp_blas_array_axpy (const INT n, const REAL a,
 *                                const REAL *x, REAL *y)
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 *
 * \note y is reused to store the resulting array.
 */
void fasp_blas_array_axpy (const INT   n,
                           const REAL  a,
                           const REAL *x,
                           REAL       *y)
{
    SHORT use_openmp = FALSE;
    INT   i;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    if ( a == 1.0 ) {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for ( i = mybegin; i < myend; ++i ) y[i] += x[i];
            }
#endif
        }
        else {
            for ( i = 0; i < n; ++i ) y[i] += x[i];
        }
    }
    
    else if ( a == -1.0 ) {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for ( i = mybegin; i < myend; ++i ) y[i] -= x[i];
            }
#endif
        }
        else {
            for ( i = 0; i < n; ++i ) y[i] -= x[i];
        }
    }
    
    else {
        if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for ( i = mybegin; i < myend; ++i ) y[i] += a*x[i];
            }
#endif
        }
        else {
            for ( i = 0; i < n; ++i ) y[i] += a*x[i];
        }
    }
}

/**
 * \fn void fasp_blas_array_axpyz (const INT n, const REAL a, const REAL *x,
 *                                 const REAL *y, REAL *z)
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
void fasp_blas_array_axpyz (const INT   n,
                            const REAL  a,
                            const REAL *x,
                            const REAL *y,
                            REAL       *z)
{
    SHORT use_openmp = FALSE;
    INT   i;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for ( i = mybegin; i < myend; ++i ) z[i] = a*x[i] + y[i];
        }
#endif
    }
    else {
        for ( i = 0; i < n; ++i ) z[i] = a*x[i] + y[i];
    }
}

/**
 * \fn void fasp_blas_array_axpby (const INT n, const REAL a, const REAL *x,
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 *
 * \note y is reused to store the resulting array.
 */
void fasp_blas_array_axpby (const INT   n,
                            const REAL  a,
                            const REAL *x,
                            const REAL  b,
                            REAL       *y)
{
    SHORT use_openmp = FALSE;
    INT   i;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
            for ( i = mybegin; i < myend; ++i ) y[i] = a*x[i] + b*y[i];
        }
#endif
    }
    else {
        for ( i = 0; i < n; ++i ) y[i] = a*x[i] + b*y[i];
    }
    
}

/**
 * \fn REAL fasp_blas_array_dotprod (const INT n, const REAL *x, const REAL *y)
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
REAL fasp_blas_array_dotprod (const INT    n,
                              const REAL  *x,
                              const REAL  *y)
{
    SHORT use_openmp = FALSE;
    INT   i;
    
    register REAL value = 0.0;
    
#ifdef _OPENMP
    if ( n > OPENMP_HOLDS ) use_openmp = TRUE;
#endif
    
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:value) private(i)
#endif
        for ( i = 0; i < n; ++i ) value += x[i]*y[i];
    }
    else {
        for ( i = 0; i < n; ++i ) value += x[i]*y[i];
    }
    
    return value;
}

/**
 * \fn REAL fasp_blas_array_norm1 (const INT n, const REAL *x)
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
REAL fasp_blas_array_norm1 (const INT    n,
                            const REAL  *x)
{
    INT   i;

    register REAL onenorm = 0.0;
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+:onenorm) private(i)
#endif
    for ( i = 0; i < n; ++i ) onenorm+=ABS(x[i]);

    return onenorm;
}

/**
 * \fn REAL fasp_blas_array_norm2 (const INT n, const REAL *x)
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
REAL fasp_blas_array_norm2 (const INT    n,
                            const REAL  *x)
{
    INT  i;
    register REAL twonorm = 0.;
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+:twonorm) private(i)
#endif
    for ( i = 0; i < n; ++i ) twonorm+=x[i]*x[i];
    
    return sqrt(twonorm);
}

/**
 * \fn REAL fasp_blas_array_norminf (const INT n, const REAL *x)
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
 * Modified by Chunsheng Feng, Zheng Li on 06/28/2012
 */
REAL fasp_blas_array_norminf (const INT    n,
                              const REAL  *x)
{
    SHORT use_openmp = FALSE;
    INT   i;
    REAL  infnorm = 0.0;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, nthreads;
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    if (use_openmp) {
#ifdef _OPENMP
        REAL infnorm_loc = 0.0;
#pragma omp parallel firstprivate(infnorm_loc) private(myid, mybegin, myend, i)
        {
            myid = omp_get_thread_num();
            fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
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
