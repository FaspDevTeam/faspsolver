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
 * \date   05/23/2012    Modified by Chunsheng Feng && Xiaoqiang Yue
 *
 * \note x is reused to store the resulting array.
 */
void fasp_blas_array_ax(const INT n, 
                        const REAL a, 
                        REAL *x)
{
	unsigned INT i;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
	    use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}

    if (a == 1.0) {
    }
    else {
        if (use_openmp) {
#pragma omp parallel for private(i) schedule(static)  
            for (i=0; i<n; ++i) x[i] *= a;
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
 * \date   05/23/2012    Modified by Chunsheng Feng && Xiaoqiang Yue
 *
 * \note y is reused to store the resulting array.
 */
void fasp_blas_array_axpy (const INT n, 
                           const REAL a, 
                           REAL *x, 
                           REAL *y)
{
	INT i, nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    if (a==1.0) {
        if (use_openmp) {
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                for (i=mybegin; i<myend; ++i) y[i] += x[i];
            }
        }
        else {
            for (i=0; i<n; ++i) y[i] += x[i];
        }
    }
    else if (a==-1.0) {
        if (use_openmp) {
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                for (i=mybegin; i<myend; ++i) y[i] -= x[i];
            }
        }
        else {
            for (i=0; i<n; ++i) y[i] -= x[i];
        }
    }
    else {
        if (use_openmp) {
            INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
            {
                myid = omp_get_thread_num();
                FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
                for (i=mybegin; i<myend; ++i) y[i] += a*x[i];
            }
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
 * \date   05/23/2012    Modified by Chunsheng Feng Xiaoqiang Yue
 */

void fasp_blas_array_axpyz (const INT n, 
                            const REAL a, 
                            REAL *x, 
                            REAL *y, 
                            REAL *z)
{
	unsigned INT i;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}

    if (use_openmp) {
        INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
            for (i=mybegin; i<myend; ++i) z[i] = a*x[i]+y[i];
        }
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
 * \date   05/23/2012    Modified by Chunsheng Feng Xiaoqiang Yue
 * 
 * \note y is reused to store the resulting array.
 */
void fasp_blas_array_axpby(const INT n, 
                           const REAL a, 
                           REAL *x, 
                           const REAL b, 
                           REAL *y) 
{
	unsigned INT i;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}

    if (use_openmp) {
        INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
            for (i=mybegin; i<myend; ++i) y[i] = a*x[i]+b*y[i];
        }
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
 * \date   05/23/2012    Modified by Chunsheng Feng Xiaoqiang Yue
 */
REAL fasp_blas_array_dotprod(const INT n, 
                             REAL *x, 
                             REAL *y)
{
	unsigned INT i;
    REAL value=0.0;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    if (use_openmp) {
#pragma omp parallel for reduction(+:value) private(i)
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
 * \date   05/23/2012    Modified by Chunsheng Feng Xiaoqiang Yue
 */

REAL fasp_blas_array_norm1 (const INT n, 
                            REAL *x)
{
	unsigned INT i;
    REAL onenorm = 0.;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}

    if (use_openmp) {
#pragma omp parallel for reduction(+:onenorm) private(i)
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
 * \date   05/23/2012    Modified by Chunsheng Feng Xiaoqiang Yue
 */

REAL fasp_blas_array_norm2 (const INT n, 
                            REAL *x) 
{
	unsigned INT i;
    REAL twonorm = 0.;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || n <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}

    if (use_openmp) {
#pragma omp parallel for reduction(+:twonorm) private(i)
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
 */
REAL fasp_blas_array_norminf (const INT n, 
                              REAL *x)
{
    unsigned INT i;
    REAL infnorm=0.0;
    for (i=0;i<n;++i) infnorm=MAX(infnorm,ABS(x[i]));
    return infnorm;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
