/*! \file array.c
 *
 *  \brief Simple array operations -- init, set, copy, etc
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
 * \fn void fasp_array_null (REAL *x)
 *
 * \brief Initialize an array
 *
 * \param x    Pointer to the vector
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
 */
void fasp_array_null (REAL *x) 
{
    x = NULL;
}

/**
 * \fn void fasp_array_set (const INT n, REAL *x, const REAL val)
 *
 * \brief Set initial value for an array to be x=val
 *
 * \param n    Number of variables
 * \param x    Pointer to the vector
 * \param val  Initial value for the REAL array
 *
 * \author Chensong Zhang 
 * \date   04/03/2010  
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012
 */
void fasp_array_set (const INT n, 
                     REAL *x, 
                     const REAL val)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE; 
    
#ifdef _OPENMP 
    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    if (val == 0.0) {
        if (use_openmp) {
            INT mybegin,myend,myid;
#ifdef _OPENMP 
#pragma omp parallel for private(myid,mybegin,myend)
#endif
            for (myid = 0; myid < nthreads; myid ++) {
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                memset(&x[mybegin], 0x0, sizeof(REAL)*(myend-mybegin));
            }
        }
        else 
            memset(x, 0x0, sizeof(REAL)*n);
    }
    else {
        if (use_openmp) {
            INT mybegin,myend,myid;
#ifdef _OPENMP 
#pragma omp parallel for private(myid,mybegin,myend,i)
#endif
            for (myid = 0; myid < nthreads; myid ++) {
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) x[i]=val;
            }
        }
        else {
            for (i=0; i<n; ++i) x[i] = val;
        }
    }
}

/**
 * \fn void fasp_iarray_set (const INT n, INT *x, const INT val)
 *
 * \brief Set initial value for an array to be x=val
 *
 * \param n    Number of variables
 * \param x    Pointer to the vector
 * \param val  Initial value for the REAL array
 *
 * \author Chensong Zhang 
 * \date   04/03/2010  
 * 
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/25/2012
 */
void fasp_iarray_set (const INT n,
                      INT *x,
                      const INT val)
{
    INT i;
    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP 
	if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif
    
    if (val == 0) {
        if (use_openmp) {
            INT mybegin,myend,myid;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend)
#endif
            for (myid = 0; myid < nthreads; myid ++) {
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                memset(&x[mybegin], 0, sizeof(INT)*(myend-mybegin));
            }
        }
        else {
            memset(x, 0, sizeof(INT)*n);
        }
    }
    else {
        if (use_openmp) {
            INT mybegin,myend,myid;
#ifdef _OPENMP 
#pragma omp parallel for private(myid, mybegin, myend,i)
#endif
            for (myid = 0; myid < nthreads; myid ++) {
                FASP_GET_START_END(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) x[i]=val;
            }
        }
        else {
            for (i=0; i<n; ++i) x[i]=val;
        }
    }
}

/**
 * \fn void fasp_array_cp (const INT n, REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x
 *
 * \param n    Number of variables
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
 */
void fasp_array_cp (const INT n, 
                    REAL *x, 
                    REAL *y)
{
    memcpy(y, x, n*sizeof(REAL));
}


/**
 * \fn void fasp_iarray_cp (const INT n, INT *x, INT *y) 
 *
 * \brief Copy an array to the other y=x
 *
 * \param n    Number of variables
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012  
 */
void fasp_iarray_cp (const INT n, 
                     INT *x, 
                     INT *y)
{
    memcpy(y, x, n*sizeof(INT));
}

/**
 * \fn void fasp_array_cp_nc3 (REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x, the length is 3
 *
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   05/01/2010
 *
 * \note Special unrolled routine designed for a specific application
 */
void fasp_array_cp_nc3 (REAL *x, 
                        REAL *y) 
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
}

/**
 * \fn void fasp_array_cp_nc5(REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x, the length is 5
 *
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   05/01/2010
 *
 * \note Special unrolled routine designed for a specific application
 */
void fasp_array_cp_nc5 (REAL *x, 
                        REAL *y) 
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];
    y[4] = x[4];
}

/**
 * \fn void fasp_array_cp_nc7(REAL *x, REAL *y) 
 *
 * \brief Copy an array to the other y=x, the length is 7
 *
 * \param x    Pointer to the original vector
 * \param y    Pointer to the destination vector
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   05/01/2010
 *
 * \note Special unrolled routine designed for a specific application
 */
void fasp_array_cp_nc7 (REAL *x, 
                        REAL *y) 
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];
    y[4] = x[4];
    y[5] = x[5];
    y[6] = x[6];
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

