/*! \file AuxArray.c
 *
 *  \brief Simple array operations -- init, set, copy, etc
 *
 *  \note This file contains Level-0 (Aux) functions
 */

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"

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
void fasp_array_set (const INT    n,
                     REAL        *x,
                     const REAL   val)
{
    INT i;
    const SHORT use_openmp = FALSE;
    
#ifdef _OPENMP 
    INT nthreads = 1;

    if ( n > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
    }
#endif
    
    if (val == 0.0) {
        if (use_openmp) {
#ifdef _OPENMP
            INT mybegin,myend,myid;
#pragma omp parallel for private(myid,mybegin,myend)
            for (myid = 0; myid < nthreads; myid ++) {
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                memset(&x[mybegin], 0x0, sizeof(REAL)*(myend-mybegin));
            }
#endif
        }
        else 
            memset(x, 0x0, sizeof(REAL)*n);
    }
    else {
        if (use_openmp) {
#ifdef _OPENMP
            INT mybegin,myend,myid;
#pragma omp parallel for private(myid,mybegin,myend,i)
            for (myid = 0; myid < nthreads; myid ++) {
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) x[i]=val;
            }
#endif
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
void fasp_iarray_set (const INT   n,
                      INT        *x,
                      const INT   val)
{
    INT i;
    const SHORT use_openmp = FALSE;
    
#ifdef _OPENMP 
    INT nthreads = 1;

    if ( n > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = fasp_get_num_threads();
	}
#endif
    
    if (val == 0) {
        if (use_openmp) {
#ifdef _OPENMP
            INT mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin, myend)
            for (myid = 0; myid < nthreads; myid ++) {
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                memset(&x[mybegin], 0, sizeof(INT)*(myend-mybegin));
            }
#endif
        }
        else {
            memset(x, 0, sizeof(INT)*n);
        }
    }
    else {
        if (use_openmp) {
#ifdef _OPENMP
            INT mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin, myend,i)
            for (myid = 0; myid < nthreads; myid ++) {
                fasp_get_start_end(myid, nthreads, n, &mybegin, &myend);
                for (i=mybegin; i<myend; ++i) x[i]=val;
            }
#endif
        }
        else {
            for (i=0; i<n; ++i) x[i]=val;
        }
    }
}

/**
 * \fn void fasp_array_cp (const INT n, const REAL *x, REAL *y)
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
void fasp_array_cp (const INT    n,
                    const REAL  *x,
                    REAL        *y)
{
    memcpy(y, x, n*sizeof(REAL));
}


/**
 * \fn void fasp_iarray_cp (const INT n, const INT *x, INT *y)
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
void fasp_iarray_cp (const INT   n,
                     const INT  *x,
                     INT        *y)
{
    memcpy(y, x, n*sizeof(INT));
}

/**
 * \fn void fasp_array_cp_nc3 (const REAL *x, REAL *y)
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
void fasp_array_cp_nc3 (const REAL  *x,
                        REAL        *y)
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
}

/**
 * \fn void fasp_array_cp_nc5 (const REAL *x, REAL *y)
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
void fasp_array_cp_nc5 (const REAL  *x,
                        REAL        *y)
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];
    y[4] = x[4];
}

/**
 * \fn void fasp_array_cp_nc7 (const REAL *x, REAL *y)
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
void fasp_array_cp_nc7 (const REAL  *x,
                        REAL        *y)
{
    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    y[3] = x[3];
    y[4] = x[4];
    y[5] = x[5];
    y[6] = x[6];
}

/**
 * \fn void fasp_array_permut_nb (const INT n, const INT nb, const REAL*x, 
 *                                const INT *p, REAL*y)
 *
 * \brief Array permutation
 *
 * \param n    Size of array
 * \param nb   Step size  
 * \param x    Pointer to the original vector 
 * \param p    Pointer to index mapping 
 * \param y    Pointer to the destination vector
 * 
 * \author Zheng Li
 * \date   12/04/2016
 */
void fasp_array_permut_nb (const INT    n,
                           const INT    nb,
                           const REAL  *x,
                           const INT   *p,
                           REAL        *y)
{
   INT i, j, indx, indy;

#ifdef _OPENMP
#pragma omp parallel for private(i, j, indx, indy)
#endif
   for (i=0; i<n; ++i) {
       indx = p[i]*nb;
       indy = i*nb;
       for (j=0; j<nb; ++j) {
           y[indy+j] = x[indx+j];
       }
   }
}

/**
 * \fn void fasp_array_invpermut_nb (const INT n, const INT nb, const REAL*x, 
 *                                   const INT *p, REAL*y)
 *
 * \brief Array inverse permutation
 *
 * \param n    Size of array
 * \param nb   Step size
 * \param x    Pointer to the original vector
 * \param p    Pointer to index mapping
 * \param y    Pointer to the destination vector
 *
 * \author Zheng Li
 * \date   12/04/2016
 */
void fasp_array_invpermut_nb (const INT    n,
                              const INT    nb,
                              const REAL  *x,
                              const INT   *p,
                              REAL        *y)
{
   INT i, j, indx, indy;

#ifdef _OPENMP
#pragma omp parallel for private(i, j, indx, indy)
#endif
   for (i=0; i<n; ++i) {
       indx = i*nb;
       indy = p[i]*nb;
       for (j=0; j<nb; ++j) {
           y[indy+j] = x[indx+j];
       }
   }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
