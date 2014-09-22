/*! \file timing.c
 *
 *  \brief Timing subroutines
 */

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn fasp_gettime (REAL *time)
 *
 * \brief Get system time
 *
 * \author Chunsheng Feng, Zheng LI
 * \date   11/10/2012
 *
 * Modified by Chensong Zhang on 09/22/2014: Use CLOCKS_PER_SEC for cross-platform
 */
void fasp_gettime (REAL *time)
{
    if ( time != NULL ) {
#ifdef _OPENMP
        *time = omp_get_wtime();
#else
        *time = (REAL) clock() / CLOCKS_PER_SEC;
#endif
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
