/*! \file timing.c
 *  \brief Timing subroutines
 */

#include <time.h>
#include <omp.h>

#include "fasp_functs.h"
#include "fasp.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn fasp_gettime (REAL *time)
 * \biref Simple timing
 *
 * \author Chunsheng Feng, Zheng LI
 * \date   11/10/2012
 *
 */
void fasp_gettime (REAL *time)
{
    if (time != NULL) {
#ifdef _OPENMP
        *time = omp_get_wtime();
#else
        *time = (REAL) clock()*1.e-6;
#endif
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
