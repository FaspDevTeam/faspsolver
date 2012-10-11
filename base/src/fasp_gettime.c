/*! \file fasp_gettime
 */

#include <time.h>
#include <omp.h>

#include "fasp_functs.h"
#include "fasp.h"


/**
 * \fn fasp_gettime(REAL *time)
 * \biref Timing 
 * 
 * \author Chunsheng Feng, Zheng LI
 * \date   11/10/2012
 *
 */

void fasp_gettime(REAL *time)
{
    if (time != NULL) {
#ifdef _OPENMP
        *time = omp_get_wtime();
#else
        *time = (REAL) clock()*1.e-6;
#endif 
    }
}
