/*! \file threads.c
 *
 *  \brief Get and set number of threads and assigne work load for each thread.
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

#ifdef _OPENMP

INT thread_ini_flag = 0;

/**
 * \fn     INT FASP_GET_NUM_THREADS ()
 *
 * \brief  Get the number of threads for thread related functions.
 *
 * \return The number of threads to run
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Zheng Li
 * \date   June/15/2012
 */
INT FASP_GET_NUM_THREADS ()
{
    static INT nthreads;
    
    if ( thread_ini_flag == 0 ) {
#pragma omp parallel
        nthreads = omp_get_num_threads();
        
        printf("\nFASP is running on %3d thread(s).\n\n",nthreads);
        thread_ini_flag = 1;
    }
    
    return nthreads;
}

/**
 * \fn     INT Fasp_Set_Num_Threads (INT nthreads)
 *
 * \brief  Set the number of threads for thread related functions.
 *
 * \param  nthreads  Desirable number of threads
 *
 * \return The number of threads to run
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Zheng Li
 * \date   June/15/2012
 */
INT Fasp_Set_Num_Threads (INT nthreads)
{
    omp_set_num_threads( nthreads );
    
    return nthreads;
}

#endif

/**
 * \fn    void FASP_GET_START_END (INT procid, INT nprocs, INT n, INT *start, INT *end)
 *
 * \brief Assign Load to each thread.
 *
 * \param procid Index of thread
 * \param nprocs Number of threads
 * \param n      Total workload
 * \param start  Pointer to the begin of each thread in total workload
 * \param end    Pointer to the end of each thread in total workload
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Zheng Li
 * \date   June/25/2012
 */
void FASP_GET_START_END (INT procid,
                         INT nprocs,
                         INT n,
                         INT *start,
                         INT *end)
{
    INT chunk_size = n / nprocs;
    INT mod =  n % nprocs;
    INT start_loc, end_loc;
    
    if ( procid < mod) {
        end_loc = chunk_size + 1;
        start_loc = end_loc * procid;
    }
    else {
        end_loc = chunk_size;
        start_loc = end_loc * procid + mod;
    }
    end_loc = end_loc + start_loc;
    
    *start = start_loc;
    *end = end_loc;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
