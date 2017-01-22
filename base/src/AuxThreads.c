/*! \file AuxThreads.c
 *
 *  \brief Get and set number of threads and assign work load for each thread
 *
 *  \note This file contains Level-0 (Aux) functions
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
 * \fn     INT fasp_get_num_threads ()
 *
 * \brief  Get the number of threads for thread related functions.
 *
 * \return The number of threads to run
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Zheng Li
 * \date   June/15/2012
 */
INT fasp_get_num_threads ()
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
 * \fn     INT fasp_set_num_threads (INT nthreads)
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
INT fasp_set_num_threads (INT nthreads)
{
    omp_set_num_threads( nthreads );
    
    return nthreads;
}

#endif

/**
 * \fn    void fasp_get_start_end (INT procid, INT nprocs, INT n, INT *start, INT *end)
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
void fasp_get_start_end (INT  procid,
                         INT  nprocs,
                         INT  n,
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

INT THDs_AMG_GS=0;  /**< AMG GS smoothing threads      */
INT THDs_CPR_lGS=0; /**< reservoir GS smoothing threads     */
INT THDs_CPR_gGS=0; /**< global matrix GS smoothing threads */

/**
 * \fn void fasp_set_GS_threads (INT threads,INT its)
 *
 * \brief  Set threads for CPR. Please add it at the begin of Krylov OpenMP method 
 *         function and after iter++.
 *
 * \param threads  Total threads of solver
 * \param its      Current its of the Krylov methods
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date   03/20/2011
 */
void fasp_set_GS_threads (INT mythreads,
                          INT its)
{
#ifdef _OPENMP
    
#if 1
    
    if (its <=8) {
        THDs_AMG_GS =  mythreads;
        THDs_CPR_lGS = mythreads ;
        THDs_CPR_gGS = mythreads ;
    }
    else if (its <=12) {
        THDs_AMG_GS =  mythreads;
        THDs_CPR_lGS = (6 < mythreads) ? 6 : mythreads;
        THDs_CPR_gGS = (4 < mythreads) ? 4 : mythreads;
    }
    else if (its <=15) {
        THDs_AMG_GS =  (3 < mythreads) ? 3 : mythreads;
        THDs_CPR_lGS = (3 < mythreads) ? 3 : mythreads;
        THDs_CPR_gGS = (2 < mythreads) ? 2 : mythreads;
    }
    else if (its <=18) {
        THDs_AMG_GS =  (2 < mythreads) ? 2 : mythreads;
        THDs_CPR_lGS = (2 < mythreads) ? 2 : mythreads;
        THDs_CPR_gGS = (1 < mythreads) ? 1 : mythreads;
    }
    else {
        THDs_AMG_GS =  1;
        THDs_CPR_lGS = 1;
        THDs_CPR_gGS = 1;
    }
    
#else
    
    THDs_AMG_GS =  mythreads;
    THDs_CPR_lGS = mythreads ;
    THDs_CPR_gGS = mythreads ;
    
#endif
    
#endif // _OPENMP
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
