#include<stdio.h>
#include<stdlib.h>
#include"fasp.h"
#include<omp.h>

/**
 * \fn     INT FASP_GET_NUM_THREADS();
 *
 * \brief  Get the number of threads for OpenMP functions.
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Zheng Li
 *
 * \data   June/15/2012
 */
 
INT thread_ini_flag = 0;

INT FASP_GET_NUM_THREADS()
{
	static INT nthreads;
	if(thread_ini_flag == 0){
       #pragma omp parallel
		nthreads = omp_get_num_threads();
		printf("\n Fasp running at %3d  threads.\n\n",nthreads);
		thread_ini_flag = 1;
	}
	return nthreads;
}

/**
 * \fn     INT Fasp_Set_Num_Threads();
 *
 * \brief  Set the number of threads for OpenMP functions.
 *
 * \author Chunsheng Feng, Xiaoqiang Yue and Zheng Li
 * \data   June/15/2012
 */

INT Fasp_Set_Num_Threads(INT nthreads )
{
	omp_set_num_threads( nthreads );

	return nthreads;
}


