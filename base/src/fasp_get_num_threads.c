#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

/**
 * \fn     FASP_GET_NUM_THREADS();
 * \brief  set the number of threads 
 * \author Chunsheng Feng && Xiaoqiang Yue 
 * \data   05/24/2012
 */
int FASP_GET_NUM_THREADS()
{
	static int flag, nthreads;
	flag = 0;
	if(flag == 0){
       #pragma omp parallel
		nthreads = omp_get_num_threads();
		flag = 1;
	}
	return nthreads;
}
