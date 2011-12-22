/*! \file vec.c
 *  \brief Simple operations for vectors (INT and REAL). 
 *
 *  \note 
 *  Every structures should be initialized before usage.
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dvector fasp_dvec_create (INT m)
 *
 * \brief Create dvector data space of REAL type
 *
 * \param m    number of rows
 *
 * \return u   the new dvector
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
dvector fasp_dvec_create (INT m)
{		
	dvector u;
	
	u.row = m;
	u.val = (REAL *)fasp_mem_calloc(m,sizeof(REAL));
	
#if CHMEM_MODE		
	total_alloc_mem += m*sizeof(REAL);
#endif
	
	return u;
}

/**
 * \fn ivector fasp_ivec_create (INT m)
 *
 * \brief Create vector data space of INT type
 *
 * \param m   number of rows
 *
 * \return u  the new ivector
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
ivector fasp_ivec_create (INT m)
{		
	ivector u;
	
	u.row = m;
	u.val = (INT *)fasp_mem_calloc(m,sizeof(INT)); 
	
#if CHMEM_MODE		
	total_alloc_mem += m*sizeof(REAL);
#endif
	
	return u;
}

/**
 * \fn void fasp_dvec_alloc (INT m, dvector *u)
 *
 * \brief Create dvector data space of REAL type
 *
 * \param m    integer, number of rows
 * \param *u   pointer to the REAL vector (OUTPUT)
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
void fasp_dvec_alloc (INT m, 
                      dvector *u)
{			
	u->row = m;
	u->val = (REAL*)fasp_mem_calloc(m,sizeof(REAL)); 
	
#if CHMEM_MODE		
	total_alloc_mem += m*sizeof(REAL);
#endif
	
	return;
}

/**
 * \fn void fasp_ivec_alloc (INT m, ivector *u)
 *
 * \brief Create vector data space of INT type
 *
 * \param m   integer, number of rows
 * \param *u  pointer to the integer vector (OUTPUT)
 *
 * \author Chensong Zhang 
 * \date 2010/04/06
 */
void fasp_ivec_alloc (INT m, 
                      ivector *u)
{		
	
	u->row = m;
	u->val = (INT*)fasp_mem_calloc(m,sizeof(INT));
	
#if CHMEM_MODE		
	total_alloc_mem += m*sizeof(INT);
#endif
	
	return;
}

/**
 * \fn void fasp_dvec_free (dvector *u)
 *
 * \brief Free vector data space of REAL type
 *
 * \param *u  pointer to the vector which needs to be deallocated
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_dvec_free (dvector *u)
{		
	if (u==NULL) return;
    
	fasp_mem_free(u->val);
	u->row = 0; u->val = NULL; 
}

/**
 * \fn void fasp_ivec_free (ivector *u)
 *
 * \brief Free vector data space of INT type
 *
 * \param *u  pointer to the vector which needs to be deallocated
 *
 * \note This function is same as fasp_dvec_free except input pointer type.
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_ivec_free (ivector *u)
{		
	if (u==NULL) return;
    
	fasp_mem_free(u->val);
	u->row = 0; u->val = NULL; 
}

/**
 * \fn void fasp_dvec_init (dvector *x) 
 *
 * \brief Initialize dvector
 *
 * \param *x  pointer to dvector which needs to be initialized
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_dvec_init (dvector *x) 
{
	x->row = 0; x->val = NULL;
}

/**
 * \fn void fasp_dvec_rand (INT n, dvector *x)
 *
 * \brief Generate random REAL vector in the range from 0 to 1
 *
 * \param x   pointer to a dvector
 * 
 * \note
 * Sample usage: 
 * \par
 *   dvector xapp;
 * \par
 *   fasp_dvec_create(100,&xapp);
 * \par
 *   fasp_dvec_rand(100,&xapp);
 * \par
 *   fasp_dvec_print(100,&xapp);
 *
 * \author Chensong Zhang
 * \date 11/16/2009
 */
void fasp_dvec_rand (INT n, 
                     dvector *x)
{
    const INT va=(REAL) 0;
    const INT vb=(REAL) n;
    
    unsigned int s=1; srand(s);
	
    INT i,j;

	x->row = n;
    for (i=0; i<n; ++i){
        j = 1 + (INT) (((REAL)n)*rand()/(RAND_MAX+1.0));
        x->val[i] = (((REAL)j)-va)/(vb-va);
    }
}

/**
 * \fn void fasp_dvec_set (INT n, dvector *x, REAL val) 
 *
 * \brief Initialize dvector x[i]=val for i=0:end-1
 *
 * \param   n    number of variables
 * \param * x    pointer to dvector
 * \param val    initial value for the dvector
 *
 * \author Chensong Zhang
 * \date 11/16/2009
 */
void fasp_dvec_set (INT n, 
                    dvector *x, 
                    REAL val)
{
	unsigned INT i;
	REAL *xpt=x->val;
	
	if (n>0) {
        x->row=n; 
    }
	else {
        n=x->row;	   
    }
    
	for (i=0; i<n; ++i) xpt[i]=val;
}

/**
 * \fn void fasp_ivec_set (INT m, ivector *u)
 *
 * \brief Set ivector value to be m
 *
 * \param   m   integer value of ivector
 * \param * u   pointer to the vector (MODIFIED)
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_ivec_set (INT m, ivector *u)
{		
	unsigned INT i;
	for (i=0; i<u->row; ++i) u->val[i] = m;
}

/**
 * \fn void fasp_dvec_cp (dvector *x, dvector *y) 
 *
 * \brief Copy dvector x to dvector y
 *
 * \param *x pointer to dvector
 * \param *y pointer to dvector (MODIFIED)
 *
 * \author Chensong Zhang
 * \date 11/16/2009
 */
void fasp_dvec_cp (dvector *x, 
                   dvector *y)
{
	y->row=x->row;
	memcpy(y->val,x->val,x->row*sizeof(REAL));
}

/**
 * \fn REAL fasp_dvec_maxdiff (dvector *x, dvector *y) 
 *
 * \brief Maximal difference of two dvector x and y
 *
 * \param * x   pointer to dvector
 * \param * y   pointer to dvector
 *
 * \return      maximal norm of x-y
 *
 * \author Chensong Zhang
 * \date 11/16/2009
 */
REAL fasp_dvec_maxdiff (dvector *x, 
                        dvector *y)
{
	const INT length=x->row;
	REAL Linf=0, diffi=0;
	REAL *xpt=x->val, *ypt=y->val;
	
	unsigned INT i;
	for (i=0; i<length; ++i) {
		if ((diffi = ABS(xpt[i]-ypt[i])) > Linf) Linf = diffi;
	}
	return Linf;
}

/**
 * \fn void fasp_dvec_symdiagscale (dvector *b, dvector *diag)
 *
 * \brief Symmetric diagonal scaling D^{-1/2}b
 *
 * \param *b      pointer to the dvector
 * \param *diag   pointer to the diagonal entries
 *
 * \author Xiaozhe Hu
 * \date 01/31/2011
 */
void fasp_dvec_symdiagscale (dvector *b, 
                             dvector *diag)
{
	// information about dvector
	const INT n = b->row;
	REAL *val = b->val;
	
	// local variables
	unsigned INT i;
	
	if (diag->row != n)
	{
		printf("### ERROR: size of diag = %d and size of dvector = %d mismatch!!", 
               diag->row, n);
		exit(ERROR_MISC);
	}
	
	// main loop
	for (i=0; i<n; i++) val[i] = val[i]/sqrt(diag->val[i]);
	
	return;
}


/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_dvec_set_omp (int n, dvector *x, double val, int nthreads, int openmp_holds)
 * \brief Initialize dvector x=val
 *
 * \param n number of variables
 * \param *x pointer to dvector
 * \param val initial value for the dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dvec_set_omp (int n, 
                        dvector *x, 
                        double val, 
                        int nthreads, 
                        int openmp_holds)
{
#if FASP_USE_OPENMP
	int i;
	double *xpt=x->val;
	
	if (n>0) x->row=n;
	else n=x->row;
	
	if (val == 0.0) {
		if (n > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
				FASP_GET_START_END(myid, nthreads, n, mybegin, myend);
				memset(&xpt[mybegin], 0x0, sizeof(double)*(myend-mybegin));
			}
		}
		else {
			memset(xpt, 0x0, sizeof(double)*n);
		}
	}
	else {
		if (n > openmp_holds) {
#pragma omp parallel for private(i) schedule(static)
			for (i=0; i<n; ++i) xpt[i]=val;
		}
		else
		{
			for (i=0; i<n; ++i) xpt[i]=val;
		}
	}
#endif
}

/**
 * \fn void fasp_dvec_cp_omp (dvector *x, dvector *y, int nthreads, int openmp_holds) 
 * \brief Copy dvector x to dvector y
 *
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_dvec_cp_omp (dvector *x, 
                       dvector *y, 
                       int nthreads, 
                       int openmp_holds) 
{
#if FASP_USE_OPENMP
	int row=x->row;
	double *x_data=x->val,*y_data=y->val;
	
	y->row=row;
	if (row > openmp_holds) {
		int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			memcpy(&y_data[mybegin], &x_data[mybegin], sizeof(double)*(myend-mybegin));
		}
	}
	else {
		memcpy(y_data,x_data,row*sizeof(double));
	}
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
