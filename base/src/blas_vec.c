/*! \file blas_vec.c
 *  \brief BLAS operations for vectors
 */

#include <math.h>
#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_dvec_axpy (const REAL a, dvector *x, dvector *y)
 *
 * \brief y = a*x + y
 *
 * \param a   REAL factor a
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 *
 */
void fasp_blas_dvec_axpy (const REAL a, 
                          dvector *x, 
                          dvector *y)
{
    unsigned INT i, m=x->row;
    REAL *xpt=x->val, *ypt=y->val;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || m <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    if ((y->row-m)!=0) {
        printf("Error: two vectors have different length!\n");
        exit(ERROR_DATA_STRUCTURE);
    }
    
    if (use_openmp) {
        INT myid, mybegin, myend;
#pragma omp parallel private(myid,mybegin,myend,i) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
            for (i=mybegin; i<myend; ++i) ypt[i] += a*xpt[i];
        }
    }
    else {
        for (i=0; i<m; ++i) ypt[i] += a*xpt[i];
    }
}

/**
 * \fn void fasp_blas_dvec_axpyz (const REAL a, dvector *x, dvector *y, dvector *z) 
 *
 * \brief z = a*x + y, z is a third vector (z is cleared)
 *
 * \param a   REAL factor a
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 * \param z   Pointer to dvector z
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */
void fasp_blas_dvec_axpyz(const REAL a, 
                          dvector *x, 
                          dvector *y, 
                          dvector *z) 
{
    const INT m=x->row;
    REAL *xpt=x->val, *ypt=y->val, *zpt=z->val;
    
    if ((y->row-m)!=0) {
        printf("Error: two vectors have different length!\n");
        exit(ERROR_DATA_STRUCTURE);
    }
    
    z->row = m;

    memcpy(ypt, zpt, m*sizeof(dvector));
    fasp_blas_array_axpy(m, a, xpt, zpt);
}

/**
 * \fn REAL fasp_blas_dvec_dotprod (dvector *x, dvector *y) 
 *
 * \brief Inner product of two vectors (x,y)
 *
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \return Inner product
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */

REAL fasp_blas_dvec_dotprod (dvector *x, 
                             dvector *y) 
{
    REAL value=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val, *ypt=y->val;    
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || length <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    if (use_openmp) {
#pragma omp parallel for reduction(+:value) private(i)
        for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
    }
    else {
        for (i=0; i<length; ++i) value+=xpt[i]*ypt[i];
    }
    return value;
}


/**
 * \fn REAL fasp_dvec_relerr (dvector *x, dvector *y) 
 *
 * \brief Relative error of two dvector x and y
 *
 * \param x   Pointer to dvector x
 * \param y   Pointer to dvector y
 *
 * \return relative error ||x-y||/||x||
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012    
 */

REAL fasp_blas_dvec_relerr (dvector *x, 
                            dvector *y)
{
    REAL diff=0, temp=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val, *ypt=y->val;
    INT nthreads, use_openmp;

    if (length!=y->row) {
        printf("Error: The lengths of vectors do not match! \n");
        exit(ERROR_DUMMY_VAR);    
    }

	if(!FASP_USE_OPENMP || length <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
        
    if (use_openmp) {
#pragma omp parallel for reduction(+:temp,diff) private(i)
        for (i=0;i<length;++i) {
            temp += xpt[i]*xpt[i];
            diff += pow(xpt[i]-ypt[i],2);
        }
    }
    else {
        for (i=0;i<length;++i) {
            temp += xpt[i]*xpt[i];
            diff += pow(xpt[i]-ypt[i],2);
        }
    }
    return sqrt(diff/temp);
}

/**
 * \fn REAL fasp_blas_dvec_norm1 (dvector *x) 
 *
 * \brief L1 norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L1 norm of x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue 
 * \date   05/23/2012    
 */

REAL fasp_blas_dvec_norm1 (dvector *x) 
{
    REAL onenorm=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || length <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    if (use_openmp) {
#pragma omp parallel for reduction(+:onenorm) private(i)
        for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
    }
    else {
        for (i=0;i<length;++i) onenorm+=ABS(xpt[i]);
    }
    return onenorm;
}

/**
 * \fn REAL fasp_blas_dvec_norm2 (dvector *x) 
 *
 * \brief L2 norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L2 norm of x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue 
 * \date   05/23/2012    
 */

REAL fasp_blas_dvec_norm2 (dvector *x)
{
    REAL twonorm=0;
    INT i;
    const INT length=x->row;
    REAL *xpt=x->val;
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || length <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    if (use_openmp) {
#pragma omp parallel for reduction(+:twonorm) private(i) 
        for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
    }
    else {
        for (i=0;i<length;++i) twonorm+=xpt[i]*xpt[i];
    }
    return sqrt(twonorm);
}

/**
 * \fn REAL fasp_blas_dvec_norminf (dvector *x) 
 *
 * \brief Linf norm of dvector x
 *
 * \param x   Pointer to dvector x
 *
 * \return L_inf norm of x
 *
 * \author Chensong Zhang
 * \date   07/01/209
 */

REAL fasp_blas_dvec_norminf (dvector *x)
{
    unsigned INT i;
    const INT length=x->row;
    REAL *xpt=x->val;
    
    register REAL infnorm=0;

    for (i=0;i<length;++i) infnorm=MAX(infnorm,ABS(xpt[i]));

    return infnorm;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
