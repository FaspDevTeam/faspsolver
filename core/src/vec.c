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
 * \param m    Number of rows
 *
 * \return u   The new dvector
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
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
 * \param m   Number of rows
 *
 * \return u  The new ivector
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
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
 * \param m    Number of rows
 * \param u    Pointer to the REAL vector (OUTPUT)
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
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
 * \param m   Number of rows
 * \param u   Pointer to the integer vector (OUTPUT)
 *
 * \author Chensong Zhang 
 * \date   2010/04/06
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
 * \param u   Pointer to the vector which needs to be deallocated
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
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
 * \param u   Pointer to the vector which needs to be deallocated
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
 *
 * \note This function is same as fasp_dvec_free except input type.
 */
void fasp_ivec_free (ivector *u)
{		
	if (u==NULL) return;
    
	fasp_mem_free(u->val);
	u->row = 0; u->val = NULL; 
}

/**
 * \fn void fasp_dvec_null (dvector *x) 
 *
 * \brief Initialize dvector
 *
 * \param x   Pointer to dvector which needs to be initialized
 *
 * \author Chensong Zhang
 * \date 2010/04/03  
 */
void fasp_dvec_null (dvector *x) 
{
	x->row = 0; x->val = NULL;
}

/**
 * \fn void fasp_dvec_rand (INT n, dvector *x)
 *
 * \brief Generate random REAL vector in the range from 0 to 1
 *
 * \param n    Size of the vector
 * \param x    Pointer to a dvector
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
 * \date   11/16/2009
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
 * \param n      Number of variables
 * \param x      Pointer to dvector
 * \param val    Initial value for the dvector
 *
 * \author Chensong Zhang
 * \date   11/16/2009
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
 * \param  m    Integer value of ivector
 * \param  u    Pointer to the vector (MODIFIED)
 *
 * \author Chensong Zhang
 * \date   2010/04/03  
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
 * \param x  Pointer to dvector
 * \param y  Pointer to dvector (MODIFIED)
 *
 * \author Chensong Zhang
 * \date   11/16/2009
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
 * \param  x    Pointer to dvector
 * \param  y    Pointer to dvector
 *
 * \return      Maximal norm of x-y
 *
 * \author Chensong Zhang
 * \date   11/16/2009
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
 * \param b       Pointer to the dvector
 * \param diag    Pointer to the diagonal entries
 *
 * \author Xiaozhe Hu
 * \date   01/31/2011
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
		printf("### ERROR: Size of diag = %d and size of dvector = %d mismatch!!", 
               diag->row, n);
		exit(ERROR_MISC);
	}
	
	// main loop
	for (i=0; i<n; i++) val[i] = val[i]/sqrt(diag->val[i]);
	
	return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
