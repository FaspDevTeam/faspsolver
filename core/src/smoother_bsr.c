/*! \file smoother_bsr.c
 *  \brief Smoothers for sparse matrix in BSR format
 *
 *  Created by Zhiyang Zhou on 10/25/2010
 *  Modified by Shiquan Zhang on 11/15/2010
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

////////////////////////////////////////  JACOBI  //////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_JACOBI (dBSRmat *A, dvector *b, dvector *u)
 * \brief Jacobi relaxation
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25
 */
void fasp_smoother_dbsr_jacobi (dBSRmat *A, 
																dvector *b, 
																dvector *u)
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	double *diaginv = NULL;
	
	int nb2 = nb*nb;
	int size = ROW*nb2;
	int i,k;
	
	//! allocate memory   
	diaginv = (double *)fasp_mem_calloc(size, sizeof(double));
	
	//! get all the diagonal sub-blocks   
	for (i = 0; i < ROW; ++i)
	{
		for (k = IA[i]; k < IA[i+1]; ++k)
		{
			if (JA[k] == i)
				memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(double));
		}
	}
	
	//! compute the inverses of all the diagonal sub-blocks   
	if (nb > 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			fasp_blas_smat_inv(diaginv+i*nb2, nb);
		}
	}
	else
	{
		for (i = 0; i < ROW; ++i)
		{  
			//! zero-diagonal should be tested previously
			diaginv[i] = 1.0 / diaginv[i];
		}
	}
	
	fasp_smoother_dbsr_jacobi1(A, b, u, diaginv);
	fasp_mem_free(diaginv);	
}


/**
 * \fn void fasp_smoother_dbsr_jacobi (dBSRmat *A, dvector *b, dvector *u, double *diaginv)
 * \brief Jacobi relaxation
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param *diaginv inverses for all the diagonal blocks of A
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_jacobi1 (dBSRmat *A, 
																 dvector *b, 
																 dvector *u, 
																 double *diaginv)	
{	
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int size = ROW*nb;
	int nb2 = nb*nb;
	int i,j,k;
	int pb;
	
	//! b_tmp = b_val
	b_tmp = (double *)fasp_mem_calloc(size, sizeof(double));
	memcpy(b_tmp, b_val, size*sizeof(double));
	
	//! It's not necessary to assign the smoothing order since the result doesn't depend on it
	if (nb == 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					b_tmp[i] -= val[k]*u_val[j];
			}
		}
		
		for (i = 0; i < ROW; ++i)
		{
			u_val[i] = b_tmp[i]*diaginv[i]; 
		}      
	}
	else if (nb > 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			pb = i*nb;
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp+pb, nb);
			}
		}
		
		for (i = 0; i < ROW; ++i)
		{
			pb = i*nb;
			fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp+pb, u_val+pb, nb);
		}     
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);   
}

////////////////////////////////////  Gauss-Seidel  ////////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_GS (dBSRmat *A, dvector *b, dvector *u, int order, int *mark)
 * \brief Gauss-Seidel relaxation
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending order
 *           DESCEND      21: in descending order
 *        when mark != NULL:  in the user-defined order
 * \param *mark pointer to NULL or to the user-defined ordering
 * \return void
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_gs (dBSRmat *A, 
														dvector *b, 
														dvector *u, 
														int order, 
														int *mark )
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	double *diaginv = NULL;
	
	int nb2 = nb*nb;
	int size = ROW*nb2;
	int i,k;
	
	//! allocate memory
	diaginv = (double *)fasp_mem_calloc(size, sizeof(double));
	
	//! get all the diagonal sub-blocks   
	for (i = 0; i < ROW; ++i)
	{
		for (k = IA[i]; k < IA[i+1]; ++k)
		{
			if (JA[k] == i)
				memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(double));
		}
	}
	
	//! compute the inverses of all the diagonal sub-blocks   
	if (nb > 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			fasp_blas_smat_inv(diaginv+i*nb2, nb);
		}
	}
	else
	{
		for (i = 0; i < ROW; ++i)
		{  
			//! zero-diagonal should be tested previously
			diaginv[i] = 1.0 / diaginv[i];
		}
	}
	
	fasp_smoother_dbsr_gs1(A, b, u, order, mark, diaginv);
	fasp_mem_free(diaginv);
}


/**
 * \fn void fasp_smoother_dbsr_gs (dBSRmat *A, dvector *b, dvector *u, int order, int *mark, double *diaginv)
 * \brief Gauss-Seidel relaxation
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending order
 *           DESCEND      21: in descending order
 *        when mark != NULL:  in the user-defined order
 * \param *mark pointer to NULL or to the user-defined ordering
 * \param *diaginv inverses for all the diagonal blocks of A 
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_gs1 (dBSRmat *A, 
														 dvector *b, 
														 dvector *u, 
														 int order, 
														 int *mark, 
														 double *diaginv )	
{	
	if (!mark)
	{
		if (order == ASCEND)       //! smooth ascendingly
		{
			fasp_smoother_dbsr_gs_ascend(A, b, u, diaginv);
		}
		else if (order == DESCEND) //! smooth descendingly
		{
			fasp_smoother_dbsr_gs_descend(A, b, u, diaginv);
		}
	}
	else //! smooth according to the order 'mark' defined by user
	{
		fasp_smoother_dbsr_gs_order1(A, b, u, diaginv, mark);
	}
}

/**
 * \fn void fasp_smoother_dbsr_gs_ascend (dBSRmat *A, dvector *b, dvector *u, double *diaginv)
 * \brief Gauss-Seidel relaxation in the ascending order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv inverses for all the diagonal blocks of A 
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_gs_ascend (dBSRmat *A, 
																	 dvector *b, 
																	 dvector *u, 
																	 double *diaginv )
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int nb2 = nb*nb;
	int i,j,k;
	int pb;
	double rhs = 0.0;
	
	//! allocate memory for b_tmp
	if (nb > 1) 
	{
		b_tmp = (double *)fasp_mem_calloc(nb, sizeof(double));
	}  
	
	if (nb == 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = rhs*diaginv[i];
		}  
	}
	else if (nb > 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);
}

/**
 * \fn void fasp_smoother_dbsr_gs_descend (dBSRmat *A, dvector *b, dvector *u, double *diaginv)
 * \brief Gauss-Seidel relaxation in the descending order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv inverses for all the diagonal blocks of A 
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_gs_descend (dBSRmat *A, 
																		dvector *b, 
																		dvector *u, 
																		double *diaginv )
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int nb2 = nb*nb;
	int i,j,k;
	int pb;
	double rhs = 0.0;
	
	//! allocate memory for b_tmp
	if (nb > 1) 
	{
		b_tmp = (double *)fasp_mem_calloc(nb, sizeof(double));
	}  
	
	if (nb == 1)
	{
		for (i = ROW-1; i >= 0; i--)
		{
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = rhs*diaginv[i];
		}  
	}
	else if (nb > 1)
	{
		for (i = ROW-1; i >= 0; i--)
		{
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);
}

/**
 * \fn void fasp_smoother_dbsr_gs_order1 (dBSRmat *A, dvector *b, dvector *u, double *diaginv, int *mark)
 * \brief Gauss-Seidel relaxation in the user-defined order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv inverses for all the diagonal blocks of A
 * \param *mark pointer to the user-defined ordering  
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_gs_order1 (dBSRmat *A, 
																	 dvector *b, 
																	 dvector *u, 
																	 double *diaginv, 
																	 int *mark )
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int nb2 = nb*nb;
	int i,j,k;
	int I,pb;
	double rhs = 0.0;
	
	//! allocate memory for b_tmp
	if (nb > 1) 
	{
		b_tmp = (double *)fasp_mem_calloc(nb, sizeof(double));
	}  
	
	if (nb == 1)
	{
		for (I = 0; I < ROW; ++I)
		{
			i = mark[I];
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = rhs*diaginv[i];
		}  
	}
	else if (nb > 1)
	{
		for (I = 0; I < ROW; ++I)
		{
			i = mark[I];
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);
}

/**
 * \fn void fasp_smoother_dbsr_gs_order2 (dBSRmat *A, dvector *b, dvector *u, int *mark)
 * \brief Gauss-Seidel relaxation in the user-defined order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *mark pointer to the user-defined ordering  
 * \return void
 *
 * \note The only difference between the functions 'fasp_smoother_dbsr_gs_order2' 
 *       and 'fasp_smoother_dbsr_gs_order1' lies in that we don't have to multiply 
 *       by the inverses of the diagonal blocks in each ROW since matrix A has  
 *       been such scaled that all the diagonal blocks become identity matrices.
 *
 * \author Zhiyang Zhou
 * \date 2010/11/08 
 */
void fasp_smoother_dbsr_gs_order2 (dBSRmat *A, 
																	 dvector *b, 
																	 dvector *u, 
																	 int *mark, 
																	 double *work)
{
	//! members of A 
	const int     ROW = A->ROW;
	const int     nb  = A->nb;
	const int     nb2 = nb*nb;
	int          *IA  = A->IA;
	int          *JA  = A->JA;   
	double       *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = work;
	
	//! local variables
	int i,j,k,I,pb;
	double rhs = 0.0;
	
	if (nb == 1)
	{
		for (I = 0; I < ROW; ++I)
		{
			i = mark[I];
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = rhs;
		}  
	}
	else if (nb > 1)
	{
		for (I = 0; I < ROW; ++I)
		{
			i = mark[I];
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			memcpy(u_val+pb, b_tmp, nb*sizeof(double));
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
}

////////////////////////////////////////  SOR  //////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_sor (dBSRmat *A, dvector *b, dvector *u, int order, int *mark, double weight)
 * \brief SOR relaxation
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending order
 *           DESCEND      21: in descending order
 *        when mark != NULL:  in the user-defined order
 * \param *mark pointer to NULL or to the user-defined ordering
 * \param weight over-relaxation parameter 
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_sor (dBSRmat *A, 
														 dvector *b, 
														 dvector *u, 
														 int order, 
														 int *mark, 
														 double weight)
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	double *diaginv = NULL;
	
	int nb2 = nb*nb;
	int size = ROW*nb2;
	int i,k;
	
	//! allocate memory
	diaginv = (double *)fasp_mem_calloc(size, sizeof(double));
	
	//! get all the diagonal sub-blocks   
	for (i = 0; i < ROW; ++i)
	{
		for (k = IA[i]; k < IA[i+1]; ++k)
		{
			if (JA[k] == i)
				memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(double));
		}
	}
	
	//! compute the inverses of all the diagonal sub-blocks   
	if (nb > 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			fasp_blas_smat_inv(diaginv+i*nb2, nb);
		}
	}
	else
	{
		for (i = 0; i < ROW; ++i)
		{  
			//! zero-diagonal should be tested previously
			diaginv[i] = 1.0 / diaginv[i];
		}
	}
	
	fasp_smoother_dbsr_sor1(A, b, u, order, mark, diaginv, weight);
	fasp_mem_free(diaginv);
}

/**
 * \fn void fasp_smoother_dbsr_sor1 (dBSRmat *A, dvector *b, dvector *u, 
 *                         int order, int *mark, double *diaginv, double weight)
 * \brief SOR relaxation
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending order
 *           DESCEND      21: in descending order
 *        when mark != NULL:  in the user-defined order
 * \param *mark pointer to NULL or to the user-defined ordering
 * \param *diaginv inverses for all the diagonal blocks of A
 * \param weight over-relaxation parameter  
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_sor1 (dBSRmat *A, 
															dvector *b, 
															dvector *u, 
															int order, 
															int *mark, 
															double *diaginv, 
															double weight)	
{	
	if (!mark)
	{
		if (order == ASCEND)       //! smooth ascendingly
		{
			fasp_smoother_dbsr_sor_ascend(A, b, u, diaginv, weight);
		}
		else if (order == DESCEND) //! smooth descendingly
		{
			fasp_smoother_dbsr_sor_descend(A, b, u, diaginv, weight);
		}
	}
	else //! smooth according to the order 'mark' defined by user
	{
		fasp_smoother_dbsr_sor_order(A, b, u, diaginv, mark, weight);
	}
}

/**
 * \fn void fasp_smoother_dbsr_sor_ascend (dBSRmat *A, dvector *b, dvector *u, double *diaginv, double weight)
 * \brief SOR relaxation in the ascending order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv inverses for all the diagonal blocks of A
 * \param weight over-relaxation parameter   
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_sor_ascend (dBSRmat *A, 
																		dvector *b, 
																		dvector *u, 
																		double *diaginv, 
																		double weight )
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int nb2 = nb*nb;
	int i,j,k;
	int pb;
	double rhs = 0.0;
	double one_minus_weight = 1.0 - weight;
	
	//! allocate memory for b_tmp
	if (nb > 1) 
	{
		b_tmp = (double *)fasp_mem_calloc(nb, sizeof(double));
	}  
	
	if (nb == 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = one_minus_weight*u_val[i] + weight*(rhs*diaginv[i]);         
		}  
	}
	else if (nb > 1)
	{
		for (i = 0; i < ROW; ++i)
		{
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			fasp_blas_smat_aAxpby(weight, diaginv+nb2*i, b_tmp, one_minus_weight, u_val+pb, nb);         
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);
}

/**
 * \fn void fasp_smoother_dbsr_sor_descend (dBSRmat *A, dvector *b, dvector *u, double *diaginv, double weight)
 * \brief SOR relaxation in the descending order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv inverses for all the diagonal blocks of A
 * \param weight over-relaxation parameter   
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_sor_descend (dBSRmat *A, 
																		 dvector *b, 
																		 dvector *u, 
																		 double *diaginv, 
																		 double weight)
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int nb2 = nb*nb;
	int i,j,k;
	int pb;
	double rhs = 0.0;
	double one_minus_weight = 1.0 - weight;
	
	//! allocate memory for b_tmp
	if (nb > 1) 
	{
		b_tmp = (double *)fasp_mem_calloc(nb, sizeof(double));
	}  
	
	if (nb == 1)
	{
		for (i = ROW-1; i >= 0; i--)
		{
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = one_minus_weight*u_val[i] + weight*(rhs*diaginv[i]);         
		}  
	}
	else if (nb > 1)
	{
		for (i = ROW-1; i >= 0; i--)
		{
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			fasp_blas_smat_aAxpby(weight, diaginv+nb2*i, b_tmp, one_minus_weight, u_val+pb, nb);         
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);
}

/**
 * \fn void fasp_smoother_dbsr_sor_order (dBSRmat *A, dvector *b, dvector *u, double *diaginv, double weight)
 * \brief SOR relaxation in the user-defined order
 * \param *A pointer to coefficient matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv inverses for all the diagonal blocks of A
 * \param *mark pointer to the user-defined ordering   
 * \param weight over-relaxation parameter   
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_sor_order (dBSRmat *A, 
																	 dvector *b, 
																	 dvector *u, 
																	 double *diaginv, 
																	 int *mark, 
																	 double weight )
{
	//! members of A 
	int     ROW = A->ROW;
	int     nb  = A->nb;
	int    *IA  = A->IA;
	int    *JA  = A->JA;   
	double *val = A->val;
	
	//! values of dvector b and u
	double *b_val = b->val;
	double *u_val = u->val;
	
	//! auxiliary array
	double *b_tmp = NULL;
	
	//! local variables
	int nb2 = nb*nb;
	int i,j,k;
	int I,pb;
	double rhs = 0.0;
	double one_minus_weight = 1.0 - weight;
	
	//! allocate memory for b_tmp
	if (nb > 1) 
	{
		b_tmp = (double *)fasp_mem_calloc(nb, sizeof(double));
	} 
	
	if (nb == 1)
	{
		for (I = 0; I < ROW; ++I)
		{
			i = mark[I];
			rhs = b_val[i];
			for (k = IA[i]; k < IA[i+1]; ++k)
			{  
				j = JA[k];
				if (j != i)
					rhs -= val[k]*u_val[j];
			}
			u_val[i] = one_minus_weight*u_val[i] + weight*(rhs*diaginv[i]);         
		}  
	}
	else if (nb > 1)
	{
		for (I = 0; I < ROW; ++I)
		{
			i = mark[I];
			pb = i*nb;
			memcpy(b_tmp, b_val+pb, nb*sizeof(double));
			for (k = IA[i]; k < IA[i+1]; ++k)
			{ 
				j = JA[k];
				if (j != i)
					fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
			}
			fasp_blas_smat_aAxpby(weight, diaginv+nb2*i, b_tmp, one_minus_weight, u_val+pb, nb);         
		}
	}
	else
	{
		printf("\n nb is illegal!\n\n");
		return;
	}
	
	fasp_mem_free(b_tmp);
}

////////////////////////////////////////  ILU  //////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_ilu (dBSRmat *A, dvector *b, dvector *x, void *data)
 * \brief ILU method as the smoother in solving Au=b with multigrid method
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param *x pointer to current solution
 * \param *data pointer to user defined data
 * \return void
 *
 * \author Zhiyang Zhou
 * \date 2010/10/25 
 */
void fasp_smoother_dbsr_ilu (dBSRmat *A, 
														 dvector *b, 
														 dvector *x, 
														 void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const unsigned int nb=iludata->nb,m=A->ROW*nb, memneed=6*m;
	double *xval = x->val, *bval = b->val;
	double *zz, *zr, *z;
	
	if (iludata->nwork<memneed) goto MEMERR; 
	
	zz = iludata->work + 3*m;
	zr = zz + m;
	z = zr + m;
	
	/** form residual zr = b - A x */
	fasp_array_cp(m,bval,zr); fasp_blas_dbsr_aAxpy(-1.0,A,xval,zr);
	
	/** solve LU z=zr */
	fasp_precond_dbsr_ilu(zr,z,iludata);
	
	/** x=x+z */
	fasp_blas_array_axpy(m,1,z,xval);		
	
	
	return;
	
MEMERR:
	printf("Error: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
