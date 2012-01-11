/*! \file blas_csr_omp.c
 *  \brief BLAS operations for sparse matrices in CSR format.
 *
 *  \note Sparse functions usually contain three runs. The three runs are all the 
 *  same but thy serve different purpose. 
 * 
 *  Example: If you do c=a+b: 
 *			- first do a dry run to find the number of non-zeroes in the result and form ic; 
 *			- allocate space (memory) for jc and form this one; If you only care about a "boolean" result of the addition, you stop here. 
 *			- you call another routine, which uses ic and jc to perform the addition. 
 *
 *		We need to redo these routines later in this way!
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*-----------------------------------omp--------------------------------------*/
/* Feng Chunsheng Yue Xiaoqiang  Expand the inner loop  /mar/14/2011/  */
/**
 * \fn void fasp_blas_dcsr_mxv_omp (dCSRmat *A, double *x, double *y, int nthreads, 
 *                                  int openmp_holds)
 * \brief Matrix-vector multiplication y = A*x
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/14/2011
 * \date Jan/11/2012 modified by  FENG Chunsheng
 */
void fasp_blas_dcsr_mxv_omp (dCSRmat *A, 
														 double *x, 
														 double *y, 
														 int nthreads, 
														 int openmp_holds)
{
#if FASP_USE_OPENMP
	const int m=A->row;
	const int *ia=A->IA, *ja=A->JA;
	const double *aj=A->val;
	int i, k, begin_row, end_row, nnz_num_row;
	register double temp;
	
	if (m > openmp_holds) {
		int myid, mybegin, myend;
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, nnz_num_row, k) 
		for (myid = 0; myid < nthreads; myid++ )
		{
			FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
			for (i=mybegin; i<myend; ++i)
			{
				temp=0.0; 
				begin_row = ia[i];
				end_row = ia[i+1];
				nnz_num_row = end_row - begin_row;
				switch(nnz_num_row)
				{
				case 3:
					k = begin_row;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					break;
				case 4:
					k = begin_row;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					break;
				case 5:
					k = begin_row;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					break;
				case 6:
					k = begin_row;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					break;
				case 7:
					k = begin_row;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					k ++;
					temp += aj[k]*x[ja[k]];
					break;
				default:
					for (k=begin_row; k<end_row; ++k)
					{
						temp += aj[k]*x[ja[k]];
					}
					break;
				}
				y[i]=temp;
			}
		}
	}
	else {
		for (i=0;i<m;++i) {
			temp=0.0;
			begin_row=ia[i];
			end_row=ia[i+1];
			nnz_num_row = end_row - begin_row;
			switch(nnz_num_row)
			{
			case 3:
				k=begin_row;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				break;
			case 4:
				k=begin_row;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				break;
			case 5:
				k=begin_row;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				break;
			case 6:
				k=begin_row;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				break;
			case 7:
				k=begin_row;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				k ++;
				temp+=aj[k]*x[ja[k]];
				break;
			default:
				for (k=begin_row; k<end_row; ++k)
				{
					temp+=aj[k]*x[ja[k]];
				}
				break;
			}
			y[i]=temp;
		}
	}
#endif
}

/**
 * \fn void fasp_blas_dcsr_aAxpy_omp (const double alpha, dCSRmat *A, double *x,
 *                                    double *y, int nthreads, int openmp_holds)
 * \brief Matrix-vector multiplication y = alpha*A*x + y
 * \param alpha real number
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_blas_dcsr_aAxpy_omp (const double alpha, 
															 dCSRmat *A, 
															 double *x, 
															 double *y, 
															 int nthreads, 
															 int openmp_holds)
{
#if FASP_USE_OPENMP
	const int  m  = A->row;
	const int *ia = A->IA, *ja = A->JA;
	const double *aj = A->val;
	int i, k, begin_row, end_row;
	register double temp;
	
	if ( alpha == 1.0 ) {
		if (m > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k) 
			for (myid = 0; myid < nthreads; myid++ )
			{
				FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
				for (i=mybegin; i<myend; ++i)
				{
					temp=0.0; 
					begin_row=ia[i]; end_row=ia[i+1]; 
					for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
					y[i]+=temp;
				}
			}
		}
		else {
			for (i=0;i<m;++i) {
				temp=0.0; 
				begin_row=ia[i]; end_row=ia[i+1]; 
				for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
				y[i]+=temp;	
			}
		}
	}
	
	else if ( alpha == -1.0 ) {
		if (m > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k) 
			for (myid = 0; myid < nthreads; myid++ )
			{
				FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
				for (i=mybegin; i<myend; ++i)
				{
					temp=0.0; 
					begin_row=ia[i]; end_row=ia[i+1]; 
					for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
					y[i]-=temp;
				}
			}
		}
		else {
			for (i=0;i<m;++i) {
				temp=0.0; 
				begin_row=ia[i]; end_row=ia[i+1]; 
				for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
				y[i]-=temp;
			}
		}
	}
	
	else {
		if (m > openmp_holds) {
			int myid, mybegin, myend;
#pragma omp parallel for private(myid, mybegin, myend, i, temp, begin_row, end_row, k) 
			for (myid = 0; myid < nthreads; myid++ )
			{
				FASP_GET_START_END(myid, nthreads, m, mybegin, myend);
				for (i=mybegin; i<myend; ++i)
				{
					temp=0.0;
					begin_row=ia[i]; end_row=ia[i+1];
					for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
					y[i]+=temp*alpha;
				}
			}
		}
		else {
			for (i=0;i<m;++i) {
				temp=0.0; 
				begin_row=ia[i]; end_row=ia[i+1]; 
				for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
				y[i]+=temp*alpha;
			}
		}
	}
#endif
}

/**
 * \fn double fasp_blas_dcsr_vmv_omp (dCSRmat *A, double *x, double *y, 
 *                                    int nthreads, int openmp_holds)
 * \brief vector-Matrix-vector multiplication alpha = y'*A*x
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
double fasp_blas_dcsr_vmv_omp (dCSRmat *A, 
															 double *x, 
															 double *y, 
															 int nthreads, 
															 int openmp_holds)
{
	register double value=0.0;
#if FASP_USE_OPENMP
	const int m=A->row;
	const int *ia=A->IA, *ja=A->JA;
	const double *aj=A->val;
	int i, k, begin_row, end_row;
	register double temp;
	
	if (m > openmp_holds) {
#pragma omp parallel for reduction(+:value) private(i,temp,begin_row,end_row,k) 
		for (i=0;i<m;++i) {
			temp=0.0;
			begin_row=ia[i]; end_row=ia[i+1];
			for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
			value+=y[i]*temp;
		}
	}
	else {
		for (i=0;i<m;++i) {
			temp=0.0;
			begin_row=ia[i]; end_row=ia[i+1];
			for (k=begin_row; k<end_row; ++k) temp+=aj[k]*x[ja[k]];
			value+=y[i]*temp;
		}
	}
#endif	
	return value;
}

/**
 * \fn void fasp_blas_dcsr_rap_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP, int nthreads, int openmp_holds)
 * \brief Triple sparse matrix multiplication RAP=R*A*P
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *RAP pointer to dCSRmat matrix equal to R*A*P
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return void
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/11/2010
 */
void fasp_blas_dcsr_rap_omp( dCSRmat  *R,
                             dCSRmat  *A,
                             dCSRmat  *P,
                             dCSRmat  *RAP,
                             int       nthreads,
                             int       openmp_holds )
{
#if FASP_USE_OPENMP
    int n_coarse = R->row;
    int *R_i = R->IA;
    int *R_j = R->JA;
    double *R_data = R->val;
    
    int n_fine = A->row;
    int *A_i = A->IA;
    int *A_j = A->JA;
    double *A_data = A->val;
    
    int *P_i = P->IA;
    int *P_j = P->JA;
    double *P_data = P->val;
    
    int RAP_size;
    int *RAP_i = NULL;
    int *RAP_j = NULL;
    double *RAP_data = NULL;
    
    int *P_marker = NULL;
    int *A_marker = NULL;
    
    int *Ps_marker = NULL;
    int *As_marker = NULL;
    int *RAP_temp = NULL;
    int *part_end = NULL;
    
    int ic, i1, i2, i3, jj1, jj2, jj3;
    int jj_counter, jj_row_begining;
    double r_entry, r_a_product, r_a_p_product;
    
    // Save the memory
    if (n_coarse <= openmp_holds)
    {
        nthreads = 1;
    }
    
    int coarse_mul_nthreads = n_coarse * nthreads;
    int fine_mul_nthreads = n_fine * nthreads;
    int coarse_add_nthreads = n_coarse + nthreads;
    int minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    int total_calloc = minus_one_length + coarse_add_nthreads + nthreads;
    
    Ps_marker = (int *)fasp_mem_calloc(total_calloc, sizeof(int));
    As_marker = Ps_marker + coarse_mul_nthreads;
    RAP_temp = As_marker + fine_mul_nthreads;
    part_end = RAP_temp + coarse_add_nthreads;
    
   /*------------------------------------------------------*
    *  First Pass: Determine size of RAP and set up RAP_i  *
    *------------------------------------------------------*/
    RAP_i = (int *)fasp_mem_calloc(n_coarse+1, sizeof(int));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        int myid, mybegin, myend, Ctemp;
#pragma omp parallel private(myid, mybegin, myend, P_marker, A_marker, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3)
        {
            myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n_coarse, mybegin, myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            jj_counter = 0;
            for (ic = mybegin; ic < myend; ic ++)
            {
                P_marker[ic] = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter ++;
                
                for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
                {
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                    {
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic)
                        {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                            {
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining)
                                {
                                    P_marker[i3] = jj_counter;
                                    jj_counter ++;
                                }
                            }
                        }
                    }
                }
                
                RAP_temp[ic+myid] = jj_row_begining;
            }
            RAP_temp[myend+myid] = jj_counter;
            
            part_end[myid] = myend + myid + 1;
        }
        fasp_iarray_cp_omp(part_end[0], RAP_temp, RAP_i, nthreads, openmp_holds);
        jj_counter = part_end[0];
        Ctemp = 0;
        for (i1 = 1; i1 < nthreads; i1 ++)
        {
            Ctemp += RAP_temp[part_end[i1-1]-1];
            for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++)
            {
                RAP_i[jj_counter] = RAP_temp[jj1] + Ctemp;
                jj_counter ++;
            }
        }
        RAP_size = RAP_i[n_coarse];
    }
    else
    {
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic ++)
        {
            Ps_marker[ic] = jj_counter;
            jj_row_begining = jj_counter;
            jj_counter ++;
            
            for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
            {
                i1 = R_j[jj1];
                
                for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                {
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic)
                    {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                        {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining)
                            {
                                Ps_marker[i3] = jj_counter;
                                jj_counter ++;
                            }
                        }
                    }
                }
            }
            
            RAP_i[ic] = jj_row_begining;
        }
        
        RAP_i[n_coarse] = jj_counter;
        RAP_size = jj_counter;
    }
    printf("A_H NNZ = %d\n", RAP_size);
    RAP_j = (int *)fasp_mem_calloc(RAP_size, sizeof(int));
    RAP_data = (double *)fasp_mem_calloc(RAP_size, sizeof(double));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, P_marker, A_marker, jj_counter, ic, jj_row_begining,\
                                      jj1, r_entry, i1, jj2, r_a_product, i2, jj3, r_a_p_product, i3)
        {
            myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n_coarse, mybegin, myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            jj_counter = RAP_i[mybegin];
            for (ic = mybegin; ic < myend; ic ++)
            {
                P_marker[ic] = jj_counter;
                jj_row_begining = jj_counter;
                RAP_j[jj_counter] = ic;
                RAP_data[jj_counter] = 0.0;
                jj_counter ++;
                for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
                {
                    r_entry = R_data[jj1];
                    
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                    {
                        r_a_product = r_entry * A_data[jj2];
                        
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic)
                        {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                            {
                                r_a_p_product = r_a_product * P_data[jj3];
                                
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining)
                                {
                                    P_marker[i3] = jj_counter;
                                    RAP_data[jj_counter] = r_a_p_product;
                                    RAP_j[jj_counter] = i3;
                                    jj_counter ++;
                                }
                                else
                                {
                                    RAP_data[P_marker[i3]] += r_a_p_product;
                                }
                            }
                        }
                        else
                        {
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                            {
                                i3 = P_j[jj3];
                                r_a_p_product = r_a_product * P_data[jj3];
                                RAP_data[P_marker[i3]] += r_a_p_product;
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic ++)
        {
            Ps_marker[ic] = jj_counter;
            jj_row_begining = jj_counter;
            RAP_j[jj_counter] = ic;
            RAP_data[jj_counter] = 0.0;
            jj_counter ++;
            
            for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
            {
                r_entry = R_data[jj1];
                
                i1 = R_j[jj1];
                for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                {
                    r_a_product = r_entry * A_data[jj2];
                    
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic)
                    {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                        {
                            r_a_p_product = r_a_product * P_data[jj3];
                            
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining)
                            {
                                Ps_marker[i3] = jj_counter;
                                RAP_data[jj_counter] = r_a_p_product;
                                RAP_j[jj_counter] = i3;
                                jj_counter ++;
                            }
                            else
                            {
                                RAP_data[Ps_marker[i3]] += r_a_p_product;
                            }
                        }
                    }
                    else
                    {
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                        {
                            i3 = P_j[jj3];
                            r_a_p_product = r_a_product * P_data[jj3];
                            RAP_data[Ps_marker[i3]] += r_a_p_product;
                        }
                    }
                }
            }
        }
    }
    
    RAP->row = n_coarse;
    RAP->col = n_coarse;
    RAP->nnz = RAP_size;
    RAP->IA = RAP_i;
    RAP->JA = RAP_j;
    RAP->val = RAP_data;
    
    fasp_mem_free(Ps_marker);
#endif
}

/**
 * \fn void fasp_blas_dcsr_rap_agg_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP, int nthreads, int openmp_holds)
 * \brief Triple sparse matrix multiplication RAP=R*A*P, where the entries of R and P are all ones.
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *RAP pointer to dCSRmat matrix equal to R*A*P
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return void
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/14/2010
 */
void fasp_blas_dcsr_rap_agg_omp( dCSRmat  *R,
                                 dCSRmat  *A,
                                 dCSRmat  *P,
                                 dCSRmat  *RAP,
                                 int       nthreads,
                                 int       openmp_holds )
{
#if FASP_USE_OPENMP
    int n_coarse = R->row;
    int *R_i = R->IA;
    int *R_j = R->JA;
    
    int n_fine = A->row;
    int *A_i = A->IA;
    int *A_j = A->JA;
    double *A_data = A->val;
    
    int *P_i = P->IA;
    int *P_j = P->JA;
    
    int RAP_size;
    int *RAP_i = NULL;
    int *RAP_j = NULL;
    double *RAP_data = NULL;
    
    int *P_marker = NULL;
    int *A_marker = NULL;
    
    int *Ps_marker = NULL;
    int *As_marker = NULL;
    int *RAP_temp = NULL;
    int *part_end = NULL;
    
    int ic, i1, i2, i3, jj1, jj2, jj3;
    int jj_counter, jj_row_begining;
    double a_entry;
    
    // Save the memory
    if (n_coarse <= openmp_holds)
    {
        nthreads = 1;
    }
    
    int coarse_mul_nthreads = n_coarse * nthreads;
    int fine_mul_nthreads = n_fine * nthreads;
    int coarse_add_nthreads = n_coarse + nthreads;
    int minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    int total_calloc = minus_one_length + coarse_add_nthreads + nthreads;
    
    Ps_marker = (int *)fasp_mem_calloc(total_calloc, sizeof(int));
    As_marker = Ps_marker + coarse_mul_nthreads;
    RAP_temp = As_marker + fine_mul_nthreads;
    part_end = RAP_temp + coarse_add_nthreads;
    
   /*------------------------------------------------------*
    *  First Pass: Determine size of RAP and set up RAP_i  *
    *------------------------------------------------------*/
    RAP_i = (int *)fasp_mem_calloc(n_coarse+1, sizeof(int));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        int myid, mybegin, myend, Ctemp;
#pragma omp parallel private(myid, mybegin, myend, P_marker, A_marker, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3)
        {
            myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n_coarse, mybegin, myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            jj_counter = 0;
            for (ic = mybegin; ic < myend; ic ++)
            {
                P_marker[ic] = jj_counter;
                jj_row_begining = jj_counter;
                jj_counter ++;
                
                for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
                {
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                    {
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic)
                        {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                            {
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining)
                                {
                                    P_marker[i3] = jj_counter;
                                    jj_counter ++;
                                }
                            }
                        }
                    }
                }
                
                RAP_temp[ic+myid] = jj_row_begining;
            }
            RAP_temp[myend+myid] = jj_counter;
            
            part_end[myid] = myend + myid + 1;
        }
        fasp_iarray_cp_omp(part_end[0], RAP_temp, RAP_i, nthreads, openmp_holds);
        jj_counter = part_end[0];
        Ctemp = 0;
        for (i1 = 1; i1 < nthreads; i1 ++)
        {
            Ctemp += RAP_temp[part_end[i1-1]-1];
            for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++)
            {
                RAP_i[jj_counter] = RAP_temp[jj1] + Ctemp;
                jj_counter ++;
            }
        }
        RAP_size = RAP_i[n_coarse];
    }
    else
    {
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic ++)
        {
            Ps_marker[ic] = jj_counter;
            jj_row_begining = jj_counter;
            jj_counter ++;
            
            for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
            {
                i1 = R_j[jj1];
                
                for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                {
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic)
                    {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                        {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining)
                            {
                                Ps_marker[i3] = jj_counter;
                                jj_counter ++;
                            }
                        }
                    }
                }
            }
            
            RAP_i[ic] = jj_row_begining;
        }
        
        RAP_i[n_coarse] = jj_counter;
        RAP_size = jj_counter;
    }
    
    RAP_j = (int *)fasp_mem_calloc(RAP_size, sizeof(int));
    RAP_data = (double *)fasp_mem_calloc(RAP_size, sizeof(double));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        int myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, P_marker, A_marker, jj_counter, ic, jj_row_begining,\
                                      jj1, i1, jj2, a_entry, i2, jj3, i3)
        {
            myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, n_coarse, mybegin, myend);
            P_marker = Ps_marker + myid * n_coarse;
            A_marker = As_marker + myid * n_fine;
            jj_counter = RAP_i[mybegin];
            for (ic = mybegin; ic < myend; ic ++)
            {
                P_marker[ic] = jj_counter;
                jj_row_begining = jj_counter;
                RAP_j[jj_counter] = ic;
                RAP_data[jj_counter] = 0.0;
                jj_counter ++;
                for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
                {
                    i1 = R_j[jj1];
                    for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                    {
                        a_entry = A_data[jj2];
                        
                        i2 = A_j[jj2];
                        if (A_marker[i2] != ic)
                        {
                            A_marker[i2] = ic;
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                            {
                                i3 = P_j[jj3];
                                if (P_marker[i3] < jj_row_begining)
                                {
                                    P_marker[i3] = jj_counter;
                                    RAP_data[jj_counter] = a_entry;
                                    RAP_j[jj_counter] = i3;
                                    jj_counter ++;
                                }
                                else
                                {
                                    RAP_data[P_marker[i3]] += a_entry;
                                }
                            }
                        }
                        else
                        {
                            for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                            {
                                i3 = P_j[jj3];
                                RAP_data[P_marker[i3]] += a_entry;
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        jj_counter = 0;
        for (ic = 0; ic < n_coarse; ic ++)
        {
            Ps_marker[ic] = jj_counter;
            jj_row_begining = jj_counter;
            RAP_j[jj_counter] = ic;
            RAP_data[jj_counter] = 0.0;
            jj_counter ++;
            
            for (jj1 = R_i[ic]; jj1 < R_i[ic+1]; jj1 ++)
            {
                i1 = R_j[jj1];
                for (jj2 = A_i[i1]; jj2 < A_i[i1+1]; jj2 ++)
                {
                    a_entry = A_data[jj2];
                    
                    i2 = A_j[jj2];
                    if (As_marker[i2] != ic)
                    {
                        As_marker[i2] = ic;
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                        {
                            i3 = P_j[jj3];
                            if (Ps_marker[i3] < jj_row_begining)
                            {
                                Ps_marker[i3] = jj_counter;
                                RAP_data[jj_counter] = a_entry;
                                RAP_j[jj_counter] = i3;
                                jj_counter ++;
                            }
                            else
                            {
                                RAP_data[Ps_marker[i3]] += a_entry;
                            }
                        }
                    }
                    else
                    {
                        for (jj3 = P_i[i2]; jj3 < P_i[i2+1]; jj3 ++)
                        {
                            i3 = P_j[jj3];
                            RAP_data[Ps_marker[i3]] += a_entry;
                        }
                    }
                }
            }
        }
    }
    
    RAP->row = n_coarse;
    RAP->col = n_coarse;
    RAP->nnz = RAP_size;
    RAP->IA = RAP_i;
    RAP->JA = RAP_j;
    RAP->val = RAP_data;
    
    fasp_mem_free(Ps_marker);
#endif
}

/**
 * \fn void fasp_blas_dcsr_rap2_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, int nthreads, int openmp_holds)
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *B   pointer to dCSRmat matrix equal to R*A*P
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return     void
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 04/25/2011
 */
void fasp_blas_dcsr_rap1_omp (dCSRmat *R, 
												 dCSRmat *A, 
												 dCSRmat *P, 
												 dCSRmat *B,
												 int nthreads,
												 int openmp_holds)
{
#if FASP_USE_OPENMP
	if (R->row > openmp_holds) {
		const int row=R->row, col=P->col;
		int *ir=R->IA, *ia=A->IA, *ip=P->IA;
		int *jr=R->JA, *ja=A->JA, *jp=P->JA;
		double *rj=R->val, *aj=A->val, *pj=P->val;
		int istart, iistart;
		int end_row, end_rowA, end_rowR;
		int i, j, jj, k, length, myid, mybegin, myend, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
		int *index = NULL;
		int *iindex = NULL;
		int *BTindex = NULL;
		double *temp = NULL;
		
		//printf(" >>> nnz = %d, row = %d, sparsity: %le\n", A->nnz, A->row, ((double)A->nnz)/A->row/A->row);
		
		int *indexs=(int *)fasp_mem_calloc(A->col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(int);
#endif
		int *iindexs=(int *)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
		int *iac=(int*)fasp_mem_calloc(row+1,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (row+1)*sizeof(int);
#endif
		int *part_end=(int*)fasp_mem_calloc(2*nthreads+row,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (2*nthreads+row)*sizeof(int);
#endif
		int *iac_temp=part_end+nthreads;
		fasp_iarray_set_omp(A->col*nthreads, indexs, -2, nthreads, openmp_holds);
		fasp_iarray_set_omp(col*nthreads, iindexs, -2, nthreads, openmp_holds);
		// First loop: form iac of R*A*P
#pragma omp parallel private(myid, mybegin, myend, index, iindex, jj_counter, ic, jj_row_begining, end_rowR, jj1, i1, end_rowA, jj2, i2, end_row, jj3, i3)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			index = indexs + myid*A->col;
			iindex = iindexs + myid*col;
			jj_counter = 0;
			for (ic = mybegin; ic < myend; ic ++)
			{
				iindex[ic] = jj_counter;
				jj_row_begining = jj_counter;
				jj_counter ++;
				end_rowR = ir[ic+1];
				for (jj1 = ir[ic]; jj1 < end_rowR; jj1 ++)
				{
					i1 = jr[jj1];
					end_rowA = ia[i1+1];
					for (jj2 = ia[i1]; jj2 < end_rowA; jj2 ++)
					{
						i2 = ja[jj2];
						if (index[i2] != ic)
						{
							index[i2] = ic;
							end_row = ip[i2+1];
							for (jj3 = ip[i2]; jj3 < end_row; jj3 ++)
							{
								i3 = jp[jj3];
								if (iindex[i3] < jj_row_begining)
								{
									iindex[i3] = jj_counter;
									jj_counter ++;
								}
							}
						}
					}
				}
				iac_temp[ic+myid] = jj_row_begining;
			}
			iac_temp[myend+myid] = jj_counter;
			part_end[myid] = myend + myid + 1;
		}
		fasp_iarray_cp_omp(part_end[0], iac_temp, iac, nthreads, openmp_holds);
		jj_counter = part_end[0];
		int Ctemp = 0;
		for (i1 = 1; i1 < nthreads; i1 ++)
		{
			Ctemp += iac_temp[part_end[i1-1]-1];
			for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++)
			{
				iac[jj_counter] = iac_temp[jj1] + Ctemp;
				jj_counter ++;
			}
		}
		//printf("A_H NNZ = %d\n", iac[row]);
		int *jac=(int*)fasp_mem_calloc(iac[row],sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(int);
#endif
		fasp_iarray_set_omp(A->col*nthreads, indexs, -2, nthreads, openmp_holds);
		fasp_iarray_set_omp(col*nthreads, iindexs, -2, nthreads, openmp_holds);
		// Second loop: form jac of R*A*P
//#pragma omp parallel private(myid, mybegin, myend, index, iindex, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, count, iistart, end_row)
#pragma omp parallel private(myid, mybegin, myend, index, iindex, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, iistart, end_row)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			index = indexs + myid*A->col;
			iindex = iindexs + myid*col;
			for (i = mybegin; i < myend; ++ i) {
				istart = -1;
				length = 0;
				i1 = i+1;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2) {
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
					}
				}
				// book-keeping [reseting length and setting iistart]
				//count = length;
				iistart = -1;
				//length = 0;
				// use each column that would have resulted from R*A
				//for (j = 0; j < count; ++ j) {
				for (j = 0; j < length; ++ j) {
					jj = istart;
					istart = index[istart];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						if (iindex[N2C(jp[N2C(k)])] == -2) {
							iindex[N2C(jp[N2C(k)])] = iistart;
							iistart = jp[N2C(k)];
							//++length;
						}
					} // end for k
				} // end for j
				// put the correct columns of p into the column list of the products
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					// put the value in B->JA
					jac[N2C(j)] = iistart;
					// set istart to the next value
					iistart = iindex[N2C(iistart)];
					// set the iindex spot to 0
					iindex[N2C(jac[j])] = -2;
				} // end j
			}
		}
		// Third loop: compute entries of R*A*P
		double *acj=(double*)fasp_mem_calloc(iac[row],sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(double);
#endif
		int *BTindexs=(int*)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
		double *temps=(double*)fasp_mem_calloc(A->col*nthreads,sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(double);
#endif
#pragma omp parallel private(myid, mybegin, myend, index, BTindex, temp, i, i1, end_row, j, istart, length, end_rowR, jj, end_rowA, k)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			index = indexs + myid*A->col;
			BTindex = BTindexs + myid*col;
			temp = temps + myid*A->col;
			for (i = mybegin; i < myend; ++ i) {
				i1 = i+1;
				// each col of B
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					BTindex[N2C(jac[N2C(j)])] = j;
				}
				// reset istart and length at the begining of each loop
				istart = -1;
				length = 0;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2){
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
						temp[N2C(ja[N2C(k)])] += rj[N2C(jj)]*aj[N2C(k)];
					}
				}
				// book-keeping [reseting length and setting iistart]
				// use each column that would have resulted from R*A
				for (j = 0; j < length; ++ j) {
					jj = N2C(istart);
					istart = index[N2C(istart)];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
					}
					temp[jj]=0.0;
				}
			}
		}
		// setup coarse matrix B
		B->row = row;
		B->col = col;
		B->IA = iac;
		B->JA = jac;
		B->val = acj;
		B->nnz = B->IA[B->row] - B->IA[0];
		fasp_mem_free(temps);
		fasp_mem_free(indexs);
		fasp_mem_free(iindexs);
		fasp_mem_free(BTindexs);
		fasp_mem_free(part_end);
	}
	else {
		fasp_blas_dcsr_rap (R, A, P, B);
	}
#endif
}

/**
 * \fn void fasp_blas_dcsr_rap2_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, int nthreads, int openmp_holds)
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *B   pointer to dCSRmat matrix equal to R*A*P
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return     void
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 04/25/2011
 */
void fasp_blas_dcsr_rap2_omp (dCSRmat *R, 
												 dCSRmat *A, 
												 dCSRmat *P, 
												 dCSRmat *B,
												 int nthreads,
												 int openmp_holds,
												 int *indexs,
												 int *iindexs,
												 int *BTindexs,
												 double *temps)
{
#if FASP_USE_OPENMP
	if (R->row > openmp_holds) {
		const int row=R->row, col=P->col;
		int *ir=R->IA, *ia=A->IA, *ip=P->IA;
		int *jr=R->JA, *ja=A->JA, *jp=P->JA;
		double *rj=R->val, *aj=A->val, *pj=P->val;
		int istart, iistart;
		int end_row, end_rowA, end_rowR;
		int i, j, jj, k, length, myid, mybegin, myend, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
		int *index = NULL;
		int *iindex = NULL;
		int *BTindex = NULL;
		double *temp = NULL;
		int Acol_mul_nt = A->col*nthreads;
		int col_mul_nt = col*nthreads;
#if 0
		int *indexs=(int *)fasp_mem_calloc(A->col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(int);
#endif
		int *iindexs=(int *)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
#endif
		int *iac=(int*)fasp_mem_calloc(row+1,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (row+1)*sizeof(int);
#endif
		int *part_end=(int*)fasp_mem_calloc(2*nthreads+row,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (2*nthreads+row)*sizeof(int);
#endif
		int *iac_temp=part_end+nthreads;
		fasp_iarray_set_omp(Acol_mul_nt, indexs, -2, nthreads, openmp_holds);
		fasp_iarray_set_omp(col_mul_nt, iindexs, -2, nthreads, openmp_holds);
		// First loop: form iac of R*A*P
#pragma omp parallel private(myid, mybegin, myend, index, iindex, jj_counter, ic, jj_row_begining, end_rowR, jj1, i1, end_rowA, jj2, i2, end_row, jj3, i3)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			index = indexs + myid*A->col;
			iindex = iindexs + myid*col;
			jj_counter = 0;
			for (ic = mybegin; ic < myend; ic ++)
			{
				iindex[ic] = jj_counter;
				jj_row_begining = jj_counter;
				jj_counter ++;
				end_rowR = ir[ic+1];
				for (jj1 = ir[ic]; jj1 < end_rowR; jj1 ++)
				{
					i1 = jr[jj1];
					end_rowA = ia[i1+1];
					for (jj2 = ia[i1]; jj2 < end_rowA; jj2 ++)
					{
						i2 = ja[jj2];
						if (index[i2] != ic)
						{
							index[i2] = ic;
							end_row = ip[i2+1];
							for (jj3 = ip[i2]; jj3 < end_row; jj3 ++)
							{
								i3 = jp[jj3];
								if (iindex[i3] < jj_row_begining)
								{
									iindex[i3] = jj_counter;
									jj_counter ++;
								}
							}
						}
					}
				}
				iac_temp[ic+myid] = jj_row_begining;
			}
			iac_temp[myend+myid] = jj_counter;
			part_end[myid] = myend + myid + 1;
		}
		fasp_iarray_cp_omp(part_end[0], iac_temp, iac, nthreads, openmp_holds);
		jj_counter = part_end[0];
		int Ctemp = 0;
		for (i1 = 1; i1 < nthreads; i1 ++)
		{
			Ctemp += iac_temp[part_end[i1-1]-1];
			for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++)
			{
				iac[jj_counter] = iac_temp[jj1] + Ctemp;
				jj_counter ++;
			}
		}
		//printf("A_H NNZ = %d\n", iac[row]);
		int *jac=(int*)fasp_mem_calloc(iac[row],sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(int);
#endif
		fasp_iarray_set_omp(Acol_mul_nt, indexs, -2, nthreads, openmp_holds);
		fasp_iarray_set_omp(col_mul_nt, iindexs, -2, nthreads, openmp_holds);
		// Second loop: form jac of R*A*P
//#pragma omp parallel private(myid, mybegin, myend, index, iindex, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, count, iistart, end_row)
#pragma omp parallel private(myid, mybegin, myend, index, iindex, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, iistart, end_row)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			index = indexs + myid*A->col;
			iindex = iindexs + myid*col;
			for (i = mybegin; i < myend; ++ i) {
				istart = -1;
				length = 0;
				i1 = i+1;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2) {
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
					}
				}
				// book-keeping [reseting length and setting iistart]
				//count = length;
				iistart = -1;
				//length = 0;
				// use each column that would have resulted from R*A
				//for (j = 0; j < count; ++ j) {
				for (j = 0; j < length; ++ j) {
					jj = istart;
					istart = index[istart];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						if (iindex[N2C(jp[N2C(k)])] == -2) {
							iindex[N2C(jp[N2C(k)])] = iistart;
							iistart = jp[N2C(k)];
							//++length;
						}
					} // end for k
				} // end for j
				// put the correct columns of p into the column list of the products
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					// put the value in B->JA
					jac[N2C(j)] = iistart;
					// set istart to the next value
					iistart = iindex[N2C(iistart)];
					// set the iindex spot to 0
					iindex[N2C(jac[j])] = -2;
				} // end j
			}
		}
		// Third loop: compute entries of R*A*P
		double *acj=(double*)fasp_mem_calloc(iac[row],sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(double);
#endif
#if 0
		int *BTindexs=(int*)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
		double *temps=(double*)fasp_mem_calloc(A->col*nthreads,sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(double);
#endif
#endif
		fasp_iarray_set_omp(col_mul_nt, BTindexs, 0, nthreads, openmp_holds);
		fasp_array_set_omp (Acol_mul_nt, temps, 0.0, nthreads, openmp_holds);
#pragma omp parallel private(myid, mybegin, myend, index, BTindex, temp, i, i1, end_row, j, istart, length, end_rowR, jj, end_rowA, k)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
			index = indexs + myid*A->col;
			BTindex = BTindexs + myid*col;
			temp = temps + myid*A->col;
			for (i = mybegin; i < myend; ++ i) {
				i1 = i+1;
				// each col of B
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					BTindex[N2C(jac[N2C(j)])] = j;
				}
				// reset istart and length at the begining of each loop
				istart = -1;
				length = 0;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2){
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
						temp[N2C(ja[N2C(k)])] += rj[N2C(jj)]*aj[N2C(k)];
					}
				}
				// book-keeping [reseting length and setting iistart]
				// use each column that would have resulted from R*A
				for (j = 0; j < length; ++ j) {
					jj = N2C(istart);
					istart = index[N2C(istart)];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
					}
					temp[jj]=0.0;
				}
			}
		}
		// setup coarse matrix B
		B->row = row;
		B->col = col;
		B->IA = iac;
		B->JA = jac;
		B->val = acj;
		B->nnz = B->IA[B->row] - B->IA[0];
#if 0
		fasp_mem_free(temps);
		fasp_mem_free(indexs);
		fasp_mem_free(iindexs);
		fasp_mem_free(BTindexs);
		fasp_mem_free(part_end);
#endif
	}
	else {
		fasp_blas_dcsr_rap (R, A, P, B);
	}
#endif
}

/**
 * \fn void fasp_blas_dcsr_rap3_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, int *icor_ysk, int nthreads, int openmp_holds)
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *B   pointer to dCSRmat matrix equal to R*A*P
 * \param *icor_ysk pointer to the array
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return     void
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 08/02/2011
 */
void fasp_blas_dcsr_rap3_omp (dCSRmat *R, 
												 dCSRmat *A, 
												 dCSRmat *P, 
												 dCSRmat *B,
												 int *icor_ysk, 
												 int nthreads, 
												 int openmp_holds)
{
#if FASP_USE_OPENMP
	if (R->row > openmp_holds) {
		const int row=R->row, col=P->col;
		int *ir=R->IA, *ia=A->IA, *ip=P->IA;
		int *jr=R->JA, *ja=A->JA, *jp=P->JA;
		double *rj=R->val, *aj=A->val, *pj=P->val;
		int istart, iistart;
		int end_row, end_rowA, end_rowR;
		int i, j, jj, k, length, myid, mybegin, myend, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
		int *index = NULL;
		int *iindex = NULL;
		int *BTindex = NULL;
		double *temp = NULL;
		int FiveMyid, min_A, min_P, A_pos, P_pos, FiveIc;
		int minus_one_length_A = icor_ysk[5*nthreads];
		int minus_one_length_P = icor_ysk[5*nthreads+1];
		int minus_one_length = minus_one_length_A + minus_one_length_P;
		
		//printf(" >>> nnz = %d, row = %d, sparsity: %le\n", A->nnz, A->row, ((double)A->nnz)/A->row/A->row);
		
		int *iindexs = (int *)fasp_mem_calloc(minus_one_length, sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += minus_one_length*sizeof(int);
#endif
		int *indexs = iindexs + minus_one_length_P;

		int *iac=(int*)fasp_mem_calloc(row+1,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (row+1)*sizeof(int);
#endif
		int *part_end=(int*)fasp_mem_calloc(2*nthreads+row,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (2*nthreads+row)*sizeof(int);
#endif
		int *iac_temp=part_end + nthreads;
		int **iindex_array = (int **)fasp_mem_calloc(nthreads, sizeof(int *));
		int **index_array = (int **)fasp_mem_calloc(nthreads, sizeof(int *));
		
		fasp_iarray_set_omp(minus_one_length, iindexs, -2, nthreads, openmp_holds);
#pragma omp parallel for private(myid, FiveMyid, mybegin, myend, min_A, min_P, index, iindex, A_pos, P_pos, ic, FiveIc, jj_counter, jj_row_begining, end_rowR, jj1, i1, end_rowA, jj2, i2, end_row, jj3, i3)
		for (myid = 0; myid < nthreads; myid ++)
		{
			FiveMyid = myid * 5;
			mybegin = icor_ysk[FiveMyid];
			if (myid == nthreads-1) {
				myend = row;
			}
			else {
				myend = icor_ysk[FiveMyid+5];
			}
			min_A = icor_ysk[FiveMyid+2];
			min_P = icor_ysk[FiveMyid+4];
			//printf(" >>> min_A = %d, min_P = %d max_A = %d, max_P = %d, in Proc %d\n", min_A, min_P, min_A+icor_ysk[FiveMyid+1]-1, min_P+icor_ysk[FiveMyid+3]-1, myid);
			A_pos = 0;
			P_pos = 0;
			for (ic = myid-1; ic >= 0; ic --) {
				FiveIc = ic * 5;
				A_pos += icor_ysk[FiveIc+1];
				P_pos += icor_ysk[FiveIc+3];
			}
			iindex_array[myid] = iindex= iindexs + P_pos - min_P;
			index_array[myid] = index = indexs + A_pos - min_A;
			jj_counter = 0;
			for (ic = mybegin; ic < myend; ic ++)
			{
				iindex[ic] = jj_counter;
				jj_row_begining = jj_counter;
				jj_counter ++;
				end_rowR = ir[ic+1];
				for (jj1 = ir[ic]; jj1 < end_rowR; jj1 ++)
				{
					i1 = jr[jj1];
					end_rowA = ia[i1+1];
					for (jj2 = ia[i1]; jj2 < end_rowA; jj2 ++)
					{
						i2 = ja[jj2];
						if (index[i2] != ic)
						{
							index[i2] = ic;
							end_row = ip[i2+1];
							for (jj3 = ip[i2]; jj3 < end_row; jj3 ++)
							{
								i3 = jp[jj3];
								if (iindex[i3] < jj_row_begining)
								{
									iindex[i3] = jj_counter;
									jj_counter ++;
								}
							}
						}
					}
				}
				iac_temp[ic+myid] = jj_row_begining;
			}
			iac_temp[myend+myid] = jj_counter;
			part_end[myid] = myend + myid + 1;
		}
		fasp_iarray_cp_omp(part_end[0], iac_temp, iac, nthreads, openmp_holds);
		jj_counter = part_end[0];
		int Ctemp = 0;
		for (i1 = 1; i1 < nthreads; i1 ++)
		{
			Ctemp += iac_temp[part_end[i1-1]-1];
			for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++)
			{
				iac[jj_counter] = iac_temp[jj1] + Ctemp;
				jj_counter ++;
			}
		}
		int *jac=(int*)fasp_mem_calloc(iac[row],sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(int);
#endif
		fasp_iarray_set_omp(minus_one_length, iindexs, -2, nthreads, openmp_holds);
#pragma omp parallel for private(myid, index, iindex, FiveMyid, mybegin, myend, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, iistart, end_row)
		for (myid = 0; myid < nthreads; myid ++)
		{
			iindex = iindex_array[myid];
			index = index_array[myid];
			FiveMyid = myid * 5;
			mybegin = icor_ysk[FiveMyid];
			if (myid == nthreads-1) {
				myend = row;
			}
			else {
				myend = icor_ysk[FiveMyid+5];
			}
			for (i = mybegin; i < myend; ++ i) {
				istart = -1;
				length = 0;
				i1 = i+1;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2) {
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
					}
				}
				// book-keeping [reseting length and setting iistart]
				//count = length;
				iistart = -1;
				//length = 0;
				// use each column that would have resulted from R*A
				//for (j = 0; j < count; ++ j) {
				for (j = 0; j < length; ++ j) {
					jj = istart;
					istart = index[istart];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						if (iindex[N2C(jp[N2C(k)])] == -2) {
							iindex[N2C(jp[N2C(k)])] = iistart;
							iistart = jp[N2C(k)];
							//++length;
						}
					} // end for k
				} // end for j
				// put the correct columns of p into the column list of the products
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					// put the value in B->JA
					jac[N2C(j)] = iistart;
					// set istart to the next value
					iistart = iindex[N2C(iistart)];
					// set the iindex spot to 0
					iindex[N2C(jac[j])] = -2;
				} // end j
			}
		}
		// Third loop: compute entries of R*A*P
		double *acj=(double*)fasp_mem_calloc(iac[row],sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(double);
#endif
		int *BTindexs=(int*)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += col*nthreads*sizeof(int);
#endif
		double *temps=(double*)fasp_mem_calloc(A->col*nthreads,sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += A->col*nthreads*sizeof(double);
#endif
#pragma omp parallel for private(myid, index, FiveMyid, mybegin, myend, BTindex, temp, i, i1, end_row, j, istart, length, end_rowR, jj, end_rowA, k)
		for (myid = 0; myid < nthreads; myid ++)
		{
			index = index_array[myid];
			FiveMyid = myid * 5;
			mybegin = icor_ysk[FiveMyid];
			if (myid == nthreads-1) {
				myend = row;
			}
			else {
				myend = icor_ysk[FiveMyid+5];
			}
			BTindex = BTindexs + myid*col;
			temp = temps + myid*A->col;
			for (i = mybegin; i < myend; ++ i) {
				i1 = i+1;
				// each col of B
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					BTindex[N2C(jac[N2C(j)])] = j;
				}
				// reset istart and length at the begining of each loop
				istart = -1;
				length = 0;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2){
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
						temp[N2C(ja[N2C(k)])] += rj[N2C(jj)]*aj[N2C(k)];
					}
				}
				// book-keeping [reseting length and setting iistart]
				// use each column that would have resulted from R*A
				for (j = 0; j < length; ++ j) {
					jj = N2C(istart);
					istart = index[N2C(istart)];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
					}
					temp[jj]=0.0;
				}
			}
		}
		// setup coarse matrix B
		B->row = row;
		B->col = col;
		B->IA = iac;
		B->JA = jac;
		B->val = acj;
		B->nnz = B->IA[B->row] - B->IA[0];
		fasp_mem_free(temps);
		fasp_mem_free(iindexs);
		fasp_mem_free(BTindexs);
		fasp_mem_free(part_end);
		fasp_mem_free(iindex_array);
		fasp_mem_free(index_array);
	}
	else {
		fasp_blas_dcsr_rap (R, A, P, B);
	}
#endif
}

/**
 * \fn void fasp_blas_dcsr_rap4_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, int *icor_ysk, int nthreads, int openmp_holds)
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * \param *R   pointer to the dCSRmat matrix
 * \param *A   pointer to the dCSRmat matrix
 * \param *P   pointer to the dCSRmat matrix
 * \param *B   pointer to dCSRmat matrix equal to R*A*P
 * \param *icor_ysk pointer to the array
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return     void
 *
 * Ref. R.E. Bank and C.C. Douglas. SMMP: Sparse Matrix Multiplication Package. 
 *      Advances in Computational Mathematics, 1 (1993), pp. 127-137.
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 08/02/2011
 */
void fasp_blas_dcsr_rap4_omp (dCSRmat *R, 
												 dCSRmat *A, 
												 dCSRmat *P, 
												 dCSRmat *B,
												 int *icor_ysk, 
												 int nthreads, 
												 int openmp_holds)
{
#if FASP_USE_OPENMP
	if (R->row > openmp_holds) {
		const int row=R->row, col=P->col;
		int *ir=R->IA, *ia=A->IA, *ip=P->IA;
		int *jr=R->JA, *ja=A->JA, *jp=P->JA;
		double *rj=R->val, *aj=A->val, *pj=P->val;
		int istart, iistart;
		int end_row, end_rowA, end_rowR;
		int i, j, jj, k, length, myid, mybegin, myend, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
		int *index = NULL;
		int *iindex = NULL;
		int *BTindex = NULL;
		double *temp = NULL;
		int FiveMyid, min_A, min_P, A_pos, P_pos, FiveIc;
		int minus_one_length_A = icor_ysk[5*nthreads];
		int minus_one_length_P = icor_ysk[5*nthreads+1];
		int minus_one_length = minus_one_length_A + minus_one_length_P;
		
		//printf(" >>> nnz = %d, row = %d, sparsity: %le\n", A->nnz, A->row, ((double)A->nnz)/A->row/A->row);
		
		int *iindexs = (int *)fasp_mem_calloc(minus_one_length+minus_one_length_P, sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += minus_one_length*sizeof(int);
#endif
		int *indexs = iindexs + minus_one_length_P;
		int *BTindexs = indexs + minus_one_length_A;

		int *iac=(int*)fasp_mem_calloc(row+1,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (row+1)*sizeof(int);
#endif
		int *part_end=(int*)fasp_mem_calloc(2*nthreads+row,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (2*nthreads+row)*sizeof(int);
#endif
		int *iac_temp=part_end + nthreads;
		int **iindex_array = (int **)fasp_mem_calloc(nthreads, sizeof(int *));
		int **index_array = (int **)fasp_mem_calloc(nthreads, sizeof(int *));
		
		fasp_iarray_set_omp(minus_one_length, iindexs, -2, nthreads, openmp_holds);
#pragma omp parallel for private(myid, FiveMyid, mybegin, myend, min_A, min_P, index, iindex, A_pos, P_pos, ic, FiveIc, jj_counter, jj_row_begining, end_rowR, jj1, i1, end_rowA, jj2, i2, end_row, jj3, i3)
		for (myid = 0; myid < nthreads; myid ++)
		{
			FiveMyid = myid * 5;
			mybegin = icor_ysk[FiveMyid];
			if (myid == nthreads-1) {
				myend = row;
			}
			else {
				myend = icor_ysk[FiveMyid+5];
			}
			min_A = icor_ysk[FiveMyid+2];
			min_P = icor_ysk[FiveMyid+4];
			A_pos = 0;
			P_pos = 0;
			for (ic = myid-1; ic >= 0; ic --) {
				FiveIc = ic * 5;
				A_pos += icor_ysk[FiveIc+1];
				P_pos += icor_ysk[FiveIc+3];
			}
			iindex_array[myid] = iindex= iindexs + P_pos - min_P;
			index_array[myid] = index = indexs + A_pos - min_A;
			jj_counter = 0;
			for (ic = mybegin; ic < myend; ic ++)
			{
				iindex[ic] = jj_counter;
				jj_row_begining = jj_counter;
				jj_counter ++;
				end_rowR = ir[ic+1];
				for (jj1 = ir[ic]; jj1 < end_rowR; jj1 ++)
				{
					i1 = jr[jj1];
					end_rowA = ia[i1+1];
					for (jj2 = ia[i1]; jj2 < end_rowA; jj2 ++)
					{
						i2 = ja[jj2];
						if (index[i2] != ic)
						{
							index[i2] = ic;
							end_row = ip[i2+1];
							for (jj3 = ip[i2]; jj3 < end_row; jj3 ++)
							{
								i3 = jp[jj3];
								if (iindex[i3] < jj_row_begining)
								{
									iindex[i3] = jj_counter;
									jj_counter ++;
								}
							}
						}
					}
				}
				iac_temp[ic+myid] = jj_row_begining;
			}
			iac_temp[myend+myid] = jj_counter;
			part_end[myid] = myend + myid + 1;
		}
		fasp_iarray_cp_omp(part_end[0], iac_temp, iac, nthreads, openmp_holds);
		jj_counter = part_end[0];
		int Ctemp = 0;
		for (i1 = 1; i1 < nthreads; i1 ++)
		{
			Ctemp += iac_temp[part_end[i1-1]-1];
			for (jj1 = part_end[i1-1]+1; jj1 < part_end[i1]; jj1 ++)
			{
				iac[jj_counter] = iac_temp[jj1] + Ctemp;
				jj_counter ++;
			}
		}
		int *jac=(int*)fasp_mem_calloc(iac[row],sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(int);
#endif
		fasp_iarray_set_omp(minus_one_length, iindexs, -2, nthreads, openmp_holds);
#pragma omp parallel for private(myid, index, iindex, FiveMyid, mybegin, myend, i, istart, length, i1, end_rowR, jj, j, end_rowA, k, iistart, end_row)
		for (myid = 0; myid < nthreads; myid ++)
		{
			iindex = iindex_array[myid];
			index = index_array[myid];
			FiveMyid = myid * 5;
			mybegin = icor_ysk[FiveMyid];
			if (myid == nthreads-1) {
				myend = row;
			}
			else {
				myend = icor_ysk[FiveMyid+5];
			}
			for (i = mybegin; i < myend; ++ i) {
				istart = -1;
				length = 0;
				i1 = i+1;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2) {
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
					}
				}
				// book-keeping [reseting length and setting iistart]
				//count = length;
				iistart = -1;
				//length = 0;
				// use each column that would have resulted from R*A
				//for (j = 0; j < count; ++ j) {
				for (j = 0; j < length; ++ j) {
					jj = istart;
					istart = index[istart];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						if (iindex[N2C(jp[N2C(k)])] == -2) {
							iindex[N2C(jp[N2C(k)])] = iistart;
							iistart = jp[N2C(k)];
							//++length;
						}
					} // end for k
				} // end for j
				// put the correct columns of p into the column list of the products
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					// put the value in B->JA
					jac[N2C(j)] = iistart;
					// set istart to the next value
					iistart = iindex[N2C(iistart)];
					// set the iindex spot to 0
					iindex[N2C(jac[j])] = -2;
				} // end j
			}
		}
		// Third loop: compute entries of R*A*P
		double *acj=(double*)fasp_mem_calloc(iac[row],sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(double);
#endif
		double *temps=(double*)fasp_mem_calloc(minus_one_length_A,sizeof(double));
#if CHMEM_MODE
		total_alloc_mem += minus_one_length_A*sizeof(double);
#endif
#pragma omp parallel for private(myid, index, FiveMyid, mybegin, myend, min_A, min_P, A_pos, P_pos, ic, FiveIc, BTindex, temp, i, i1, end_row, j, istart, length, end_rowR, jj, end_rowA, k)
		for (myid = 0; myid < nthreads; myid ++)
		{
			index = index_array[myid];
			FiveMyid = myid * 5;
			mybegin = icor_ysk[FiveMyid];
			if (myid == nthreads-1) {
				myend = row;
			}
			else {
				myend = icor_ysk[FiveMyid+5];
			}
			min_A = icor_ysk[FiveMyid+2];
			min_P = icor_ysk[FiveMyid+4];
			A_pos = 0;
			P_pos = 0;
			for (ic = myid-1; ic >= 0; ic --) {
				FiveIc = ic * 5;
				A_pos += icor_ysk[FiveIc+1];
				P_pos += icor_ysk[FiveIc+3];
			}
			BTindex = BTindexs + P_pos - min_P;
			temp = temps + A_pos - min_A;
			for (i = mybegin; i < myend; ++ i) {
				i1 = i+1;
				// each col of B
				end_row = iac[i1];
				for (j = iac[i]; j < end_row; ++ j) {
					BTindex[N2C(jac[N2C(j)])] = j;
				}
				// reset istart and length at the begining of each loop
				istart = -1;
				length = 0;
				// go across the rows in R
				end_rowR = ir[i1];
				for (jj = ir[i]; jj < end_rowR; ++ jj) {
					j = jr[N2C(jj)];
					// for each column in A
					end_rowA = ia[j+1];
					for (k = ia[j]; k < end_rowA; ++ k) {
						if (index[N2C(ja[N2C(k)])] == -2){
							index[N2C(ja[N2C(k)])] = istart;
							istart = ja[N2C(k)];
							++ length;
						}
						temp[N2C(ja[N2C(k)])] += rj[N2C(jj)]*aj[N2C(k)];
					}
				}
				// book-keeping [reseting length and setting iistart]
				// use each column that would have resulted from R*A
				for (j = 0; j < length; ++ j) {
					jj = N2C(istart);
					istart = index[N2C(istart)];
					index[N2C(jj)] = -2;
					// go across the row of P
					end_row = ip[jj+1];
					for (k = ip[jj]; k < end_row; ++ k) {
						// pull out the appropriate columns of P
						acj[BTindex[N2C(jp[N2C(k)])]]+=temp[jj]*pj[k];
					}
					temp[jj]=0.0;
				}
			}
		}
		// setup coarse matrix B
		B->row = row;
		B->col = col;
		B->IA = iac;
		B->JA = jac;
		B->val = acj;
		B->nnz = B->IA[B->row] - B->IA[0];
		fasp_mem_free(temps);
		fasp_mem_free(iindexs);
		fasp_mem_free(part_end);
		fasp_mem_free(iindex_array);
		fasp_mem_free(index_array);
	}
	else {
		fasp_blas_dcsr_rap (R, A, P, B);
	}
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
