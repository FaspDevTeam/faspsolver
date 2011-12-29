/*! \file blas_csr_omp.c
 *  \brief BLAS operations for sparse matrices in CSR format.
 *
 *  \note Sparse functions usually contain three runs. The three runs are all the 
 *        same but thy serve different purpose. 
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
/*-----------------------------------omp--------------------------------------*/
/* Feng Chunsheng Yue Xiaoqiang  Expand the inner loop  /mar/14/2011/  */
/**
 * \fn void fasp_blas_dcsr_mxv_omp (dCSRmat *A, REAL *x, REAL *y, INT nthreads, 
 *                                  INT openmp_holds)
 * \brief Matrix-vector multiplication y = A*x
 * \param *A pointer to dCSRmat CSR matrix
 * \param *x pointer to dvector
 * \param *y pointer to dvector
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/14/2011
 */
void fasp_blas_dcsr_mxv_omp (dCSRmat *A, 
                             REAL *x, 
                             REAL *y, 
                             INT nthreads, 
                             INT openmp_holds)
{
#if FASP_USE_OPENMP
	const INT m=A->row;
	const INT *ia=A->IA, *ja=A->JA;
	const REAL *aj=A->val;
	INT i, k, begin_row, end_row, nnz_num_row;
	register REAL temp;
	
	if (m > openmp_holds) {
		INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, temp, begin_row, end_row, nnz_num_row, k) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
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
 * \fn void fasp_blas_dcsr_aAxpy_omp (const REAL alpha, dCSRmat *A, REAL *x,
 *                                    REAL *y, INT nthreads, INT openmp_holds)
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
void fasp_blas_dcsr_aAxpy_omp (const REAL alpha, 
                               dCSRmat *A, 
                               REAL *x, 
                               REAL *y, 
                               INT nthreads, 
                               INT openmp_holds)
{
#if FASP_USE_OPENMP
	const INT  m  = A->row;
	const INT *ia = A->IA, *ja = A->JA;
	const REAL *aj = A->val;
	INT i, k, begin_row, end_row;
	register REAL temp;
	
	if ( alpha == 1.0 ) {
		if (m > openmp_holds) {
			INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, temp, begin_row, end_row, k) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
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
			INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, temp, begin_row, end_row, k) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
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
			INT myid, mybegin, myend;
#pragma omp parallel private(myid, mybegin, myend, i, temp, begin_row, end_row, k) ////num_threads(nthreads)
			{
				myid = omp_get_thread_num();
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
 * \fn REAL fasp_blas_dcsr_vmv_omp (dCSRmat *A, REAL *x, REAL *y, 
 *                                    INT nthreads, INT openmp_holds)
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
REAL fasp_blas_dcsr_vmv_omp (dCSRmat *A, 
                             REAL *x, 
                             REAL *y, 
                             INT nthreads, 
                             INT openmp_holds)
{
	register REAL value=0.0;
#if FASP_USE_OPENMP
	const INT m=A->row;
	const INT *ia=A->IA, *ja=A->JA;
	const REAL *aj=A->val;
	INT i, k, begin_row, end_row;
	register REAL temp;
	
	if (m > openmp_holds) {
#pragma omp parallel for reduction(+:value) private(i,temp,begin_row,end_row,k) ////num_threads(nthreads)
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
 * \fn void fasp_blas_dcsr_rap_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP, INT nthreads, INT openmp_holds)
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
void fasp_blas_dcsr_rap_omp (dCSRmat  *R,
                             dCSRmat  *A,
                             dCSRmat  *P,
                             dCSRmat  *RAP,
                             INT       nthreads,
                             INT       openmp_holds )
{
#if FASP_USE_OPENMP
    INT n_coarse = R->row;
    INT *R_i = R->IA;
    INT *R_j = R->JA;
    REAL *R_data = R->val;
    
    INT n_fine = A->row;
    INT *A_i = A->IA;
    INT *A_j = A->JA;
    REAL *A_data = A->val;
    
    INT *P_i = P->IA;
    INT *P_j = P->JA;
    REAL *P_data = P->val;
    
    INT RAP_size;
    INT *RAP_i = NULL;
    INT *RAP_j = NULL;
    REAL *RAP_data = NULL;
    
    INT *P_marker = NULL;
    INT *A_marker = NULL;
    
    INT *Ps_marker = NULL;
    INT *As_marker = NULL;
    INT *RAP_temp = NULL;
    INT *part_end = NULL;
    
    INT ic, i1, i2, i3, jj1, jj2, jj3;
    INT jj_counter, jj_row_begining;
    REAL r_entry, r_a_product, r_a_p_product;
    
    // Save the memory
    if (n_coarse <= openmp_holds)
    {
        nthreads = 1;
    }
    
    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;
    
    Ps_marker = (INT *)fasp_mem_calloc(total_calloc, sizeof(int));
    As_marker = Ps_marker + coarse_mul_nthreads;
    RAP_temp = As_marker + fine_mul_nthreads;
    part_end = RAP_temp + coarse_add_nthreads;
    
    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
    RAP_i = (INT *)fasp_mem_calloc(n_coarse+1, sizeof(int));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        INT myid, mybegin, myend, Ctemp;
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
    
    RAP_j = (INT *)fasp_mem_calloc(RAP_size, sizeof(int));
    RAP_data = (REAL *)fasp_mem_calloc(RAP_size, sizeof(REAL));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        INT myid, mybegin, myend;
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
 * \fn void fasp_blas_dcsr_rap_agg_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *RAP, INT nthreads, INT openmp_holds)
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
void fasp_blas_dcsr_rap_agg_omp (dCSRmat  *R,
                                 dCSRmat  *A,
                                 dCSRmat  *P,
                                 dCSRmat  *RAP,
                                 INT       nthreads,
                                 INT       openmp_holds )
{
#if FASP_USE_OPENMP
    INT n_coarse = R->row;
    INT *R_i = R->IA;
    INT *R_j = R->JA;
    
    INT n_fine = A->row;
    INT *A_i = A->IA;
    INT *A_j = A->JA;
    REAL *A_data = A->val;
    
    INT *P_i = P->IA;
    INT *P_j = P->JA;
    
    INT RAP_size;
    INT *RAP_i = NULL;
    INT *RAP_j = NULL;
    REAL *RAP_data = NULL;
    
    INT *P_marker = NULL;
    INT *A_marker = NULL;
    
    INT *Ps_marker = NULL;
    INT *As_marker = NULL;
    INT *RAP_temp = NULL;
    INT *part_end = NULL;
    
    INT ic, i1, i2, i3, jj1, jj2, jj3;
    INT jj_counter, jj_row_begining;
    REAL a_entry;
    
    // Save the memory
    if (n_coarse <= openmp_holds)
    {
        nthreads = 1;
    }
    
    INT coarse_mul_nthreads = n_coarse * nthreads;
    INT fine_mul_nthreads = n_fine * nthreads;
    INT coarse_add_nthreads = n_coarse + nthreads;
    INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    INT total_calloc = minus_one_length + coarse_add_nthreads + nthreads;
    
    Ps_marker = (INT *)fasp_mem_calloc(total_calloc, sizeof(int));
    As_marker = Ps_marker + coarse_mul_nthreads;
    RAP_temp = As_marker + fine_mul_nthreads;
    part_end = RAP_temp + coarse_add_nthreads;
    
    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_i  *
     *------------------------------------------------------*/
    RAP_i = (INT *)fasp_mem_calloc(n_coarse+1, sizeof(int));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        INT myid, mybegin, myend, Ctemp;
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
    
    RAP_j = (INT *)fasp_mem_calloc(RAP_size, sizeof(int));
    RAP_data = (REAL *)fasp_mem_calloc(RAP_size, sizeof(REAL));
    
    fasp_iarray_set_omp(minus_one_length, Ps_marker, -1, nthreads, openmp_holds);
    
    if (n_coarse > openmp_holds)
    {
        INT myid, mybegin, myend;
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
 * \fn void fasp_blas_dcsr_rap1_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, INT nthreads, INT openmp_holds)
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
                              INT nthreads,
                              INT openmp_holds)
{
#if FASP_USE_OPENMP
	if (R->row > openmp_holds) {
		const INT row=R->row, col=P->col;
		INT *ir=R->IA, *ia=A->IA, *ip=P->IA;
		INT *jr=R->JA, *ja=A->JA, *jp=P->JA;
		REAL *rj=R->val, *aj=A->val, *pj=P->val;
		INT istart, iistart;
		INT end_row, end_rowA, end_rowR;
		INT i, j, jj, k, length, myid, mybegin, myend, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
		INT *index = NULL;
		INT *iindex = NULL;
		INT *BTindex = NULL;
		REAL *temp = NULL;
		
		INT *indexs=(INT *)fasp_mem_calloc(A->col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(int);
#endif
		INT *iindexs=(INT *)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
		INT *iac=(int*)fasp_mem_calloc(row+1,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (row+1)*sizeof(int);
#endif
		INT *part_end=(int*)fasp_mem_calloc(2*nthreads+row,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (2*nthreads+row)*sizeof(int);
#endif
		INT *iac_temp=part_end+nthreads;
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
		INT Ctemp = 0;
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
		INT *jac=(int*)fasp_mem_calloc(iac[row],sizeof(int));
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
		REAL *acj=(REAL*)fasp_mem_calloc(iac[row],sizeof(REAL));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(REAL);
#endif
		INT *BTindexs=(int*)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
		REAL *temps=(REAL*)fasp_mem_calloc(A->col*nthreads,sizeof(REAL));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(REAL);
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
 * \fn void fasp_blas_dcsr_rap2_omp (dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B, INT nthreads, INT openmp_holds)
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
                              INT nthreads,
                              INT openmp_holds,
                              INT *indexs,
                              INT *iindexs,
                              INT *BTindexs,
                              REAL *temps)
{
#if FASP_USE_OPENMP
	if (R->row > openmp_holds) {
		const INT row=R->row, col=P->col;
		INT *ir=R->IA, *ia=A->IA, *ip=P->IA;
		INT *jr=R->JA, *ja=A->JA, *jp=P->JA;
		REAL *rj=R->val, *aj=A->val, *pj=P->val;
		INT istart, iistart;
		INT end_row, end_rowA, end_rowR;
		INT i, j, jj, k, length, myid, mybegin, myend, jj_counter, ic, jj_row_begining, jj1, i1, jj2, i2, jj3, i3;
		INT *index = NULL;
		INT *iindex = NULL;
		INT *BTindex = NULL;
		REAL *temp = NULL;
		INT Acol_mul_nt = A->col*nthreads;
		INT col_mul_nt = col*nthreads;
#if 0
		INT *indexs=(INT *)fasp_mem_calloc(A->col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(int);
#endif
		INT *iindexs=(INT *)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
#endif
		INT *iac=(int*)fasp_mem_calloc(row+1,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (row+1)*sizeof(int);
#endif
		INT *part_end=(int*)fasp_mem_calloc(2*nthreads+row,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (2*nthreads+row)*sizeof(int);
#endif
		INT *iac_temp=part_end+nthreads;
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
		INT Ctemp = 0;
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
		INT *jac=(int*)fasp_mem_calloc(iac[row],sizeof(int));
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
		REAL *acj=(REAL*)fasp_mem_calloc(iac[row],sizeof(REAL));
#if CHMEM_MODE
		total_alloc_mem += iac[row]*sizeof(REAL);
#endif
#if 0
		INT *BTindexs=(int*)fasp_mem_calloc(col*nthreads,sizeof(int));
#if CHMEM_MODE
		total_alloc_mem += (col*nthreads)*sizeof(int);
#endif
		REAL *temps=(REAL*)fasp_mem_calloc(A->col*nthreads,sizeof(REAL));
#if CHMEM_MODE
		total_alloc_mem += (A->col*nthreads)*sizeof(REAL);
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
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
