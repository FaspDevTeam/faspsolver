/*! \file smoother.c
 *  \brief Smoothers for sparse matrix in CSR format
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dcsr_jacobi (dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Jacobi method as the smoother in solving Au=b with multigrid method
 *
 * \param u    initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1  the index to begin with
 * \param i_n  the index to end
 * \param s    the step
 * \param *A   pointer to stiffness matrix
 * \param *b   pointer to right hand side
 * \param L    number of iterations
 * \return     void
 *
 * \author     Xuehai Huang, Chensong Zhang
 * \date       09/26/2009
 */
void fasp_smoother_dcsr_jacobi (dvector *u, 
																int i_1, 
																int i_n, 
																int s, 
																dCSRmat *A, 
																dvector *b, 
																int L)
{
	const int N = ABS(i_n - i_1)+1;
	const int *ia=A->IA, *ja=A->JA;
	const double *aj=A->val, *bval=b->val, *uval=u->val;
	int i,j,k,begin_row,end_row;
	
	// Checks should be outside of for; t,d can be allocated before calling!!! --Chensong
	double *t = (double *)fasp_mem_calloc(N,sizeof(double));	
#if CHMEM_MODE
	total_alloc_mem += (N)*sizeof(double);
#endif	
	double *d = (double *)fasp_mem_calloc(N,sizeof(double));
#if CHMEM_MODE
	total_alloc_mem += (N)*sizeof(double);
#endif	
	
	while (L--) {
		
		if (s>0) {
			for (i=i_1;i<=i_n;i+=s) {
				t[i]=bval[i];
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					if (i!=j) t[i]-=aj[k]*uval[j];
					else d[i]=aj[k];
				}				
			}
			
			for (i=i_1;i<=i_n;i+=s) {	
				if (ABS(d[i])>SMALLREAL) u->val[i]=t[i]/d[i];
			}			
		} 
		
		else {
			
			for (i=i_1;i>=i_n;i+=s) {
				t[i]=bval[i];
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					if (i!=j) t[i]-=aj[k]*uval[j];
					else d[i]=aj[k];
				}
				
			}
			
			for (i=i_1;i>=i_n;i+=s) {
				if (ABS(d[i])>SMALLREAL) u->val[i]=t[i]/d[i];
			}
		} 
		
	} // end while
	
	fasp_mem_free(t);
	fasp_mem_free(d);
	
	return;	
}

/**
 * \fn void fasp_smoother_dcsr_gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Gauss-Seidel method  as the smoother in solving Au=b with multigrid method
 *
 * \param u    initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1  the index to begin with
 * \param i_n  the index to end
 * \param s    the step (s=1: forward, s=-1: backward)
 * \param *A   pointer to stiffness matrix
 * \param *b   pointer to right hand side
 * \param L    number of iterations
 * \return     void
 *
 * \author     Xuehai Huang, Chensong Zhang
 * \date       09/26/2009
 */
void fasp_smoother_dcsr_gs (dvector *u, 
														int i_1, 
														int i_n, 
														int s, 
														dCSRmat *A, 
														dvector *b, 
														int L)
{
	const int *ia=A->IA,*ja=A->JA;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	
	int i,j,k,begin_row,end_row;
	double t,d=0.0;
	
	if (s > 0) {
		
		while (L--) {
			for (i=i_1;i<=i_n;i+=s) {
				t = bval[i];
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					if (i!=j) 
						t-=aj[k]*uval[j]; 
					else if (ABS(aj[k])>SMALLREAL) d=1.e+0/aj[k];					
				} // end for k
					uval[i]=t*d;
			} // end for i
		} // end while		
		
	} // if s
	else {		
		
		while (L--) {
			for (i=i_1;i>=i_n;i+=s) {
				t=bval[i];
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					if (i!=j) 
						t-=aj[k]*uval[j];
					else if (ABS(aj[k])>SMALLREAL) d=1e+0/aj[k];
				} // end for k
					uval[i]=t*d;
			} // end for i
		} // end while		
		
  } // end if
	return;
}

/**
 * \fn void fasp_smoother_dcsr_gs_cf (dvector *u, dCSRmat *A, dvector *b, int L, int *mark, int order)
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \param u     initial guess and the new approximation to the solution obtained after L GS steps
 * \param *A    pointer to stiffness matrix
 * \param *b    pointer to right hand side
 * \param L     number of iterations
 * \param *mark C/F marker array
 * \param order C/F ordering: -1: F-first; 1: C-first
 * \return      void
 *
 * \author Zhiyang Zhou
 * \date 2010/11/12 
 */
void fasp_smoother_dcsr_gs_cf (dvector *u, 
															 dCSRmat *A, 
															 dvector *b, 
															 int L, 
															 int *mark, 
															 int order )
{
	const int *ia=A->IA,*ja=A->JA;
	
	int i,j,k,begin_row,end_row;
	int size = b->row;
	
	double *aj=A->val,*bval=b->val,*uval=u->val;
	double t,d=0.0;
	
	if (order == -1) // F-point first
	{	
		while (L--) 
		{
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] != 1)
				{
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d) > SMALLREAL) uval[i] = t/d;
				}
			} // end for i
			
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] == 1)
				{
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d) > SMALLREAL) uval[i] = t/d;
				}
			} // end for i
		} // end while		
	}
	else
	{
		while (L--) 
		{
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] == 1)
				{
					t = bval[i];
					begin_row = ia[i],end_row = ia[i+1]-1;
					for (k = begin_row; k <= end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d) > SMALLREAL) uval[i] = t/d;
				}
			} // end for i
			
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] != 1)
				{
					t = bval[i];
					begin_row = ia[i],end_row = ia[i+1]-1;
					for (k = begin_row; k <= end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d) > SMALLREAL) uval[i] = t/d;
				}
			} // end for i
		} // end while	
	}		
	
	return;
}

/**
 * \fn void fasp_smoother_dcsr_sgs (dvector *u, dCSRmat *A, dvector *b, int L)
 * \brief Symmetric Gauss-Seidel method  as the smoother in solving Au=b with multigrid method
 *
 * \param u    initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param *A   pointer to stiffness matrix
 * \param *b   pointer to right hand side
 * \param L    number of iterations
 * \return     void
 *
 * \author     Xiaozhe Hu
 * \date       10/26/2010
 */
void fasp_smoother_dcsr_sgs (dvector *u, 
														 dCSRmat *A, 
														 dvector *b, 
														 int L)
{
	const int *ia=A->IA,*ja=A->JA;
	const int nm1=b->row-1;
	
	int     i,j,k,begin_row,end_row;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	double  t,d;
	
	while (L--) {
		
		// forward sweep
		for (i=0;i<=nm1;++i) {
			t=bval[i];
			begin_row=ia[i], end_row=ia[i+1];
			for (k=begin_row;k<end_row;++k) {
				j=ja[k];
				if (i!=j) t-=aj[k]*uval[j];
				else d=aj[k];
			} // end for k
			if (ABS(d)>SMALLREAL) uval[i]=t/d;
		} // end for i
		
		// backward sweep
		for (i=nm1-1;i>=0;--i) {
			t=bval[i];
			begin_row=ia[i], end_row=ia[i+1];
			for (k=begin_row;k<end_row;++k) {
				j=ja[k];
				if (i!=j) t-=aj[k]*uval[j];
				else d=aj[k];
			} // end for k
			if (ABS(d)>SMALLREAL) uval[i]=t/d;
		} // end for i
		
	} // end while		
	
	return;
}

/**
 * \fn void fasp_smoother_dcsr_sor (dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L, double w)
 * \brief Successive Overrelaxation method as the smoother in solving Au=b with multigrid method
 *
 * \param u     initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1   the index to begin with
 * \param i_n   the index to end
 * \param s     the step
 * \param *A    pointer to stiffness matrix
 * \param *b    pointer to right hand side
 * \param L     number of iterations
 * \param w     over-relaxation parameter
 * \return      void
 *
 * \author      Xiaozhe Hu
 * \date        10/26/2010
 */
void fasp_smoother_dcsr_sor (dvector *u, 
														 int i_1, 
														 int i_n, 
														 int s, 
														 dCSRmat *A, 
														 dvector *b, 
														 int L, 
														 double w)
{
	const int *ia=A->IA,*ja=A->JA;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	int i,j,k,begin_row,end_row;
	double t, d;
	
	while (L--) {
		if (s>0) {
			for (i=i_1;i<=i_n;i+=s)
			{
				t=bval[i];
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k)
				{
					j=ja[k];
					if (i!=j)
						t-=aj[k]*uval[j];
					else
						d=aj[k];
				}
				if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
			}
		} 
		else {
			for (i=i_1;i>=i_n;i+=s)
			{
				t=bval[i];
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k)
				{
					j=ja[k];
					if (i!=j)
						t-=aj[k]*uval[j];
					else
						d=aj[k];
				}
				if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
			}
		} 
	}  // end while
	
	return;
}

/**
 * \fn void fasp_smoother_dcsr_sor_cf ( dvector *u, dCSRmat *A, dvector *b, int L, double w, int *mark, int order )
 * \brief SOR smoother with C/F ordering for Au=b
 *
 * \param u      initial guess and the new approximation to the solution obtained after L SOR steps
 * \param *A     pointer to stiffness matrix
 * \param *b     pointer to right hand side
 * \param L      number of iterations
 * \param *mark  C/F marker array
 * \param order  C/F ordering: -1: F-first; 1: C-first
 * \return       void
 *
 * \author Zhiyang Zhou
 * \date 2010/11/12 
 */
void fasp_smoother_dcsr_sor_cf (dvector *u, 
																dCSRmat *A, 
																dvector *b, 
																int L, 
																double w, 
																int *mark, 
																int order )
{
	const int *ia=A->IA,*ja=A->JA;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	
	int i,j,k,begin_row,end_row;
	int size = b->row;
	double t,d=0.0;
	
	if (order == -1) // F-point first
	{	
		while (L--) 
		{
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] == 0 || mark[i] == 2)
				{
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
				}
			} // end for i
			
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] == 1)
				{
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
				}
			} // end for i
		} // end while		
	}
	else
	{
		while (L--) 
		{
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] == 1)
				{
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
				}
			} // end for i
			
			for (i = 0; i < size; i ++) 
			{
				if (mark[i] != 1)
				{
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aj[k]*uval[j]; 
						else d = aj[k];					
					} // end for k
					if (ABS(d)>SMALLREAL) uval[i]=w*(t/d)+(1-w)*uval[i];
				}
			} // end for i
		} // end while	
	}		
	
	return;
}

/**
 * \fn void fasp_smoother_dcsr_ilu (dCSRmat *A, dvector *b, dvector *x, void *data)
 * \brief ILU method as the smoother in solving Au=b with multigrid method
 *
 * \param *A    pointer to stiffness matrix
 * \param *b    pointer to right hand side
 * \param *x    pointer to current solution
 * \param *data pointer to user defined data
 * \return void
 *
 * \author Shiquan Zhang, Xiaozhe Hu
 * \date 2010/11/12 
 */
void fasp_smoother_dcsr_ilu (dCSRmat *A, 
														 dvector *b, 
														 dvector *x, 
														 void *data)
{
	ILU_data *iludata=(ILU_data *)data;
	const unsigned int m=A->row, m2=2*m, memneed=3*m;
	double *zz, *zr, *z;
	
	if (iludata->nwork<memneed) goto MEMERR; 
	
	zz = iludata->work; 
	zr = iludata->work+m;
	z  = iludata->work+m2;
	
	{
		int i, j, jj, begin_row, end_row;
		double *lu = iludata->luval;
		int *ijlu = iludata->ijlu;
		double *xval = x->val, *bval = b->val;
		
		/** form residual zr = b - A x */
		fasp_array_cp(m,bval,zr); fasp_blas_dcsr_aAxpy(-1.0,A,xval,zr);
		
		// forward sweep: solve unit lower matrix equation L*zz=zr
		zz[0]=zr[0];
		for (i=1;i<m;++i) {
			begin_row=ijlu[i]; end_row=ijlu[i+1];
			for (j=begin_row;j<end_row;++j) {
				jj=ijlu[j];
				if (jj<i) zr[i]-=lu[j]*zz[jj];
				else break;
			}
			zz[i]=zr[i];
		}
		
		// backward sweep: solve upper matrix equation U*z=zz
		z[m-1]=zz[m-1]*lu[m-1];
		for (i=m-2;i>=0;--i) {
			begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
			for (j=end_row;j>=begin_row;--j) {
				jj=ijlu[j];
				if (jj>i) zz[i]-=lu[j]*z[jj];
				else break;
			} 
			z[i]=zz[i]*lu[i];
		}
		
		fasp_blas_array_axpy(m,1,z,xval);		
	}
	
	return;
	
MEMERR:
	printf("Error: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
	exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_smoother_dcsr_kaczmarz (dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief kaczmarz method as the smoother for solving Au=b
 *
 * \param u    initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1  the index to begin with
 * \param i_n  the index to end
 * \param s    the step (s=1: forward, s=-1: backward)
 * \param *A   pointer to stiffness matrix
 * \param *b   pointer to right hand side
 * \param L    number of iterations
 * \return void
 *
 * \author Xiaozhe Hu
 * \date 2010/11/12 
 */
void fasp_smoother_dcsr_kaczmarz (dvector *u, 
																	int i_1, 
																	int i_n, 
																	int s, 
																	dCSRmat *A, 
																	dvector *b, 
																	int L, 
																	double w)
{
	const int *ia=A->IA,*ja=A->JA;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	
	int i,j,k,begin_row,end_row;
	double temp1,temp2,alpha;
	
	if (s > 0) {
		
		while (L--) {
			for (i=i_1;i<=i_n;i+=s) {
				temp1 = 0; temp2 = 0;
				
				begin_row=ia[i], end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					temp1 += aj[k]*aj[k];
					temp2 += aj[k]*uval[j];
				} // end for k
				
				alpha = (bval[i] - temp2)/temp1;
				
				for (k=begin_row;k<end_row;++k){
					j = ja[k];
					uval[j] += w*alpha*aj[k];
				}// end for k
			} // end for i
		} // end while	
		
	} // if s
	
	else {		
		
		while (L--) {
			for (i=i_1;i>=i_n;i+=s) {
				temp1 = 0; temp2 = 0;
				
				begin_row=ia[i], end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					temp1 += aj[k]*aj[k];
					temp2 += aj[k]*uval[j];
				} // end for k
				
				alpha = (bval[i] - temp2)/temp1;
				
				for (k=begin_row;k<end_row;++k){
					j = ja[k];
					uval[j] += w*alpha*aj[k];
				}// end for k
			} // end for i
		} // end while		
		
	} // end if
	
	return;
}

/**
 * \fn void fasp_smoother_dcsr_L1diag (dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Diagonal scaling (using L1 norm) as the smoother in solving Au=b with multigrid method
 *
 * \param u    initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1  the index to begin with
 * \param i_n  the index to end
 * \param s    the step
 * \param *A   pointer to stiffness matrix
 * \param *b   pointer to right hand side
 * \param L    number of iterations
 * \return void
 *
 * \author Xiaozhe Hu, James Brannick
 * \date 01/26/2011
 */
void fasp_smoother_dcsr_L1diag (dvector *u, 
																int i_1, 
																int i_n, 
																int s, 
																dCSRmat *A, 
																dvector *b, 
																int L)
{
	const int N = ABS(i_n - i_1)+1;
	const int *ia=A->IA, *ja=A->JA;
	const double *aj=A->val, *bval=b->val, *uval=u->val;
	int i,j,k,begin_row,end_row;
	
	// Checks should be outside of for; t,d can be allocated before calling!!! --Chensong
	double *t = (double *)fasp_mem_calloc(N,sizeof(double));	
#if CHMEM_MODE
	total_alloc_mem += (N)*sizeof(double);
#endif	
	double *d = (double *)fasp_mem_calloc(N,sizeof(double));
#if CHMEM_MODE
	total_alloc_mem += (N)*sizeof(double);
#endif	
	
	while (L--) {
		
		if (s>0) {
			for (i=i_1;i<=i_n;i+=s) {
				t[i]=bval[i]; d[i]=0.0;
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					t[i]-=aj[k]*uval[j];
					d[i]+=ABS(aj[k]);
				}				
			}
			
			for (i=i_1;i<=i_n;i+=s) {	
				if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
			}			
		} 
		
		else {
			
			for (i=i_1;i>=i_n;i+=s) {
				t[i]=bval[i];d[i]=0.0;
				begin_row=ia[i],end_row=ia[i+1];
				for (k=begin_row;k<end_row;++k) {
					j=ja[k];
					t[i]-=aj[k]*uval[j];
					d[i]+=ABS(aj[k]);
				}
				
			}
			
			for (i=i_1;i>=i_n;i+=s) {
				if (ABS(d[i])>SMALLREAL) u->val[i]+=t[i]/d[i];
			}
		} 
		
	} // end while
	
	fasp_mem_free(t);
	fasp_mem_free(d);
	
	return;	
}

#if 0
/**
 * \fn dCSRmat static form_contractor(dCSRmat *A, int smoother, int steps, int ndeg, double relax, double dtol)
 * \brief form contractor I-BA
 *
 * \param A          pointer to the dCSRmat
 * \param smoother   the smoother type
 * \param steps      smoothing steps
 * \param ndeg       degree of the polynomial smoother
 * \param relax      relaxation parameter for SOR smoother
 * \param dtol       drop tplerance for droping small entries in matrix
 * \return dCSRmat 
 *
 * \author Xiaozhe Hu, James Brannick
 * \date   01/26/2011
 *
 * \note: not a O(N) algorithm, need to be modified!!!!
 */
static dCSRmat form_contractor(dCSRmat *A, int smoother, int steps, int ndeg, double relax, double dtol)
{
	const int n=A->row;
	unsigned int i;
	
	double *work= (double *)fasp_mem_calloc(2*n,sizeof(double));
	
#if CHMEM_MODE
	total_alloc_mem += (2*n)*sizeof(double);
#endif	
	
	dvector b, x;
	b.row=x.row=n;
	b.val=work; x.val=work+n;
	
	int *index=(int *)fasp_mem_calloc(n,sizeof(int));	
	
#if CHMEM_MODE
	total_alloc_mem += (n)*sizeof(int);
#endif
	
	for (i=0; i<n; ++i) index[i]=i;
	
	dCSRmat B=fasp_dcsr_create(n, n, n*n); // too much memory required, need to change!!
	dCSRmat C, D;
	
	for (i=0; i<n; ++i){
		
		// get i-th column
		fasp_dcsr_getcol(i, A, b.val);
		
		// set x =0.0 
		fasp_dvec_set(n, &x, 0.0);
		
		// smooth
		switch (smoother) {
			case GS:
				fasp_smoother_dcsr_gs(&x, 0, n-1, 1, A, &b, steps);
				break;
			case POLY:
				fasp_smoother_dcsr_poly(A, &b, &x, n, ndeg, steps); 
				break;
			case JACOBI:
				fasp_smoother_dcsr_jacobi(&x, 0, n-1, 1, A, &b, steps);
				break;
			case SGS:
				fasp_smoother_dcsr_sgs(&x, A, &b, steps);
				break;
			case SOR:
				fasp_smoother_dcsr_sor(&x, 0, n-1, 1, A, &b, steps, relax);
				break;
			case SSOR:
				fasp_smoother_dcsr_sor(&x, 0, n-1, 1, A, &b, steps, relax);
				fasp_smoother_dcsr_sor(&x, n-1, 0,-1, A, &b, steps, relax);
				break;
			case GSOR:
				fasp_smoother_dcsr_gs(&x, 0, n-1, 1, A, &b, steps);
				fasp_smoother_dcsr_sor(&x, n-1, 0, -1, A, &b, steps, relax);
				break;
			case SGSOR:
				fasp_smoother_dcsr_gs(&x, 0, n-1, 1, A, &b, steps);
				fasp_smoother_dcsr_gs(&x, n-1, 0,-1, A, &b, steps);
				fasp_smoother_dcsr_sor(&x, 0, n-1, 1, A, &b, steps, relax);
				fasp_smoother_dcsr_sor(&x, n-1, 0,-1, A, &b, steps, relax);
				break;
			default:
				printf("Error: wrong smoother type!\n"); exit(ERROR_INPUT_PAR);
		} 
		
		// store to B
		B.IA[i] = i*n;
		memcpy(&(B.JA[i*n]), index, n*sizeof(int));
		memcpy(&(B.val[i*n]), x.val, x.row*sizeof(double));
		
	}
	
	B.IA[n] = n*n;
	
	// drop small entries
	compress_dCSRmat(&B, &D, dtol);
	
	// get contractor
	fasp_dcsr_trans(&D, &C);
	
	// clean up
	fasp_mem_free(work);
	fasp_dcsr_free(&B);
	fasp_dcsr_free(&D);
	
	return C;
}
#endif

/*-----------------------------------omp--------------------------------------*/

/**
 * \fn void fasp_smoother_gs_cf_omp( dvector *u, dCSRmat *A, dvector *b, int L, int *mark, int order, int nthreads, int openmp_holds )
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \param u     initial guess and the new approximation to the solution obtained after L GS steps
 * \param *A    pointer to stiffness matrix
 * \param *b    pointer to right hand side
 * \param L     number of iterations
 * \param *mark C/F marker array
 * \param order C/F ordering: -1: F-first; 1: C-first
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return      void
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
void fasp_smoother_dcsr_gs_cf_omp (dvector *u, 
																	 dCSRmat *A, 
																	 dvector *b, 
																	 int L, 
																	 int *mark, 
																	 int order, 
																	 int nthreads, 
																	 int openmp_holds)
{
#if FASP_USE_OPENMP
	const int *ia=A->IA,*ja=A->JA;
	
	int i,j,k,begin_row,end_row;
	int size = b->row;
	
	double *aj=A->val,*bval=b->val,*uval=u->val;
	double t,d=0.0;
	int myid, mybegin, myend;
	
	if (order == -1) // F-point first
	{	
		while (L--) 
		{
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] != 1)
						{
							t = bval[i];
							begin_row = ia[i], end_row = ia[i+1];
							for (k = begin_row; k < end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] != 1)
					{
						t = bval[i];
						begin_row = ia[i], end_row = ia[i+1];
						for (k = begin_row; k < end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
			
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] == 1)
						{
							t = bval[i];
							begin_row = ia[i], end_row = ia[i+1];
							for (k = begin_row; k < end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];					
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] == 1)
					{
						t = bval[i];
						begin_row = ia[i], end_row = ia[i+1];
						for (k = begin_row; k < end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
		} // end while		
	}
	else
	{
		while (L--) 
		{
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] == 1)
						{
							t = bval[i];
							begin_row = ia[i],end_row = ia[i+1]-1;
							for (k = begin_row; k <= end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] == 1)
					{
						t = bval[i];
						begin_row = ia[i],end_row = ia[i+1]-1;
						for (k = begin_row; k <= end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
			
			if (size > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d) //num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] != 1)
						{
							t = bval[i];
							begin_row = ia[i],end_row = ia[i+1]-1;
							for (k = begin_row; k <= end_row; k ++) 
							{
								j = ja[k];
								if (i!=j) t -= aj[k]*uval[j]; 
								else d = aj[k];
							} // end for k
							if (ABS(d) > SMALLREAL) uval[i] = t/d;
						}
					} // end for i
				}
			}
			else {
				for (i = 0; i < size; i ++) 
				{
					if (mark[i] != 1)
					{
						t = bval[i];
						begin_row = ia[i],end_row = ia[i+1]-1;
						for (k = begin_row; k <= end_row; k ++) 
						{
							j = ja[k];
							if (i!=j) t -= aj[k]*uval[j]; 
							else d = aj[k];					
						} // end for k
						if (ABS(d) > SMALLREAL) uval[i] = t/d;
					}
				} // end for i
			}
		} // end while	
	}		
#endif
	return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
