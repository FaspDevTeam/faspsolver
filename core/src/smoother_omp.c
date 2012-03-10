/*! \file smoother_omp.c
 *  \brief Smoothers for sparse matrix in CSR format
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/*
 The following code is based on SiPSMG (Simple Poisson Solver based on MultiGrid)
 
 (c) 2008 Johannes Kraus, Jinchao Xu, Yunrong Zhu, Ludmil Zikatanov
 */

/*==========================================================================*/
/* swep3db (backward)                                                       */
/*==========================================================================*/
void swep3db(int *ia, 
             int *ja, 
             double *aa,
             double *u, 
             double *f,
             int nbegx,
             int nbegy,
             int nbegz,
             int *mark,
             int nx, int ny, int nz)
{
    int nxy,k,k0,j,j0,i,i0;
    int begin_row,end_row,ii,jj;
    double t,d;
    
    nxy=nx*ny;
    nbegx = nx + nbegx;
    nbegy = ny + nbegy;
    nbegz = nz + nbegz;
    
    for (k=nbegz; k >=0; k-=2)
    {
        k0= k*nxy;
        for (j = nbegy; j >=0; j-=2)
        {
            j0= j*nx;
            
            for (i = nbegx; i >=0; i-=2)			
            {
                i0 = i   +  j0    + k0;
                i0 = mark[i0]-1;  //Fortran to C
                if (i0>=0 ) {
                    t = f[i0];
                    begin_row = ia[i0], end_row = ia[i0+1];
                    for (ii = begin_row; ii < end_row; ii ++) 
                    {
                        jj = ja[ii];
                        if (i0!=jj) t -= aa[ii]*u[jj]; 
                        else d = aa[ii];					
                    } // end for ii
                    
                    if (ABS(d) > SMALLREAL) u[i0] = t/d;
                } //if (i0>=0 ) 
            }
        }
    }
}

/*==========================================================================*/
/* rb0b3d                                                                   */
/*==========================================================================*/
/**
 * \fn rb0b3d(int *ia, int *ja, double *aa,double *u, double *f, int *mark, int nx, int ny, int nz, int nsweeps)
 * \brief Colores Gauss-Seidel smoother  for Au=b
 *  fix the bug ,let nx,ny,nz fit any even or odd.
 * * \param mark     maping from geometry number to algebraical dof number
 * \author FENG Chunsheng
 * \date 2012/Feb/06 
 */
void rb0b3d(int *ia, int *ja, double *aa,double *u, double *f, int *mark, int nx, int ny, int nz, int nsweeps)
{
    /*
     n1 = n+1
     ... Backward GAUSS - SEIDEL. ON INPUT THE INITIAL GUESS IS IN U
     ... On OUT the new iterate is also in U
     ...  OUTPUT AND THE DIFFERENCE BETWEEN THE INITIAL GUESS AND
     ...  THIS ITERATE IS IN U
     ...  the ordering is as follows (e for even o for odd):
     ...  e-e-e
     ...  e-e-o
     ...  e-o-e
     ...  e-o-o
     ...  o-e-e
     ...  o-e-o
     ...  o-o-e
     ...  o-o-o
     ... There are 8 loops that look the same.
     */
    int n0e,n0o,isweep;
    int ex,ey, ez;
    int ox,oy, oz;
        
    n0e= -1;
    n0o= -2;
    
    for (isweep = 1; isweep <= nsweeps; isweep++)
    {
        if ((nx%2==0) &&(ny%2 ==0)  &&(nz%2==0)) {          //000
            /*...  e-e-e (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
        } else if ((nx%2==0) &&(ny%2 ==0)  &&(nz%2==1)) {   //001
            /*...  e-e-o (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
        }  else if ((nx%2==0)&&(ny%2 ==1)&&(nz%2==0)) {     //010
            /*...  e-o-e (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);	
        }  else if ((nx%2==0)&&(ny%2 ==1)&&(nz%2==1)) {     //011
            /*...  e-o-o (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);	
        }  else if ((nx%2==1)&&(ny%2 ==0)&&(nz%2==0)) {     //100
            /*...  o-e-e (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);	
        }  else if ((nx%2==1)&&(ny%2 ==0)&&(nz%2==1)) {   //101
            /*...  o-e-o (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);	
        }  else if ((nx%2==1)&&(ny%2 ==1)&&(nz%2==0)) {     //110
            /*...  o-o-e (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  o-o-o */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);	
        }  else if ((nx%2==1)&&(ny%2 ==1)&&(nz%2==1)) {   //111
            /*...  o-o-o (and going backwards) */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
            /*...  o-o-e */
            swep3db(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
            /*...  o-e-o */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
            /*...  o-e-e */
            swep3db(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
            /*...  e-o-o */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
            /*...  e-o-e */
            swep3db(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
            /*...  e-e-o */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
            /*...  e-e-e */
            swep3db(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);	
        }
    }
}


/*==========================================================================*/
/* swep3df (forward)                                                        */
/*==========================================================================*/
void swep3df(int *ia, 
             int *ja, 
             double *aa,
             double *u, 
             double *f,
             int nbegx, 
             int nbegy, 
             int nbegz, 
             int *mark,
             int nx, int ny, int nz)
{
    int nxy,k,k0,j,j0,i,i0;
    int begin_row,end_row,ii,jj;
    double t,d;
    nxy=nx*ny;
    for (k=nbegz; k < nz; k+=2)
    {
        k0= k*nxy;
        
        for (j = nbegy; j < ny; j+=2)
        {
            j0= j*nx;
            
            for (i = nbegx; i < nx; i+=2)			/*!*/
            {
                i0 = i   +  j0    + k0;
                i0 = mark[i0]-1; //Fortran to C
                
                //		printf("%d %d %d %d\n",i,j0,k0,i0);
                if (i0>=0 ) {
                    
                    t = f[i0];
                    begin_row = ia[i0], end_row = ia[i0+1];
                    for (ii = begin_row; ii < end_row; ii ++) 
                    {
                        jj = ja[ii];
                        if (i0!=jj) t -= aa[ii]*u[jj]; 
                        else d = aa[ii];					
                    } // end for ii
                    
                    if (ABS(d) > SMALLREAL) u[i0] = t/d;
                } //	if (i0>=0 ) 
            }
        }
    }
    
}

/*==========================================================================*/
/* rb0f3d                                                                   */
/*==========================================================================*/
void rb0f3d(int *ia, int *ja, double *aa,double *u, double *f, int *mark, int nx, int ny, int nz, int nsweeps)
{
    /*
     n1 = n+1
     ... FORWARD GAUSS - SEIDEL. ON INPUT THE INITIAL GUESS IS IN U
     ... THE NEW ITERATE IS IN V ON
     ...  OUTPUT AND THE DIFFERENCE BETWEEN THE INITIAL GUESS AND
     ...  THIS ITERATE IS IN U
     ...  the ordering is as follows (e for even o for odd):
     ...  o-o-o
     ...  o-o-e
     ...  o-e-o
     ...  o-e-e
     ...  e-o-o
     ...  e-o-e
     ...  e-e-o
     ...  e-e-e
     */
    int n0e,n0o,isweep;
    
    n0e=0;
    n0o=1;
    
    
    for (isweep = 1; isweep <= nsweeps; isweep++)
    {
        /*...  o-o-o */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0o,mark,nx,ny,nz);
        
        /*...  o-o-e */
        swep3df(ia,ja,aa,u,f,n0o,n0o,n0e,mark,nx,ny,nz);
        
        /*...  o-e-o */
        swep3df(ia,ja,aa,u,f,n0o,n0e,n0o,mark,nx,ny,nz);
        
        /*...  o-e-e */
        swep3df(ia,ja,aa,u,f,n0o,n0e,n0e,mark,nx,ny,nz);
        
        /*...  e-o-o */
        swep3df(ia,ja,aa,u,f,n0e,n0o,n0o,mark,nx,ny,nz);
        
        /*...  e-o-e */
        swep3df(ia,ja,aa,u,f,n0e,n0o,n0e,mark,nx,ny,nz);
        
        /*...  e-e-o */
        swep3df(ia,ja,aa,u,f,n0e,n0e,n0o,mark,nx,ny,nz);
        
        /*...  e-e-e */
        swep3df(ia,ja,aa,u,f,n0e,n0e,n0e,mark,nx,ny,nz);
        
    }
}




/**
 * \fn void fasp_smoother_dcsr_gs_rb3d (dvector *u, dCSRmat *A, dvector *b, int L,int order,int nx,int ny,int nz)
 * \brief Colores Gauss-Seidel smoother  for Au=b
 *
 * \param u     initial guess and the new approximation to the solution obtained after L GS steps
 * \param A    pointer to stiffness matrix
 * \param b    pointer to right hand side
 * \param L     number of iterations
 * \param order ordering: -1: Forward; 1: Backward
 * \param nx     number vertex of X diricter
 * \param ny     number vertex of Y diricter
 * \param nz     number vertex of Z diricter
 * \return      void
 *
 * \author FENG Chunsheng
 * \date 2011/11/23 
 */
void fasp_smoother_dcsr_gs_rb3d (dvector *u, 
                                 dCSRmat *A, 
                                 dvector *b, 
                                 int L, 
                                 int order,
                                 int *mark,
                                 int maximap,
                                 int nx,
                                 int ny,
                                 int nz )
{
	int *ia=A->IA,*ja=A->JA;
	double *aa=A->val, *bval=b->val, *uval=u->val;
    
	int i,ii,j,k,begin_row,end_row;
	int size = b->row;
    double t,d=0.0;
    // L =10;
	if (order == 1) // forward
	{
		while (L--) 
		{ 
		    rb0f3d( ia, ja, aa, uval, bval, mark, nx,  ny,  nz, 1);	
#if 1
			for (ii =0;ii <10;ii++)
                for (i = maximap; i < size; i ++) 
                {
					t = bval[i];
					begin_row = ia[i], end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aa[k]*uval[j]; 
						else d = aa[k];					
					} // end for k
					if (ABS(d) > SMALLREAL) uval[i] = t/d;
                } // end for i
#endif 
		} // end while		
	}
	else    // backward
	{
        while (L--) 
		{
#if 1		
			for (ii =0;ii <10;ii++)
                for (i = size-1; i >= maximap; i --) 
                {
					t = bval[i];
					begin_row = ia[i],end_row = ia[i+1];
					for (k = begin_row; k < end_row; k ++) 
					{
						j = ja[k];
						if (i!=j) t -= aa[k]*uval[j]; 
						else d = aa[k];					
					} // end for k
					if (ABS(d) > SMALLREAL) uval[i] = t/d;
                } // end for i
#endif 
	        rb0b3d( ia, ja, aa, uval, bval, mark, nx,  ny,  nz, 1);	
		} // end while	
	}
	return;
}


/**
 * \fn void fasp_smoother_gs_cf_omp( dvector *u, dCSRmat *A, dvector *b, int L, int *mark, int order, int nthreads, int openmp_holds )
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \param u     initial guess and the new approximation to the solution obtained after L GS steps
 * \param A    pointer to stiffness matrix
 * \param b    pointer to right hand side
 * \param L     number of iterations
 * \param mark C/F marker array
 * \param order C/F ordering: -1: F-first; 1: C-first
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return      void
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
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
	int myid,mybegin,myend;
    
	
	if (order == -1) // F-point first
	{	
		while (L--) 
		{
			if (size > openmp_holds) {
                
#pragma omp parallel for private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++)
                {
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
#pragma omp parallel for private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++)
                {
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
#pragma omp parallel for private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++)
                {
                    FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] == 1)
						{
							t = bval[i];
							begin_row = ia[i],end_row = ia[i+1];
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
						begin_row = ia[i],end_row = ia[i+1];
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
#pragma omp parallel for private(myid, mybegin, myend, i,t,begin_row,end_row,k,j,d)
                for (myid = 0; myid < nthreads; myid ++)
                {
                    FASP_GET_START_END(myid, nthreads, size, mybegin, myend);
					for (i=mybegin; i<myend; i++)
					{
						if (mark[i] != 1)
						{
							t = bval[i];
							begin_row = ia[i],end_row = ia[i+1];
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
						begin_row = ia[i],end_row = ia[i+1];
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
#endif
	return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
