/**
 *  interface_mumps.c		
 *  MUMPS Interface 
 *
 *------------------------------------------------------
 *
 *		Created by Chunsheng Feng on 02/27/2013.
 *
 *------------------------------------------------------
 *
 */

/*! \file interface_mumps.c
 *  \brief Call MUMPS direct solver. 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if WITH_DISOLVE
#include "dmumps_c.h"
#endif

/**
 * \fn int fasp_mumps_jobs(int *nrow, int *nnz, int *IA,int *JA, double *AA,double *b, double *x,int *job)
 * \brief Solve Ax=b by MUMPS 
 *
 * \param *nrow   pointer to the number of rows of A
 * \param *nnz    pointer to the number of nonzero entries of A
 * \param *IA     pointer to the integer array of row pointers, the size is nrow+1
 * \param *JA     pointer to the integer array of column indexes, the size is nnz
 * \param *AA     pointer to the non zero entris of A, the size is nnz
 * \param *b      pointer to the double array of right hand side term
 * \param *x      pointer to the double array of dofs
 * \param *job    1: Setup, 2: Sovle, 3 Destory
 *
 * \author Chunsheng Feng
 * \data 02/27/2013
 */

int fasp_mumps_jobs(int *nrow, int *nnz, int *IA,int *JA, double *AA,double *b, double *x,int *job)
{

#if WITH_DISOLVE
  static  DMUMPS_STRUC_C id;
  int i,j;
  int n =  *nrow;
  int nz = *nnz;
  int *irn = id.irn;
  int *jcn = id.jcn;
  double *a = id.a;
  double *rhs = id.rhs;
  
  switch( *job ) {
  case 1:
	{
  /* Define A and rhs */
  irn = (int *)malloc( sizeof(int)*nz );
  jcn = (int *)malloc( sizeof(int)*nz );
  a   = (double *)malloc( sizeof(double)*nz );
  rhs = (double *)malloc( sizeof(double)*n );

  
  int begin_row,end_row;

  if (IA[0]== 0) {
    for (i=0; i<n; i++) {
        begin_row = IA[i]; end_row = IA[i+1];
     for (j=begin_row; j< end_row; j++)     { 
        irn[j] = i + 1;
        jcn[j] = JA[j]+1;
        a[j] =   AA[j];
     }
    }
  } else if (IA[0] ==1 ){
    for (i=0; i<n; i++) {
      begin_row = IA[i]-1; end_row = IA[i+1]-1;
     for (j=begin_row; j< end_row; j++)     { 
        irn[j] = i + 1;
        jcn[j] = JA[j];
        a[j] =   AA[j];
     }
    }
  } else 
  printf("Matrix is wrong format IA[0] = %d \n", IA[0] );


  /* Initialize a MUMPS instance. */
  id.job = -1; id.par=1; id.sym=0;id.comm_fortran=0;
  dmumps_c(&id);
  /* Define the problem on the host */
    id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
    id.a = a; id.rhs = rhs;
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  /* No outputs */
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;	
    id.job=4;
    dmumps_c(&id);  
	}
	break;


  case 2:
	{
/* Call the MUMPS package. */
    for(i=0; i<n; i++) rhs[i] = b[i];
  
    id.job=3;  dmumps_c(&id);
  
    for(i=0; i<n; i++) x[i] = id.rhs[i];	
	}
	break;

  case 3:
	{
     id.job = -2; 
     dmumps_c(&id); /* Terminate instance */
     free(irn);
     free(jcn);
     free(a);
     free(rhs);	
	}
	break;
	
  default:
    printf("Warning job should be 1,2,3 !!!!!!!!!!!!!!!!!!\n");
    break;  
	
  }
#endif
  return 0;
}


/**
 * \fn int fasp_mumps_jobs(int *nrow, int *nnz, int *IA,int *JA, double *AA,double *b, double *x,int *job)
 * \brief Solve Ax=b by MUMPS 
 *
 * \param *nrow   pointer to the number of rows of A
 * \param *nnz    pointer to the number of nonzero entries of A
 * \param *IA     pointer to the integer array of row pointers, the size is nrow+1
 * \param *JA     pointer to the integer array of column indexes, the size is nnz
 * \param *AA     pointer to the non zero entris of A, the size is nnz
 * \param *b      pointer to the double array of right hand side term
 * \param *x      pointer to the double array of dofs
 *
 * \author Chunsheng Feng
 * \data 02/27/2013
 */

int fasp_mumps_direct(int *nrow, int *nnz, int *IA,int *JA, double *AA,double *b, double *x)
{
#if WITH_DISOLVE
  static  DMUMPS_STRUC_C id;
  int i,j;
  int n =  *nrow;
  int nz = *nnz;
  int *irn;
  int *jcn;
  double *a;
  double *rhs;
  
  /* Define A and rhs */
  irn = (int *)malloc( sizeof(int)*nz );
  jcn = (int *)malloc( sizeof(int)*nz );
  a   = (double *)malloc( sizeof(double)*nz );
  rhs = (double *)malloc( sizeof(double)*n );
  
  int begin_row,end_row;
  if (IA[0]== 0) {
    for (i=0; i<n; i++) {
        begin_row = IA[i]; end_row = IA[i+1];
     for (j=begin_row; j< end_row; j++)     { 
        irn[j] = i + 1;
        jcn[j] = JA[j]+1;
        a[j] =   AA[j];
     }
    }
  } else if (IA[0] ==1 ){
    for (i=0; i<n; i++) {
      begin_row = IA[i]-1; end_row = IA[i+1]-1;
     for (j=begin_row; j< end_row; j++)     { 
        irn[j] = i + 1;
        jcn[j] = JA[j];
        a[j] =   AA[j];
     }
    }
  } else 
  printf("Matrix is wrong format IA[0] = %d \n", IA[0] );

  /* Initialize a MUMPS instance. */
  id.job=-1; id.par=1; id.sym=0;id.comm_fortran=0;
  dmumps_c(&id);
  /* Define the problem on the host */
    id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
    id.a = a; id.rhs = rhs;

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
/* No outputs */
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;

/* Call the MUMPS package. */
  for(i=0; i<n; i++) rhs[i] = b[i];
 
  id.job=6;
  dmumps_c(&id);
  
  for(i=0; i<n; i++) x[i] = id.rhs[i];

  id.job=-2; 
  dmumps_c(&id); /* Terminate instance */

  free(irn);
  free(jcn);
  free(a);
  free(rhs);
#endif

  return 0;
}
