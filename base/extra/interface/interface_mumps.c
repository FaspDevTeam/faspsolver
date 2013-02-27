/*! \file interface_mumps.c
 *  \brief Call MUMPS direct solver.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if WITH_MUMPS
#include "dmumps_c.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_solver_mumps (int *nrow, int *nnz, int *IA, int *JA,
 *                            double *AA, double *b, double *x)
 *
 * \brief Solve Ax=b by MUMPS directly
 *
 * \param nrow   pointer to the number of rows of A
 * \param nnz    pointer to the number of nonzero entries of A
 * \param IA     pointer to the integer array of row pointers, the size is nrow+1
 * \param JA     pointer to the integer array of column indexes, the size is nnz
 * \param AA     pointer to the non zero entris of A, the size is nnz
 * \param b      pointer to the double array of right hand side term
 * \param x      pointer to the double array of dofs
 *
 * \author Chunsheng Feng
 * \data   02/27/2013
 *
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 *
 */
int fasp_solver_mumps (int *nrow,
                       int *nnz,
                       int *IA,
                       int *JA,
                       double *AA,
                       double *b,
                       double *x)
{
    
#if WITH_MUMPS
    
    static DMUMPS_STRUC_C id;
    int i,j;
    int n =  *nrow;
    int nz = *nnz;
    int *irn;
    int *jcn;
    double *a;
    double *rhs;
    
    int begin_row, end_row;
    
    // First check the matrix format
    if ( IA[0] != 0 && IA[0] != 1 ) {
        printf("### ERROR: Matrix format is wrong -- IA[0] = %d\n", IA[0]);
        return ERROR_SOLVER_EXIT;
    }
    
    /* Define A and rhs */
    irn = (int *)malloc( sizeof(int)*nz );
    jcn = (int *)malloc( sizeof(int)*nz );
    a   = (double *)malloc( sizeof(double)*nz );
    rhs = (double *)malloc( sizeof(double)*n );
    
    if ( IA[0] == 0 ) { // C-convention
        for (i=0; i<n; i++) {
            begin_row = IA[i]; end_row = IA[i+1];
            for (j=begin_row; j< end_row; j++)     {
                irn[j] = i + 1;
                jcn[j] = JA[j]+1;
                a[j] =   AA[j];
            }
        }
    }
    else { // For-convention
        for (i=0; i<n; i++) {
            begin_row = IA[i]-1; end_row = IA[i+1]-1;
            for (j=begin_row; j< end_row; j++)     {
                irn[j] = i + 1;
                jcn[j] = JA[j];
                a[j] =   AA[j];
            }
        }
    }
    
    /* Initialize a MUMPS instance. */
    id.step=-1; id.par=1; id.sym=0;id.comm_fortran=0;
    dmumps_c(&id);
    /* Define the problem on the host */
    id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
    id.a = a; id.rhs = rhs;
    
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    /* No outputs */
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
    
    /* Call the MUMPS package. */
    for(i=0; i<n; i++) rhs[i] = b[i];
    
    id.step=6;
    dmumps_c(&id);
    
    for(i=0; i<n; i++) x[i] = id.rhs[i];
    
    id.step=-2;
    dmumps_c(&id); /* Terminate instance */
    
    free(irn);
    free(jcn);
    free(a);
    free(rhs);
    
    return SUCCESS;
    
#else
    
	printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

/**
 * \fn int fasp_solver_mumps_steps (int *nrow, int *nnz, int *IA, int *JA,
 *                                  double *AA, double *b, double *x, int *step)
 *
 * \brief Solve Ax=b by MUMPS in three steps
 *
 * \param nrow   pointer to the number of rows of A
 * \param nnz    pointer to the number of nonzero entries of A
 * \param IA     pointer to the integer array of row pointers, the size is nrow+1
 * \param JA     pointer to the integer array of column indexes, the size is nnz
 * \param AA     pointer to the non zero entris of A, the size is nnz
 * \param b      pointer to the double array of right hand side term
 * \param x      pointer to the double array of dofs
 * \param step   1: Setup, 2: Sovle, 3 Destory
 *
 * \author Chunsheng Feng
 * \data   02/27/2013
 * 
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 */
int fasp_solver_mumps_steps (int *nrow,
                             int *nnz,
                             int *IA,
                             int *JA,
                             double *AA,
                             double *b,
                             double *x,
                             int *step)
{
    
#if WITH_MUMPS
    static  DMUMPS_STRUC_C id;
    int i,j;
    int n =  *nrow;
    int nz = *nnz;
    int *irn = id.irn;
    int *jcn = id.jcn;
    double *a = id.a;
    double *rhs = id.rhs;
    
    switch( *step ) {
            
        case 1:
        {
            
            int begin_row, end_row;
            
            // First check the matrix format
            if ( IA[0] != 0 && IA[0] != 1 ) {
                printf("### ERROR: Matrix format is wrong -- IA[0] = %d\n", IA[0]);
                return ERROR_SOLVER_EXIT;
            }
            
            // Define A and rhs
            irn = (int *)malloc( sizeof(int)*nz );
            jcn = (int *)malloc( sizeof(int)*nz );
            a   = (double *)malloc( sizeof(double)*nz );
            rhs = (double *)malloc( sizeof(double)*n );
            
            if ( IA[0] == 0 ) { // C-convention
                for (i=0; i<n; i++) {
                    begin_row = IA[i]; end_row = IA[i+1];
                    for (j=begin_row; j< end_row; j++)     {
                        irn[j] = i + 1;
                        jcn[j] = JA[j]+1;
                        a[j] =   AA[j];
                    }
                }
            }
            else { // For-convention
                for (i=0; i<n; i++) {
                    begin_row = IA[i]-1; end_row = IA[i+1]-1;
                    for (j=begin_row; j< end_row; j++)     {
                        irn[j] = i + 1;
                        jcn[j] = JA[j];
                        a[j] =   AA[j];
                    }
                }
            }
            
            /* Initialize a MUMPS instance. */
            id.step = -1; id.par=1; id.sym=0;id.comm_fortran=0;
            dmumps_c(&id);
            /* Define the problem on the host */
            id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
            id.a = a; id.rhs = rhs;
            
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
            /* No outputs */
            id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
            id.step=4;
            dmumps_c(&id);
        }
            break;
            
        case 2:
        {
            /* Call the MUMPS package. */
            for(i=0; i<n; i++) rhs[i] = b[i];
            
            id.step=3;  dmumps_c(&id);
            
            for(i=0; i<n; i++) x[i] = id.rhs[i];
        }
            break;
            
        case 3:
        {
            id.step = -2;
            dmumps_c(&id); /* Terminate instance */
            free(irn);
            free(jcn);
            free(a);
            free(rhs);
        }
            break;
            
        default:
            printf("### ERROR: Parameter step should be 1, 2, or 3!\n");
            return ERROR_SOLVER_EXIT;
            
    }
    
    return SUCCESS;
    
#else
    
	printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
