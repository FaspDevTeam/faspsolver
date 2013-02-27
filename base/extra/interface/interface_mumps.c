/*! \file interface_mumps.c
 *  \brief Call MUMPS direct solver.
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if WITH_MUMPS
#include "dmumps_c.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_solver_mumps (dCSRmat *ptrA, dvector *b, dvector *u, 
 *                            const int print_level)
 *
 * \brief Solve Ax=b by MUMPS directly
 *
 * \param ptrA         Pointer to a dCSRmat matrix
 * \param b            Pointer to the dvector of right-hand side term
 * \param u            Pointer to the dvector of solution
 * \param print_level  Output level
 *
 * \author Chunsheng Feng
 * \data   02/27/2013
 *
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 *
 */
int fasp_solver_mumps ( dCSRmat *ptrA,
                        dvector *b,
                        dvector *u,
                        const int print_level)
{
    
#if WITH_MUMPS
    
    DMUMPS_STRUC_C id;

    const  int n =  ptrA->row;
    const  int nz = ptrA->nnz;
    int *IA = ptrA->IA;
    int *JA = ptrA->JA;
    double *AA =  ptrA->val;
    double *b1 = b->val;
    double *x  = u->val;
	
    int *irn;
    int *jcn;
    double *a;
    double *rhs;
    int i,j;
    int begin_row, end_row;

#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_mumps ...... [Start]\n");
	printf("### DEBUG: nr=%d,  nnz=%d\n",  n, nz);
#endif
    
    // First check the matrix format
    if ( IA[0] != 0 && IA[0] != 1 ) {
        printf("### ERROR: Matrix format is wrong -- IA[0] = %d\n", IA[0]);
        return ERROR_SOLVER_EXIT;
    }

	clock_t start_time = clock();
    
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
    id.job=-1; id.par=1; id.sym=0; id.comm_fortran=0;
    dmumps_c(&id);
    /* Define the problem on the host */
    id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
    id.a = a; id.rhs = rhs;
    
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    /* No outputs */
    id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
    
    /* Call the MUMPS package. */
    for(i=0; i<n; i++) rhs[i] = b1[i];
    
    id.job=6; dmumps_c(&id);
    
    for(i=0; i<n; i++) x[i] = id.rhs[i];
    
    id.job=-2;
    dmumps_c(&id); /* Terminate instance */
    
    free(irn);
    free(jcn);
    free(a);
    free(rhs);
	
	if ( print_level > PRINT_MIN ) {
		clock_t end_time = clock();
		double solve_duration = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
		printf("MUMPS costs %f seconds.\n", solve_duration);
	}  

#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_mumps ...... [Finish]\n");
#endif	
    return SUCCESS;
#else
    
	printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

/**
 * \fn int fasp_solver_mumps_steps (dCSRmat *ptrA, dvector *b, dvector *u, 
 *                                  const int job)
 *
 * \brief Solve Ax=b by MUMPS in three steps
 *
 * \param ptrA   Pointer to a dCSRmat matrix
 * \param b      Pointer to the dvector of right-hand side term
 * \param u      Pointer to the dvector of solution
 * \param job    1: Setup, 2: Sovle, 3 Destory
 *
 * \author Chunsheng Feng
 * \data   02/27/2013
 * 
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 */
int fasp_solver_mumps_steps ( dCSRmat *ptrA,
                              dvector *b,
                              dvector *u,
                              const int job)
{
    
#if WITH_MUMPS
    static  DMUMPS_STRUC_C id;
    int i,j;
	
    const  int n =  ptrA->row;
    const  int nz = ptrA->nnz;
    int *IA = ptrA->IA;
    int *JA = ptrA->JA;
    double *AA =  ptrA->val;
    double *b1 = b->val;
    double *x  = u->val;

    int *irn = id.irn;
    int *jcn = id.jcn;
    double *a = id.a;
    double *rhs = id.rhs;
    
    switch ( job ) {
            
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
            id.job = -1; id.par=1; id.sym=0; id.comm_fortran=0;
            dmumps_c(&id);
            /* Define the problem on the host */
            id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
            id.a = a; id.rhs = rhs;
            
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
            /* No outputs */
            id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;

            id.job=4; dmumps_c(&id);
        }
            break;
            
        case 2:
        {
            /* Call the MUMPS package. */
            for(i=0; i<id.n; i++) rhs[i] = b1[i];
            
            id.job=3; dmumps_c(&id);
            
            for(i=0; i<id.n; i++) x[i] = id.rhs[i];
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
            printf("### ERROR: Parameter job should be 1, 2, or 3!\n");
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
