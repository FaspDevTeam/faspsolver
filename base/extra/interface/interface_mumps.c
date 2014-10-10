/*! \file interface_mumps.c
 *  \brief Interface to MUMPS direct solvers
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
 * \date   02/27/2013
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
    return FASP_SUCCESS;
#else
    
	printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

/**
 * \fn int fasp_solver_mumps_steps (dCSRmat *ptrA, dvector *b, dvector *u, 
 *                                  Mumps_param *mumps)
 *
 * \brief Solve Ax=b by MUMPS in three steps
 *
 * \param ptrA   Pointer to a dCSRmat matrix
 * \param b      Pointer to the dvector of right-hand side term
 * \param u      Pointer to the dvector of solution
 * \param mumps  Pointer to MUMPS parameters
 *
 * \author Chunsheng Feng
 * \date   02/27/2013
 * 
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 * Modified by Zheng Li on 10/10/2014 to adjust input parameters.
 * 
 */
int fasp_solver_mumps_steps ( dCSRmat *ptrA,
                              dvector *b,
                              dvector *u,
                              Mumps_data *mumps) 
{
    
#if WITH_MUMPS

    DMUMPS_STRUC_C id;

	int job = mumps->job;

	static int job_stat = 0;
    int i,j;

	int *irn;
    int *jcn;
    double*a;
    int *rhs;
    
#if DEBUG_MODE
	printf("### DEBUG: %s job_stat = %d\n", __FUNCTION__, job_stat);
#endif

    switch ( job ) {
            
        case 1:
        {
#if DEBUG_MODE
           printf("### DEBUG: %s Step 1 ...... [Start]\n", __FUNCTION__);         
#endif
            int begin_row, end_row;
            const  int n =  ptrA->row;
            const  int nz = ptrA->nnz;
            int *IA = ptrA->IA;
            int *JA = ptrA->JA;
            double *AA =  ptrA->val;

            irn = id.irn;
            jcn = id.jcn;
            a   = id.a;
            rhs = id.rhs;

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
            job_stat = 1;

			mumps->id = id;

#if DEBUG_MODE
			printf("### DEBUG: %s, Step 1 ...... [Finish]\n", __FUNCTION__);   
#endif
        }
            break;
            
        case 2:
        {
#if DEBUG_MODE
			printf("### DEBUG: %s, Step 2 ...... [Start]\n", __FUNCTION__);   
#endif
			id = mumps->id;

            if ( job_stat != 1 )
                printf("### ERROR: fasp_solver_mumps_steps has not finish Setup...... [2]\n");

            /* Call the MUMPS package. */
            for(i=0; i<id.n; i++) id.rhs[i] = b->val[i];
                    
            id.job=3; dmumps_c(&id);
                    
            for(i=0; i<id.n; i++) u->val[i] = id.rhs[i];
#if DEBUG_MODE
			printf("### DEBUG: %s, Step 2 ...... [Finish]\n", __FUNCTION__);   
#endif
        }
            break;
            
        case 3:
        {
			id = mumps->id;

            if ( job_stat !=1 )
                printf("### ERROR: %s has not been setted up!\n", __FUNCTION__);
            
            free(id.irn);
            free(id.jcn);
            free(id.a);
            free(id.rhs);
            //id.job = -2;
            //dmumps_c(&id); /* Terminate instance */
        }
            break;
            
        default:
            printf("### ERROR: Parameter job should be 1, 2, or 3!\n");
            return ERROR_SOLVER_EXIT;
            
    }
    
    return FASP_SUCCESS;
    
#else
    
    printf("### ERROR: MUMPS is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

#if WITH_MUMPS
/**
 ** \fn DMUMPS_STRUC_C fasp_mumps_factorize (dCSRmat *ptrA, dvector *b, dvector *u,
 **                                   const INT print_level)
 ** \brief factorize A by MUMPS
 **
 ** \param ptrA      pointer to stiffness matrix of levelNum levels
 ** \param b         pointer to the dvector of right hand side term
 ** \param u         pointer to the dvector of dofs
 **
 ** \author Zheng Li
 ** \date   10/09/2014
 **/
Mumps_data fasp_mumps_factorize (dCSRmat *ptrA,
                                 dvector *b,
                                 dvector *u,
								 const INT print_level)
{
	Mumps_data mumps;
    DMUMPS_STRUC_C id;

    int i,j;	
    const  int m =  ptrA->row;
    const  int n =  ptrA->col;
    const  int nz = ptrA->nnz;
    int *IA = ptrA->IA;
    int *JA = ptrA->JA;
    double *AA =  ptrA->val;

    int *irn = id.irn;
    int *jcn = id.jcn;
    double *a = id.a;
    double *rhs = id.rhs;
    
    int begin_row, end_row;

#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nz);
#endif

    clock_t start_time = clock();

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

   if ( print_level > PRINT_MIN ) {
       clock_t end_time = clock();
       double fac_duration = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
       printf("UMFPACK factorize costs %f seconds.\n", fac_duration);
   }

#if DEBUG_MODE
   printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);   
#endif

   mumps.id = id;

   return mumps;
}

#endif 



#if WITH_MUMPS
/**
 ** \fn void fasp_mumps_solve( dCSRmat *ptrA, dvector *b, dvector *u, DMUMPS_STRUC_C id)
 **                            const INT print_level)
 ** \brief solve A by MUMPS
 **
 ** \param ptrA      pointer to stiffness matrix of levelNum levels
 ** \param b         pointer to the dvector of right hand side term
 ** \param u         pointer to the dvector of dofs
 ** \param id        pointer to the numerical factorization
 **
 ** \author Zheng Li
 ** \date   10/09/2014
 **/
void fasp_mumps_solve ( dCSRmat *ptrA,
                        dvector *b,
                        dvector *u,
                        Mumps_data mumps,
						const INT print_level)
{
    int i,j;

	DMUMPS_STRUC_C id = mumps.id;
  
    const  int m =  ptrA->row;
    const  int n =  ptrA->row;
    const  int nz = ptrA->nnz;
    int *IA = ptrA->IA;
    int *JA = ptrA->JA;
    double *AA =  ptrA->val;

    int *irn = id.irn;
    int *jcn = id.jcn;
    double *a = id.a;
    double *rhs = id.rhs;
    
#if DEBUG_MODE
	printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nz);
#endif

	clock_t start_time = clock();

    double *b1 = b->val;
    double *x  = u->val;

    /* Call the MUMPS package. */
    for(i=0; i<id.n; i++) rhs[i] = b1[i];
                    
    id.job=3; dmumps_c(&id);
                    
    for(i=0; i<id.n; i++) x[i] = id.rhs[i];

    if (print_level>0) {
        clock_t end_time = clock();
        double solve_duration = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
        printf("UMFPACK costs %f seconds.\n", solve_duration);
    }

#if DEBUG_MODE
	printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);   
#endif
    
}
#endif

#if WITH_MUMPS
/**
 ** \fn void fasp_mumps_free (Mumps_data *mumps)
 **                                 
 ** \brief free memory
 ** 
 ** \param mumps Pointer to mumps data
 **
 ** \author Zheng Li
 ** \date   10/09/2014
 **/
void fasp_mumps_free ( Mumps_data *mumps )
{
    DMUMPS_STRUC_C id = mumps->id;

    free(id.irn);
    free(id.jcn);
    free(id.a);
    free(id.rhs);
}

#endif 
/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
