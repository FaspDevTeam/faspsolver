/*! \file interface_umfpack.c
 *  \brief Interface to UMFPACK direct solvers
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if WITH_UMFPACK
#include "umfpack.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_solver_umfpack (dCSRmat *ptrA, dvector *b, dvector *u, 
 *                              const int print_level)
 *
 * \brief Solve Au=b by UMFpack
 *
 * \param ptrA         Pointer to a dCSRmat matrix
 * \param b            Pointer to the dvector of right-hand side term
 * \param u            Pointer to the dvector of solution
 * \param print_level  Output level
 *
 * \author Chensong Zhang
 * \date   05/20/2010
 *
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 */
int fasp_solver_umfpack (dCSRmat *ptrA,
                         dvector *b,
                         dvector *u,
                         const int print_level)
{
    
#if WITH_UMFPACK
    
	const int n = ptrA->col;
	const int m = ptrA->row;
	const int nnz = ptrA->nnz;
	
	int *Ap = ptrA->IA;
	int *Ai = ptrA->JA;
	double *Ax = ptrA->val;
	void *Symbolic, *Numeric;
	int status=SUCCESS;
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_umfpack ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
	clock_t start_time = clock();
    
	status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
	status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
	umfpack_di_free_symbolic (&Symbolic);
	status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, u->val, b->val, Numeric, NULL, NULL);
	umfpack_di_free_numeric (&Numeric);
	
	if ( print_level > PRINT_MIN ) {
		clock_t end_time = clock();
		double solve_duration = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
		printf("UMFPACK costs %f seconds.\n", solve_duration);
	}
	
#if DEBUG_MODE
	printf("### DEBUG: fasp_solver_umfpack ...... [Finish]\n");
#endif
	
	return status;
	
#else
    
	printf("### ERROR: UMFPACK is not available!\n");
    return ERROR_SOLVER_EXIT;
	
#endif
    
}

#if WITH_UMFPACK

/**
 * \fn int umfpack_factorize (dCSRmat *ptrA, dvector *b, dvector *u,
 *                            const int print_level)
 * \brief factorize A by UMFpack
 *
 * \param ptrA      pointer to stiffness matrix of levelNum levels
 * \param Numeric   pointer to the numerical factorization
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 */
int umfpack_factorize (dCSRmat *ptrA,
                       void *Numeric,
                       const int print_level)
{
	const int n = ptrA->col;
	const int m = ptrA->row;
	const int nnz = ptrA->nnz;
	
	int *Ap = ptrA->IA;
	int *Ai = ptrA->JA;
	double *Ax = ptrA->val;
	void *Symbolic;
	int status=SUCCESS;
	
#if DEBUG_MODE
	printf("### DEBUG: umfpack_factorize ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
	clock_t start_time = clock();
    
	status = umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
	status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
	umfpack_di_free_symbolic (&Symbolic);
	
	if ( print_level > PRINT_MIN ) {
		clock_t end_time = clock();
		double fac_duration = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
		printf("UMFPACK factorize costs %f seconds.\n", fac_duration);
	}
	
#if DEBUG_MODE
	printf("### DEBUG: umfpack_factorize ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn int umfpack_solve (dCSRmat *ptrA, dvector *b, dvector *u,
 *                        void *Numeric, const int print_level)
 * \brief Solve Au=b by UMFpack, numerical factorization is given
 *
 * \param ptrA      pointer to stiffness matrix of levelNum levels
 * \param b         pointer to the dvector of right hand side term
 * \param u         pointer to the dvector of dofs
 * \param Numeric   pointer to the numerical factorization
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 */

int umfpack_solve (dCSRmat *ptrA,
                   dvector *b,
                   dvector *u,
                   void *Numeric,
                   const int print_level)
{
	const int n = ptrA->col;
	const int m = ptrA->row;
	const int nnz = ptrA->nnz;
	
	int *Ap = ptrA->IA;
	int *Ai = ptrA->JA;
	double *Ax = ptrA->val;
	int status=SUCCESS;
	
#if DEBUG_MODE
	printf("### DEBUG: umfpack_solve ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
	clock_t start_time = clock();
    
	status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, u->val, b->val, Numeric, NULL, NULL);
	
	if (print_level>0) {
		clock_t end_time = clock();
		double solve_duration = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
		printf("UMFPACK costs %f seconds.\n", solve_duration);
	}
	
#if DEBUG_MODE
	printf("### DEBUG: umfpack_solve ...... [Finish]\n");
#endif
	
	return status;
}

/**
 * \fn int umfpack_free_numeric (void *Numeric)
 * \brief Solve Au=b by UMFpack
 *
 * \param ptrA     pointer to stiffness matrix of levelNum levels
 * \param b        pointer to the dvector of right hand side term
 * \param u        pointer to the dvector of dofs
 *
 * \author Xiaozhe Hu
 * \date   10/20/2010
 */
int umfpack_free_numeric (void *Numeric)
{
	int status=SUCCESS;
	
#if DEBUG_MODE
	printf("### DEBUG: umfpack_free_numeric ...... [Start]\n");
#endif
 	
	umfpack_di_free_numeric (&Numeric);
	
#if DEBUG_MODE
	printf("### DEBUG: umfpack_free_numeric ...... [Finish]\n");
#endif
	
	return status;
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
