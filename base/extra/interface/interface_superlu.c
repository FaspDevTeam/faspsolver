/*! \file interface_superlu.c
 *  \brief Interface to SuperLU direct solvers
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if WITH_SuperLU
#include "slu_ddefs.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_solver_superlu (dCSRmat *ptrA, dvector *b, dvector *u,
 *                              const int print_level)
 *
 * \brief Solve Au=b by SuperLU
 *
 * \param ptrA         Pointer to a dCSRmat matrix
 * \param b            Pointer to the dvector of right-hand side term
 * \param u            Pointer to the dvector of solution
 * \param print_level  Output level
 *
 * \author Xiaozhe Hu
 * \data   11/05/09
 *
 * Modified by Chensong Zhang on 11/01/2012 for new FASP function names.
 * Modified by Chensong Zhang on 02/27/2013 for new FASP function names.
 *
 */
int fasp_solver_superlu (dCSRmat *ptrA,
                         dvector *b,
                         dvector *u,
                         const int print_level)
{
    
#if WITH_SuperLU
	
	SuperMatrix A, L, U, B;
	
	int *perm_r; /* row permutations from partial pivoting */
	int *perm_c; /* column permutation vector */
	int  nrhs=1, info, m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;
	
	if (print_level>0) printf("superlu: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	clock_t LU_start=clock();
	
	dCSRmat tempA=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(ptrA,&tempA);
    
	dvector tempb=fasp_dvec_create(n);       fasp_dvec_cp(b, &tempb);
	
	/* Create matrix A in the format expected by SuperLU. */
	dCreate_CompCol_Matrix(&A, m, n, nnz, tempA.val, tempA.JA, tempA.IA, SLU_NR, SLU_D, SLU_GE);
	
	/* Create right-hand side B. */
	dCreate_Dense_Matrix(&B, m, nrhs, tempb.val, m, SLU_DN, SLU_D, SLU_GE);
	
	if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
	
	/* Set the default input options. */
	superlu_options_t options;
	set_default_options(&options);
	//options.PrintStat = NO;
	options.ColPerm = NATURAL;
	
	/* Initialize the statistics variables. */
	SuperLUStat_t stat;
	StatInit(&stat);
	/* SuperLU */
	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
	
	DNformat *BB = (DNformat *) B.Store;
	u->val = (double *) BB->nzval;
	u->row = n;
	
	if ( print_level > PRINT_MIN ) {
		clock_t LU_end=clock();
		double LUduration = (double)(LU_end - LU_start)/(double)(CLOCKS_PER_SEC);
		printf("SuperLU totally costs %f seconds.\n", LUduration);
	}
	
	/* De-allocate storage */
	SUPERLU_FREE(perm_r);
	SUPERLU_FREE(perm_c);
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	StatFree(&stat);
	
	return SUCCESS;
    
#else
    
	printf("### ERROR: SuperLU is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
