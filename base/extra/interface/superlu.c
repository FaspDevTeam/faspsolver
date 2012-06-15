/*
 *  superlu.c
 *
 *------------------------------------------------------
 *  Created by Xiaozhe Hu on 11/05/09.
 *------------------------------------------------------
 *
 */

/*! \file superlu.c
 *  \brief Direct solver on the coarest level (call subroutines from SuperLU)
 *
 *  How to use SuperLU as a coarset level solver:
 *
 *	  - Download SuperLU and Install it on your own computer
 *	  - Modify msc/csrc/multigrid.c 
 *	    (On coarset level, comment pcg() and uncomment superlu())
 *	  - Modify msc/Makefile
 *	    (Define SUPERLULIB and BLASLIB according to where the libraries of SuperLU and Blas are on your computer)
 *	  - Then use "make" to compile 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if With_SuperLU
#include "slu_ddefs.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int superlu(dCSRmat *ptrA, dvector *b, dvector *u, const int print_level)
 * \brief Solve Au=b by SuperLU 
 *
 * \param *ptrA   pointer to stiffness matrix of levelNum levels
 * \param *b      pointer to the dvector of right hand side term
 * \param *u      pointer to the dvector of dofs
 *
 * \author Xiaozhe Hu
 * \data 11/05/09
 */
int superlu(dCSRmat *ptrA, dvector *b, dvector *u, const int print_level)
{

#if With_SuperLU
	
	SuperMatrix A, L, U, B;
	
	int *perm_r; /* row permutations from partial pivoting */
	int *perm_c; /* column permutation vector */
	int  nrhs=1, info, m=ptrA->row, n=ptrA->col, nnz=ptrA->nnz;
	
	if (print_level>0) printf("superlu: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	clock_t LU_start=clock();
	
	dCSRmat tempA=fasp_dcsr_create(m,n,nnz); copy_dCSRmat(ptrA,&tempA);
	dvector tempb=fasp_dvec_create(n); copy_dvector(b, &tempb);
	
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
	
	if (print_level>0) {
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

	printf("Error: SuperLU is not available!\n");
	exit(ERROR_SOLVER_TYPE);

#endif

}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
