/*! \file InterfacePardiso.c
 *
 *  \brief Interface to Intel MKL PARDISO direct solvers
 *
 *  Reference for Intel MKL PARDISO:
 *  https://software.intel.com/en-us/node/470282
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#if WITH_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_solver_pardiso (dCSRmat *ptrA, dvector *b, dvector *u,
 *                            const SHORT prtlvl)
 *
 * \brief Solve Ax=b by PARDISO directly. Each row of A should be
 *        in ascending order w.r.t. column indices.
 *
 * \param ptrA      Pointer to a dCSRmat matrix
 * \param b         Pointer to the dvector of right-hand side term
 * \param u         Pointer to the dvector of solution
 * \param prtlvl    Output level
 *
 * \author Hongxuan Zhang
 * \date   11/28/2015
 */
INT fasp_solver_pardiso (dCSRmat * ptrA,
                         dvector *b,
                         dvector *u,
                         const SHORT prtlvl)
{
#if WITH_PARDISO
    INT status = FASP_SUCCESS;
    
    MKL_INT n = ptrA->col;
    MKL_INT m = ptrA->row;
    MKL_INT nnz = ptrA->nnz;
    MKL_INT *ia = ptrA->IA;
    MKL_INT *ja = ptrA->JA;
    REAL *a = ptrA->val;
    
    MKL_INT mtype = 11;    /* Real unsymmetric matrix */
    MKL_INT nrhs = 1;      /* Number of right hand sides. */
    MKL_INT idum;          /* Integer dummy. */
    MKL_INT iparm[64];     /* Pardiso control parameters. */
    MKL_INT maxfct, mnum, phase, error, msglvl;    /* Auxiliary variables. */
    
    REAL * f = b->val;     /* RHS vector. */
    REAL * x = u->val;     /* Solution vector. */
    void *pt[64];          /* Internal solver memory pointer pt. */
    double ddum;           /* Double dummy */
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    PARDISOINIT(pt, &mtype, iparm); /* Initialize. */
    iparm[34] = 1;        /* Use 0-based indexing. */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 0;           /* Do not print statistical information in file */
    error = 0;            /* Initialize error flag */
    
    clock_t start_time = clock();
    
    phase = 11; /* Reordering and symbolic factorization. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("### ERROR: Symbolic factorization failed %d!\n", error);
        exit (1);
    }
    
    phase = 22; /* Numerical factorization. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("\n### ERROR: Numerical factorization failed %d!\n", error);
        exit (2);
    }
    
    phase = 33; /* Back substitution and iterative refinement. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, f, x, &error);
    
    if ( error != 0 ) {
        printf ("\n### ERROR: Solution failed %d!\n", error);
        exit (3);
    }
    
    if ( prtlvl > PRINT_MIN ) {
        clock_t end_time = clock();
        double solve_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
        printf("PARDISO costs %f seconds.\n", solve_time);
    }
    
    phase = -1; /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
    
#else
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    printf("### ERROR: PARDISO is not available!\n");
    return ERROR_SOLVER_EXIT;
    
#endif
    
}

#if WITH_PARDISO
/**
 * \fn INT fasp_pardiso_factorize (dCSRmat *ptrA, Pardiso_data *pdata,
 *                                 const SHORT prtlvl)
 * \brief factorize A by PARDISO. Each row of A should be
 *        in ascending order w.r.t. column indices.
 *
 * \param ptrA      Pointer to matrix A
 * \param pdata     Pointer to numerical factorization data
 * \param prtlvl    Output level
 *
 * \author Hongxuan Zhang
 * \date   11/28/2015
 */
INT fasp_pardiso_factorize (dCSRmat *ptrA,
                            Pardiso_data *pdata,
                            const SHORT prtlvl)
{
    INT status = FASP_SUCCESS;
    
    MKL_INT n = ptrA->col;
    MKL_INT m = ptrA->row;
    MKL_INT nnz = ptrA->nnz;
    MKL_INT *ia = ptrA->IA;
    MKL_INT *ja = ptrA->JA;
    REAL *a = ptrA->val;
    
    pdata->mtype = 11;    /* Real unsymmetric matrix */
    
    double  ddum;         /* Double dummy */
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    MKL_INT idum;         /* Integer dummy. */
    MKL_INT phase, error, msglvl;    /* Auxiliary variables. */
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    PARDISOINIT(pdata->pt, &(pdata->mtype), pdata->iparm); /* Initialize. */
    pdata->iparm[34] = 1;  /* Use 0-based indexing. */
    pdata->maxfct = 1;     /* Maximum number of numerical factorizations. */
    pdata->mnum = 1;       /* Which factorization to use. */
    msglvl = 0;            /* Do not print statistical information in file */
    error = 0;             /* Initialize error flag */
    
    clock_t start_time = clock();
    
    phase = 11; /* Reordering and symbolic factorization. */
    PARDISO (pdata->pt, &(pdata->maxfct), &(pdata->mnum), &(pdata->mtype), &phase,
             &n, a, ia, ja, &idum, &nrhs, pdata->iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ) {
        printf ("\n### ERROR: Symbolic factorization failed %d!\n", error);
        exit (1);
    }
    
    phase = 22; /* Numerical factorization. */
    PARDISO (pdata->pt, &(pdata->maxfct), &(pdata->mnum), &(pdata->mtype), &phase,
             &n, a, ia, ja, &idum, &nrhs, pdata->iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\n### ERROR: Numerical factorization failed %d!\n", error);
        exit (2);
    }
    
    if ( prtlvl > PRINT_MIN ) {
        clock_t end_time = clock();
        double solve_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
        printf("PARDISO factoization costs %f seconds.\n", solve_time);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    return status;
}

/**
 * \fn INT fasp_pardiso_solve (dCSRmat *ptrA, dvector *b, dvector *u,
 *                             Pardiso_data *pdata, const SHORT prtlvl)
 * \brief Solve Au=b by Intel MKL PARDISO, numerical factorization is given.
 *        Each row of A should be in ascending order w.r.t. column indices.
 *
 * \param ptrA      Pointer to stiffness matrix of levelNum levels
 * \param b         Pointer to the dvector of right hand side term
 * \param u         Pointer to the dvector of dofs
 * \param pdata     Pointer to the numerical factorization data
 * \param prtlvl    Output level
 *
 * \author Hongxuan Zhang
 * \date   11/28/2015
 */
INT fasp_pardiso_solve (dCSRmat *ptrA,
                        dvector *b,
                        dvector *u,
                        Pardiso_data *pdata,
                        const SHORT prtlvl)
{
    INT status = FASP_SUCCESS;
    
    MKL_INT n = ptrA->col;
    MKL_INT m = ptrA->row;
    MKL_INT nnz = ptrA->nnz;
    MKL_INT *ia = ptrA->IA;
    MKL_INT *ja = ptrA->JA;
    
    REAL *a = ptrA->val;
    REAL * f = b->val;     /* RHS vector. */
    REAL * x = u->val;     /* Solution vector. */
    MKL_INT mtype = 11;    /* Real unsymmetric matrix */
    MKL_INT nrhs = 1;      /* Number of right hand sides. */
    MKL_INT idum;          /* Integer dummy. */
    MKL_INT phase, error, msglvl;    /* Auxiliary variables. */
    
    msglvl = 0; /* Do not print statistical information in file */
    
    clock_t start_time = clock();
    
    phase = 33; /* Back substitution and iterative refinement. */
    PARDISO (pdata->pt, &(pdata->maxfct), &(pdata->mnum), &(pdata->mtype), &phase,
             &n, a, ia, ja, &idum, &nrhs, pdata->iparm, &msglvl, f, x, &error);
    
    if ( error != 0 ) {
        printf ("### ERROR: Solution failed %d!\n", error);
        exit (3);
    }
    
    if ( prtlvl > PRINT_MIN ) {
        clock_t end_time = clock();
        double solve_time = (double)(end_time - start_time)/(double)(CLOCKS_PER_SEC);
        printf("PARDISO costs %f seconds.\n", solve_time);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

/**
 * \fn INT fasp_pardiso_free_internal_mem (Pardiso_data *pdata)
 * \brief Free internal solver memory for PARDISO
 *
 * \param  pdata  Pointer to the numerical factorization data
 *
 * \author Hongxuan Zhang
 * \date   11/28/2015
 */
INT fasp_pardiso_free_internal_mem (Pardiso_data *pdata)
{
    INT status = FASP_SUCCESS;
    
    MKL_INT n, m, nnz;
    MKL_INT *ia = NULL;
    MKL_INT *ja = NULL;
    REAL *a = NULL;
    
    REAL * f = NULL;      /* RHS vector. */
    REAL * x = NULL;      /* Solution vector. */
    double ddum;          /* Double dummy */
    
    MKL_INT idum;         /* Integer dummy. */
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    MKL_INT phase, error, msglvl;    /* Auxiliary variables. */
    
    msglvl = 0; /* Do not print statistical information in file */
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif
    
    phase = -1;           /* Release internal memory. */
    PARDISO (pdata->pt, &(pdata->maxfct), &(pdata->mnum), &(pdata->mtype), &phase,
             &idum, a, ia, ja, &idum, &nrhs, pdata->iparm, &msglvl, &ddum, &ddum, &error);
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
    return status;
}

#endif

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
