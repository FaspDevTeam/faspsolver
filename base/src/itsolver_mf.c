/*! \file itsolver_mf.c
 *
 *  \brief Iterative solvers with matrix-free spmv
 */

#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"
#include "fasp_block.h"

static void fasp_blas_mxv_csr (void *A, REAL *x, REAL *y);
static void fasp_blas_mxv_bsr (void *A, REAL *x, REAL *y);
static void fasp_blas_mxv_str (void *A, REAL *x, REAL *y);
static void fasp_blas_mxv_bcsr (void *A, REAL *x, REAL *y);
static void fasp_blas_mxv_bbsr (void *A, REAL *x, REAL *y);
static void fasp_blas_mxv_csrl (void *A, REAL *x, REAL *y);

#include "itsolver_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_solver_itsolver (mxv_matfree *mf, dvector *b, dvector *x, 
 *                               precond *pc, itsolver_param *itparam)
 *
 * \brief Solve Ax=b by preconditioned Krylov methods for CSR matrices
 *
 * \param mf       Pointer to mxv_matfree matrix-free spmv operation
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param pc       Pointer to the preconditioning action
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang
 * \date   09/25/2009 
 *
 * \note This is an abstract interface for iterative methods.
 *
 * Modified by Feiteng Huang on 09/19/2012: matrix free
 */
INT fasp_solver_itsolver (mxv_matfree *mf, 
                          dvector *b, 
                          dvector *x, 
                          precond *pc, 
                          itsolver_param *itparam)
{
    const SHORT  print_level   = itparam->print_level;    
    const SHORT  itsolver_type = itparam->itsolver_type;
    const SHORT  stop_type     = itparam->stop_type;
    const SHORT  restart       = itparam->restart;
    const INT    MaxIt         = itparam->maxit;
    const REAL   tol           = itparam->tol; 
    
    /* Local Variables */
    REAL solver_start, solver_end, solver_duration;
    INT iter;
    
    fasp_gettime(&solver_start);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_itsolver ...... [Start]\n");
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    /* Safe-guard checks on parameters */
    ITS_CHECK ( MaxIt, tol );
    
    /* Choose a desirable Krylov iterative solver */
    switch ( itsolver_type ) {
        case SOLVER_CG:
            if (print_level>PRINT_NONE) printf("\nCalling PCG solver (matrix-free) ...\n");
            iter = fasp_solver_pcg(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        case SOLVER_BiCGstab:
            if (print_level>PRINT_NONE) printf("\nCalling BiCGstab solver (matrix-free) ...\n");
            iter = fasp_solver_pbcgs(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        case SOLVER_MinRes:
            if (print_level>PRINT_NONE) printf("\nCalling MinRes solver (matrix-free) ...\n");
            iter = fasp_solver_pminres(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        case SOLVER_GMRES:
            if (print_level>PRINT_NONE) printf("\nCalling GMRes solver (matrix-free) ...\n");
            iter = fasp_solver_pgmres(mf, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_VGMRES: 
            if (print_level>PRINT_NONE) printf("\nCalling vGMRes solver (matrix-free) ...\n");
            iter = fasp_solver_pvgmres(mf, b, x, pc, tol, MaxIt, restart, stop_type, print_level);    
            break;
            
        case SOLVER_VFGMRES: 
            if (print_level>PRINT_NONE) printf("\nCalling vFGMRes solver (matrix-free) ...\n");
            iter = fasp_solver_pvfgmres(mf, b, x, pc, tol, MaxIt, restart, stop_type, print_level);
            break;
            
        case SOLVER_GCG:
            if (print_level>PRINT_NONE) printf("\nCalling GCG solver (matrix-free) ...\n");
            iter = fasp_solver_pgcg(mf, b, x, pc, tol, MaxIt, stop_type, print_level); 
            break;
            
        default:
            printf("### ERROR: Unknown itertive solver type %d!\n", itsolver_type);
            return ERROR_SOLVER_TYPE;
            
    } 
    
    if ( (print_level>=PRINT_SOME) && (iter >= 0) ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Iterative method", solver_duration);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_itsolver ...... [Finish]\n");
#endif
    
    return iter;
}    

/**
 * \fn INT fasp_solver_krylov (mxv_matfree *mf, dvector *b, dvector *x, 
 *                             itsolver_param *itparam)
 *
 * \brief Solve Ax=b by standard Krylov methods -- without preconditioner 
 *
 * \param mf       Pointer to mxv_matfree matrix-free spmv operation
 * \param b        Pointer to the right hand side in dvector format
 * \param x        Pointer to the approx solution in dvector format
 * \param itparam  Pointer to parameters for iterative solvers
 *
 * \return         Number of iterations if succeed
 *
 * \author Chensong Zhang, Shiquan Zhang
 * \date   09/25/2009 
 *
 * Modified by Feiteng Huang on 09/20/2012: matrix free
 */
INT fasp_solver_krylov (mxv_matfree *mf, 
                        dvector *b, 
                        dvector *x, 
                        itsolver_param *itparam)
{
    const SHORT print_level = itparam->print_level;
    
    /* Local Variables */
    INT      status = SUCCESS;
    REAL     solver_start, solver_end, solver_duration;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_krylov ...... [Start]\n");
    printf("### DEBUG: matrix size: %d %d %d\n", A->row, A->col, A->nnz);
    printf("### DEBUG: rhs/sol size: %d %d\n", b->row, x->row);    
#endif
    
    fasp_gettime(&solver_start);
    
    status = fasp_solver_itsolver(mf,b,x,NULL,itparam);
    
    if ( print_level>=PRINT_MIN ) {
        fasp_gettime(&solver_end);
        solver_duration = solver_end - solver_start;
        print_cputime("Krylov method totally", solver_duration);
    }    
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_krylov ...... [Finish]\n");
#endif
    
    return status;
}

/**
 * \fn void fasp_solver_itsolver_init (INT matrix_format, mxv_matfree *mf, void *A)
 *
 * \brief Initialize itsovlers
 *
 * \param matrix_format    matrix format
 * \param mf               Pointer to mxv_matfree matrix-free spmv operation
 * \param A                void pointer to matrix
 *
 * \author Feiteng Huang
 * \date   09/18/2012
 *
 * Modified by Chensong Zhang on 05/10/2013: Change interface of mat-free mv
 */
void fasp_solver_itsolver_init (INT matrix_format,
                                mxv_matfree *mf,
                                void *A)
{
    switch ( matrix_format ) {
            
        case MAT_CSR:
            mf->fct = fasp_blas_mxv_csr;
            break;
            
        case MAT_BSR:
            mf->fct = fasp_blas_mxv_bsr;
            break;
            
        case MAT_STR:
            mf->fct = fasp_blas_mxv_str;
            break;
            
        case MAT_bCSR:
            mf->fct = fasp_blas_mxv_bcsr;
            break;
            
        case MAT_bBSR:
            mf->fct = fasp_blas_mxv_bbsr;
            break;
            
        case MAT_CSRL:
            mf->fct = fasp_blas_mxv_csrl;
            break;
            
        default:
            printf("### ERROR: Wrong matrix format %d!\n", matrix_format);
            exit(ERROR_DATA_STRUCTURE);
            
    }
    
    mf->data = A;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void fasp_blas_mxv_csr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to CSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_csr (void *A,
                               REAL *x,
                               REAL *y)
{
    dCSRmat *a = (dCSRmat *)A;
    fasp_blas_dcsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_bsr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to BSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_bsr (void *A,
                               REAL *x,
                               REAL *y)
{
    dBSRmat *a = (dBSRmat *)A;
    fasp_blas_dbsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_str (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to STR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_str (void *A,
                               REAL *x,
                               REAL *y)
{
    dSTRmat *a = (dSTRmat *)A;
    fasp_blas_dstr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_bcsr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to bCSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_bcsr (void *A,
                                REAL *x,
                                REAL *y)
{
    block_dCSRmat *a = (block_dCSRmat *)A;
    fasp_blas_bdcsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_bbsr (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to bBSR matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_bbsr (void *A,
                                REAL *x,
                                REAL *y)
{
    block_BSR *a = (block_BSR *)A;
    fasp_blas_bdbsr_mxv(a, x, y);
}

/**
 * \fn static void fasp_blas_mxv_csrl (void *A, REAL *x, REAL *y)
 *
 * \brief Matrix-vector multiplication y = A*x
 *
 * \param A               Pointer to CSRL matrix A
 * \param x               Pointer to array x
 * \param y               Pointer to array y
 *
 * \author Feiteng Huang
 * \date   09/19/2012
 */
static void fasp_blas_mxv_csrl (void *A,
                                REAL *x,
                                REAL *y)
{
    dCSRLmat *a = (dCSRLmat *)A;
    fasp_blas_dcsrl_mxv(a, x, y);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
