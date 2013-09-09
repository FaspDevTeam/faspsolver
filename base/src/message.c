/*! \file message.c
 *
 *  \brief Output some useful messages
 *
 *  \note These routines are meant for internal use only.
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void print_itinfo (const INT ptrlvl, const INT stop_type, const INT iter,
 *                        const REAL relres, const REAL absres, const REAL factor)
 *
 * \brief Print out iteration information for iterative solvers
 *
 * \param ptrlvl     Level for output
 * \param stop_type  Type of stopping criteria
 * \param iter       Number of iterations
 * \param relres     Relative residual of different kinds
 * \param absres     Absolute residual of different kinds
 * \param factor     Contraction factor
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 *
 * Modified by Chensong Zhang on 03/28/2013: Output initial guess
 * Modified by Chensong Zhang on 04/05/2013: Fix a typo
 */
void print_itinfo (const INT ptrlvl,
                   const INT stop_type,
                   const INT iter,
                   const REAL relres,
                   const REAL absres,
                   const REAL factor)
{
    if ( ptrlvl > PRINT_SOME ) {
        
        if ( iter > 0 ) {
            printf("%6d | %13.6e   | %13.6e  | %10.4f\n",iter,relres,absres,factor);
        }
        else { // iter = 0: initial guess
            printf("-----------------------------------------------------------\n");
            switch (stop_type) {
                case STOP_REL_RES:
                    printf("It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
                    break;
                case STOP_REL_PRECRES:
                    printf("It Num | ||r||_B/||b||_B |    ||r||_B     |  Conv. Factor\n");
                    break;
                case STOP_MOD_REL_RES:
                    printf("It Num |   ||r||/||x||   |     ||r||      |  Conv. Factor\n");
                    break;
            }
            printf("-----------------------------------------------------------\n");
            printf("%6d | %13.6e   | %13.6e  |     -.-- \n",iter,relres,absres);
        } // end if iter
        
    } // end if ptrlvl
}

/**
 * \fn void void print_amgcomplexity (AMG_data *mgl, const SHORT prtlvl)
 *
 * \brief Print complexities of AMG method
 *
 * \param mgl      Multilevel hierachy for AMG
 * \param prtlvl   How much information to print
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
void print_amgcomplexity (AMG_data *mgl,
                          const SHORT prtlvl)
{
    const SHORT   max_levels=mgl->num_levels;
    SHORT         level;
    REAL          gridcom=0.0, opcom=0.0;
    
    if ( prtlvl >= PRINT_SOME ) {
        printf("-----------------------------------------------\n");
        printf("  Level     Num of rows      Num of nonzeros\n");
        printf("-----------------------------------------------\n");
        for ( level = 0; level < max_levels; ++level) {
            printf("%5d  %14d  %16d\n", level, mgl[level].A.row, mgl[level].A.nnz);
            gridcom += mgl[level].A.row;
            opcom   += mgl[level].A.nnz;
        }
        printf("-----------------------------------------------\n");
        
        gridcom /= mgl[0].A.row;
        opcom   /= mgl[0].A.nnz;
        printf("AMG grid complexity     = %.3f\n", gridcom);
        printf("AMG operator complexity = %.3f\n", opcom);
    }
}

/**
 * \fn void void print_amgcomplexity_bsr (AMG_data_bsr *mgl, const SHORT prtlvl)
 *
 * \brief Print complexities of AMG method for BSR matrices
 *
 * \param mgl      Multilevel hierachy for AMG
 * \param prtlvl   How much information to print
 *
 * \author Chensong Zhang
 * \date   05/10/2013
 */
void print_amgcomplexity_bsr (AMG_data_bsr *mgl,
                              const SHORT prtlvl)
{
    const SHORT  max_levels = mgl->num_levels;
    SHORT        level;
    REAL         gridcom = 0.0, opcom = 0.0;
    
    if ( prtlvl >= PRINT_SOME ) {
        printf("-----------------------------------------------\n");
        printf("  Level     Num of rows      Num of nonzeros\n");
        printf("-----------------------------------------------\n");
        for ( level = 0; level < max_levels; ++level ) {
            printf("%5d  %14d  %16d\n", level,mgl[level].A.ROW, mgl[level].A.NNZ);
            gridcom += mgl[level].A.ROW;
            opcom   += mgl[level].A.NNZ;
        }
        printf("-----------------------------------------------\n");
        
        gridcom /= mgl[0].A.ROW;
        opcom   /= mgl[0].A.NNZ;
        printf("AMG (BSR) grid complexity     = %.3f\n", gridcom);
        printf("AMG (BSR) operator complexity = %.3f\n", opcom);
    }
}

/**
 * \fn void void print_cputime (const char *message, const REAL cputime)
 *
 * \brief Print CPU walltime
 *
 * \param message   Some string to print out
 * \param cputime   Walltime since start to end
 *
 * \author Chensong Zhang
 * \date   04/10/2012
 */
void print_cputime (const char *message,
                    const REAL cputime)
{
    printf("%s costs %.4f seconds.\n", message, cputime);
}

/**
 * \fn void print_message (const INT ptrlvl, const char *message)
 *
 * \brief Print output information if necessary
 *
 * \param ptrlvl   Level for output
 * \param message  Error message to print
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
void print_message (const INT ptrlvl,
                    const char *message)
{
    if ( ptrlvl > PRINT_NONE ) printf("%s", message);
}

/**
 * \fn void fasp_chkerr (const SHORT status, const char *fctname)
 *
 * \brief Check error status and print out error messages before quit
 *
 * \param status   Error status
 * \param fctname  Function name where this routine is called
 *
 * \author Chensong Zhang
 * \date   01/10/2012
 */
void fasp_chkerr (const SHORT status,
                  const char *fctname)
{
    if ( status >= 0 ) return; // No error at all
    
    switch ( status ) {
        case ERROR_OPEN_FILE:
            printf("### ERROR: %s -- Cannot open file!!!\n", fctname);
            break;
        case ERROR_WRONG_FILE:
            printf("### ERROR: %s -- Wrong file format!!!\n", fctname);
            break;
        case ERROR_INPUT_PAR:
            printf("### ERROR: %s -- Wrong input arguments!!!\n", fctname);
            break;
        case ERROR_REGRESS:
            printf("### ERROR: %s -- Regression test failed!!!\n", fctname);
            break;
        case ERROR_ALLOC_MEM:
            printf("### ERROR: %s -- Cannot allocate memory!!!\n", fctname);
            break;
        case ERROR_NUM_BLOCKS:
            printf("### ERROR: %s -- Wrong number of blocks!!!\n", fctname);
            break;
        case ERROR_DATA_STRUCTURE:
            printf("### ERROR: %s -- Data structure mismatch!!!\n", fctname);
            break;
        case ERROR_DATA_ZERODIAG:
            printf("### ERROR: %s -- Matrix has zero diagonal entries!!!\n", fctname);
            break;
        case ERROR_DUMMY_VAR:
            printf("### ERROR: %s -- Unexpected input argument!!!\n", fctname);
            break;
        case ERROR_AMG_INTERP_TYPE:
            printf("### ERROR: %s -- Unknown AMG interpolation type!!!\n", fctname);
            break;
        case ERROR_AMG_COARSE_TYPE:
            printf("### ERROR: %s -- Unknown AMG coarsening type!!!\n", fctname);
            break;
        case ERROR_AMG_SMOOTH_TYPE:
            printf("### ERROR: %s -- Unknown AMG smoother type!!!\n", fctname);
            break;
        case ERROR_SOLVER_TYPE:
            printf("### ERROR: %s -- Unknown solver type!!!\n", fctname);
            break;
        case ERROR_SOLVER_PRECTYPE:
            printf("### ERROR: %s -- Unknown preconditioner type!!!\n", fctname);
            break;
        case ERROR_SOLVER_STAG:
            printf("### ERROR: %s -- Solver stagnation error!!!\n", fctname);
            break;
        case ERROR_SOLVER_SOLSTAG:
            printf("### ERROR: %s -- Solution is close to zero!!!\n", fctname);
            break;
        case ERROR_SOLVER_TOLSMALL:
            printf("### ERROR: %s -- Tol is too small for the solver!!!\n", fctname);
            break;
        case ERROR_SOLVER_ILUSETUP:
            printf("### ERROR: %s -- ILU setup failed!!!\n", fctname);
            break;
        case ERROR_SOLVER_MAXIT:
            printf("### ERROR: %s -- Max iteration number reached!!!\n", fctname);
            break;
        case ERROR_SOLVER_EXIT:
            printf("### ERROR: %s -- Solver exited unexpected!!!\n", fctname);
            break;
        case ERROR_SOLVER_MISC:
            printf("### ERROR: %s -- Unknown solver runtime error!!!\n", fctname);
            break;
        case ERROR_MISC:
            printf("### ERROR: %s -- Unknown error occurred!!!\n", fctname);
            break;
        case ERROR_QUAD_TYPE:
            printf("### ERROR: %s -- Unknown quadrature rules!!!\n", fctname);
            break;
        case ERROR_QUAD_DIM:
            printf("### ERROR: %s -- Num of quad points is not supported!!!\n", fctname);
            break;
        case RUN_FAIL:
            printf("### ERROR: %s -- Function does not exit successfully!!!\n", fctname);
            break;
        default:
            break;
    }
    
    exit(status);
}    

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
