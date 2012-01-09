/*! \file message.c
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
 * \param ptrlvl     level for output
 * \param stop_type  type of stopping criteria
 * \param relres     relative residual of different kinds
 * \param absres     absolute residual of different kinds
 * \param factor     contraction factor
 *
 * \author Chensong Zhang
 * \date 11/16/2009
 */
void print_itinfo (const INT ptrlvl, 
                   const INT stop_type, 
                   const INT iter, 
                   const REAL relres, 
                   const REAL absres, 
                   const REAL factor)
{												
	if (ptrlvl>PRINT_MIN) {
        
		if (iter>1) {
			printf("%6d | %15.6e   | %13.6e  | %10.3f\n",iter,relres,absres,factor);
		}
		else { 
            // iter = 0 means initial guess, iter = 1 is the first iteration
			printf("---------------------------------------------------------------\n");
			switch (stop_type) {
				case STOP_REL_RES:
					printf("It Num |    ||r||/||b||    |     ||r||      |  Conv. Factor\n");
					break;
				case STOP_REL_PRECRES:
					printf("It Num |  ||r||_B/||b||_B  |     ||r||      |  Conv. Factor\n");
					break;
				case STOP_MOD_REL_RES:
					printf("It Num |    ||r||/||x||    |     ||r||      |  Conv. Factor\n");
					break;
				default:
					printf("Error: wrong stopping criteria!\n");
					exit(ERROR_INPUT_PAR);
			}
			printf("---------------------------------------------------------------\n");
			printf("%6d | %15.6e   | %13.6e  | \n",iter,relres,absres);
		} // end if iter
        
	} // end if ptrlvl
}		

/**
 * \fn void print_message (const INT ptrlvl, const char *message)
 * \brief Print output information if necessary 
 *
 * \param ptrlvl   level for output
 * \param message  error message to print
 *
 * \author Chensong Zhang
 * \date 11/16/2009
 */
void print_message (const INT ptrlvl, 
                    const char *message)
{												
	if (ptrlvl>PRINT_NONE) printf("%s",message);
}		

/**
 * \fn void fasp_chkerr (const SHORT status, const char *fctname)
 * \brief Check error status and print out error messages before quit 
 *
 * \param status   error status
 * \param fctname  function name where this routine is called
 *
 * \author Chensong Zhang
 * \date 01/10/2012
 */
void fasp_chkerr (const SHORT status, 
                  const char *fctname)
{												
    if (status>=0) return;

    switch (status) {
        case ERROR_OPEN_FILE:
            printf("### ERROR: %s -- Cannot open file!!!\n", fctname);
            break;            
        case ERROR_WRONG_FILE:
            printf("### ERROR: %s -- Wrong input file!!!\n", fctname);
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
        case ERROR_DATA_STRUCTURE:
            printf("### ERROR: %s -- Data structure mismatch!!!\n", fctname);
            break;
        case ERROR_DATA_ZERODIAG:
            printf("### ERROR: %s -- Matrix has zero diagonal entries!!!\n", fctname);
            break;
        case ERROR_DUMMY_VAR:
            printf("### ERROR: %s -- Unexpected input argument!!!\n", fctname);
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
            printf("### ERROR: %s -- Tolerance is too small for the solver!!!\n", fctname);
            break;    
        case ERROR_SOLVER_ILUSETUP:
            printf("### ERROR: %s -- ILU setup failed!!!\n", fctname);
            break;     
        case ERROR_SOLVER_MAXIT:
            printf("### ERROR: %s -- Maximal iteration number reached!!!\n", fctname);
            break;    
        case ERROR_SOLVER_EXIT:
            printf("### ERROR: %s -- Solver exited unexpected!!!\n", fctname);
            break;        
        case ERROR_SOLVER_MISC:
            printf("### ERROR: %s -- Unknown solver runtime error occurred!!!\n", fctname);
            break;   
        case ERROR_MISC:
            printf("### ERROR: %s -- Unknown error occurred!!!\n", fctname);
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
