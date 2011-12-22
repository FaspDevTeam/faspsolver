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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
