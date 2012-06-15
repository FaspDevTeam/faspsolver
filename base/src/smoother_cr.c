/*! \file smoother_cr.c
 *  \brief Smoothers for sparse matrix in CSR format for CR
 *
 *  \note Restricted-smoothers for compatible relaxation, C/F smoothing, etc.
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dcsr_gscr (INT pt, INT n, REAL *u, INT *ia, INT *ja, REAL *a, REAL *b, INT L, INT *CF)
 *
 * \brief Gauss Seidel method restriced to a block 
 *
 * \param pt   Relax type, e.g., cpt, fpt, etc.. 
 * \param n    Number of variables
 * \param u    Iterated solution 
 * \param ia   Row pointer
 * \param ja   Column index
 * \param a    Pointers to sparse matrix values in CSR format
 * \param b    Pointer to right hand side -- remove later also as MG relaxation on error eqn
 * \param L    Number of iterations
 * \param CF   Marker for C, F points 
 *
 * \author James Brannick
 * \date   09/07/2010
 *
 * \note Gauss Seidel CR smoother (Smoother_Type = 99)
 */
void fasp_smoother_dcsr_gscr (INT pt, 
                              INT n,
                              REAL *u,
                              INT *ia,
                              INT *ja, 
                              REAL *a, 
                              REAL *b, 
                              INT L, 
                              INT *CF)
{ 
    INT i,j,k,l;
    REAL t, d;
    
    for (l=0;l<L;++l) {
        for (i=0;i<n;++i) {
            if (CF[i] == pt) { 
                t=b[i];
                for (k=ia[i];k<ia[i+1];++k) {
                    j=ja[k];
                    if (CF[j] == pt) {
                        if (i!=j) {
                            t-=a[k]*u[j]; 
                        }
                        else { 
                            d=a[k];
                        }
                        if (ABS(d)>SMALLREAL) { 
                            u[i]=t/d;
                        }
                        else {
                            printf("### ERROR: Diagonal entry (%d,%e) is zero!\n",i,d);
                            exit(ERROR_MISC);
                        }
                    }
                }
            }
            else {
                u[i]=0.e0;
            }
        } 
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
