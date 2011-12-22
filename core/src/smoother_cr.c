/*! \file smoother_cr.c
 *  \brief Smoothers for sparse matrix in CSR format for CR
 *
 *  \note Restricted-smoothers for compatible relaxation,
 *         C/F smoothing, etc...
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn void fasp_smoother_dcsr_gscr (int pt, int n, double *u, int *ia, int *ja, double *a, double *b, int L, int *CF)
 * \brief Gauss Seidel method restriced to a block 
 * \param pt gives relax type, e.g., cpt, fpt, etc.. 
 * \param u iterated solution 
 * \param *ia, *ja, *a pointers to stiffness matrix in CSR format
 * \param *b pointer to right hand side -- remove later 
 *  also as MG relaxation on error eqn
 * \param L number of iterations
 * \return void
 *
 * \author James Brannick
 * \date   09/07/2010
 *
 * Gauss Seidel CR smoother (Smoother_Type = 99)
 */
void fasp_smoother_dcsr_gscr (int pt, 
															int n,
															double *u,
															int *ia,
															int *ja, 
															double *a, 
															double *b, 
															int L, 
															int *CF)
{ 
  int i,j,k,l;
  double t, d;
	
  for (l=0;l<L;++l){
    for (i=0;i<n;++i){
      if (CF[i] == pt) { 
        t=b[i];
       	for (k=ia[i];k<ia[i+1];++k){
					j=ja[k];
					if (CF[j] == pt){
						if (i!=j){
							t-=a[k]*u[j]; 
						}else{ 
							d=a[k];
						}
						if (ABS(d)>SMALLREAL){ 
							u[i]=t/d;
						}else{
							printf("Error: diagonal entry (%d,%e) is zero!\n",i,d);
							exit(ERROR_MISC);
						}
					}
        }
      }else{
        u[i]=0.e0;
      }
    } 
  }
}

