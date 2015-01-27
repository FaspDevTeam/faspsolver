/*! \file mgrecur.c
 *
 *  \brief Abstract multigrid cycle -- recursive version
 * 
 *  \note  Not used any more. Will be removed! --Chensong
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "mg_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_mgrecur (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive multigrid K-cycle
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 * \param level  Index of the current level
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   04/06/2010
 *
 * Modified by Chensong Zhang on 01/10/2012
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 */
void fasp_solver_mgrecur (AMG_data *mgl,
                          AMG_param *param,
                          INT level)
{
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  smooth_order = param->smooth_order;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A; // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    const INT m0 = A0->row, m1 = A1->row;
    
    ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
    REAL *r = mgl[level].w.val; // for residual
    INT *ordering = mgl[level].cfmark.val; // for smoother ordering
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: n=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.nnz);
#endif
    
    if ( prtlvl >= PRINT_MOST )
        printf("AMG level %d, smoother %d.\n", level, smoother);
    
    if ( level < mgl[level].num_levels-1 ) {
        
        // pre smoothing
        if (level<param->ILU_levels) {
            fasp_smoother_dcsr_ilu(A0, b0, e0, LU_level);
        }
        else {
            fasp_dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,ndeg,smooth_order,ordering);
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dcsr_aAxpy(-1.0,A0,e0->val,r);
        
        // restriction r1 = R*r0
        fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val);
        
        { // call MG recursively: type = 1 for V cycle, type = 2 for W cycle
            SHORT i;
            fasp_dvec_set(m1,e1,0.0);
            for (i=0; i<cycle_type; ++i) fasp_solver_mgrecur (mgl, param, level+1);
        }
        
        // prolongation e0 = e0 + P*e1
        fasp_blas_dcsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
        
        // post smoothing
        if (level < param->ILU_levels) {
            fasp_smoother_dcsr_ilu(A0, b0, e0, LU_level);
        }
        else {
            fasp_dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,ndeg,smooth_order,ordering);
        }
        
    }
    
    else { // coarsest level solver
        
        switch (coarse_solver) {
                
#if WITH_SuperLU
            /* use SuperLU direct solver on the coarsest level */
            case SOLVER_SUPERLU:
                fasp_solver_superlu(A0, b0, e0, 0);
                break;
#endif
                
#if WITH_UMFPACK
            /* use UMFPACK direct solver on the coarsest level */
            case SOLVER_UMFPACK:
                //fasp_solver_umfpack(A0, b0, e0, 0);
                fasp_umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif
                
#if WITH_MUMPS
            /* use MUMPS direct solver on the coarsest level */
            case SOLVER_MUMPS:
                mgl[level].mumps.job = 2;
                fasp_solver_mumps_steps(A0, b0, e0, &mgl[level].mumps);
                break;
#endif
                
            /* use iterative solver on the coarsest level */
            default:
                fasp_coarse_itsolver(A0, b0, e0, tol, prtlvl);
                
        }
        
    }
    
#if DEBUG_MODE
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
