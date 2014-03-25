/*! \file amlirecur.c
 *
 *  \brief Abstract AMLI multilevel iteration -- recursive version
 *
 *  \note Contains AMLI and nonlinear AMLI cycles
 *
 *  TODO: Need to add a non-recursive version! --Chensong
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "forts_ns.h"

#include "mg_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_amli (AMG_data *mgl, AMG_param *param, INT level)
 *
 * \brief Solve Ax=b with recursive AMLI-cycle
 *
 * \param mgl    Pointer to AMG data: AMG_data
 * \param param  Pointer to AMG parameters: AMG_param
 * \param level  Current level
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 *
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 */
void fasp_solver_amli (AMG_data *mgl,
                       AMG_param *param,
                       INT level)
{
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const SHORT  degree= param->amli_degree;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    // local variables
    REAL   alpha  = 1.0;
    REAL * coef   = param->amli_coef;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A_level0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A_level0->row, m1 = A_level1->row;
    
    INT      *ordering = mgl[level].cfmark.val; // smoother ordering
    ILU_data *LU_level = &mgl[level].LU;        // fine level ILU decomposition
    REAL     *r        = mgl[level].w.val;      // work array for residual
    REAL     *r1       = mgl[level+1].w.val+m1; // work array for residual
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amli ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    if (prtlvl>=PRINT_MOST)
        printf("AMLI level %d, pre-smoother %d.\n", level, smoother);
    
    if (level < mgl[level].num_levels-1) {
        
        // pre smoothing
        if (level<param->ILU_levels) {
            fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
        }
        else if (level<mgl->schwarz_levels){
            switch (mgl[level].schwarz.schwarz_type){
                case 3: // TODO: Need to give a name instead of 3 --Chensong
                    fbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
                    bbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
                    break;
                default:
                    fbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
                    break;
            }
        }
        
        else {
            fasp_dcsr_presmoothing(smoother,A_level0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,ndeg,smooth_order,ordering);
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dcsr_aAxpy(-1.0,A_level0,e0->val,r);
        
        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                fasp_blas_dcsr_mxv_agg(&mgl[level].R, r, b1->val); break;
            default:
                fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val); break;
        }
        
        // coarse grid correction
        {
            fasp_array_cp(m1,b1->val,r1);
            
            INT i;
            for (i=1; i<=degree; i++) {
                fasp_dvec_set(m1,e1,0.0);
                fasp_solver_amli(mgl, param, level+1);
                
                // b1 = (coef[degree-i]/coef[degree])*r1 + A_level1*e1;
                // First, compute b1 = A_level1*e1
                fasp_blas_dcsr_mxv(A_level1, e1->val, b1->val);
                // Then, compute b1 = b1 + (coef[degree-i]/coef[degree])*r1
                fasp_blas_array_axpy(m1, coef[degree-i]/coef[degree], r1, b1->val);
            }
            
            fasp_dvec_set(m1,e1,0.0);
            fasp_solver_amli(mgl, param, level+1);
        }
        
        // find the optimal scaling factor alpha
        fasp_blas_array_ax(m1, coef[degree], e1->val);
        if ( param->coarse_scaling == ON ) {
            alpha = fasp_blas_array_dotprod(m1, e1->val, r1)
            / fasp_blas_dcsr_vmv(A_level1, e1->val, e1->val);
        }
        
        // prolongation e0 = e0 + alpha * P * e1
        switch (amg_type) {
            case UA_AMG:
                fasp_blas_dcsr_aAxpy_agg(alpha, &mgl[level].P, e1->val, e0->val); break;
            default:
                fasp_blas_dcsr_aAxpy(alpha, &mgl[level].P, e1->val, e0->val); break;
        }
        
        // post smoothing
        if (level < param->ILU_levels) {
            fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
        }
        else if (level<mgl->schwarz_levels) {
			switch (mgl[level].schwarz.schwarz_type) {
                case 3:
                    bbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
                    fbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
                    break;
                default:
                    bbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
                    break;
			}
		}
        
        else {
            fasp_dcsr_postsmoothing(smoother,A_level0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,ndeg,smooth_order,ordering);
        }
        
    }
    
    // coarsest level solver
    else {
#if   WITH_MUMPS
        /* use MUMPS direct solver on the coarsest level */
        fasp_solver_mumps(A_level0, b0, e0, 0);
#elif WITH_UMFPACK
        /* use UMFPACK direct solver on the coarsest level */
        fasp_solver_umfpack(A_level0, b0, e0, 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(A_level0, b0, e0, 0);
#else
        /* use iterative solver on the coarest level */
        fasp_coarse_itsolver(A_level0, b0, e0, tol, prtlvl);
#endif
    }
    
    if (prtlvl>=PRINT_MOST)
        printf("AMLI level %d, post-smoother %d.\n", level, smoother);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_amli ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_solver_nl_amli (AMG_data *mgl, AMG_param *param,
 *                               INT level, INT num_levels)
 *
 * \brief Solve Ax=b with recursive nonlinear AMLI-cycle
 *
 * \param mgl         Pointer to AMG_data data
 * \param param       Pointer to AMG parameters
 * \param level       Current level
 * \param num_levels  Total numebr of levels
 *
 * \author Xiaozhe Hu
 * \date   04/06/2010
 *
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 */
void fasp_solver_nl_amli (AMG_data *mgl,
                          AMG_param *param,
                          INT level,
                          INT num_levels)
{
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A_level0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A_level0->row, m1 = A_level1->row;
    
    INT      *ordering = mgl[level].cfmark.val; // smoother ordering
    ILU_data *LU_level = &mgl[level].LU;        // fine level ILU decomposition
    REAL     *r        = mgl[level].w.val;      // work array for residual
    
    dvector uH, bH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;
    bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    if (prtlvl>=PRINT_MOST)
        printf("Nonlinear AMLI level %d, pre-smoother %d.\n", num_levels, smoother);
    
    if (level < num_levels-1) {
        
        // pre smoothing
        if (level<param->ILU_levels) {
            fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
        }
        else if (level<mgl->schwarz_levels){
			switch (mgl[level].schwarz.schwarz_type){
				case 3:
					fbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
					bbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
					break;
				default:
					fbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
					break;
			}
		}
        
        else {
            fasp_dcsr_presmoothing(smoother,A_level0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,ndeg,smooth_order,ordering);
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dcsr_aAxpy(-1.0,A_level0,e0->val,r);
        
        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                fasp_blas_dcsr_mxv_agg(&mgl[level].R, r, b1->val);
                break;
            default:
                fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val);
                break;
        }
        
        // call nonlinear AMLI-cycle recursively
        {
            fasp_dvec_set(m1,e1,0.0);
            
            // The coarsest problem is solved exactly.
            // No need to call krylov method on second coarest level
            if (level == num_levels-2) {
                fasp_solver_nl_amli(&mgl[level+1], param, 0, num_levels-1);
            }
            else {  // recursively call preconditioned Krylov method on coarse grid
                precond_data pcdata;
                
                fasp_param_amg_to_prec(&pcdata, param);
                pcdata.maxit = 1;
                pcdata.max_levels = num_levels-1;
                pcdata.mgl_data = &mgl[level+1];
                
                precond pc;
                pc.data = &pcdata;
                pc.fct = fasp_precond_nl_amli;
                
                fasp_array_cp (m1, b1->val, bH.val);
                fasp_array_cp (m1, e1->val, uH.val);
                
                const INT maxit = param->amli_degree+1;
                
                switch (param->nl_amli_krylov_type) {
                    case SOLVER_GCG: // Use GCG
                        fasp_solver_dcsr_pgcg(A_level1,&bH,&uH,&pc,param->tol,
                                              maxit,1,PRINT_NONE);
                        break;
                    default: // Use FGMRES
                        fasp_solver_dcsr_pvfgmres(A_level1,&bH,&uH,&pc,param->tol,
                                                  maxit,30,1,PRINT_NONE);
                        break;
                }
                
                fasp_array_cp (m1, bH.val, b1->val);
                fasp_array_cp (m1, uH.val, e1->val);
            }
            
        }
        
        // prolongation e0 = e0 + P*e1
        switch (amg_type) {
            case UA_AMG:
                fasp_blas_dcsr_aAxpy_agg(1.0, &mgl[level].P, e1->val, e0->val);
                break;
            default:
                fasp_blas_dcsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
                break;
        }
        
        // post smoothing
        if (level < param->ILU_levels) {
            fasp_smoother_dcsr_ilu(A_level0, b0, e0, LU_level);
        }
        else if (level<mgl->schwarz_levels){
			switch (mgl[level].schwarz.schwarz_type){
				case 3:
					bbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
					fbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
					break;
				default:
					bbgs2ns_(&(mgl[level].schwarz.A.row),
                             mgl[level].schwarz.A.IA,
                             mgl[level].schwarz.A.JA,
                             mgl[level].schwarz.A.val,
                             mgl[level].x.val,
                             mgl[level].b.val,
                             &(mgl[level].schwarz.nblk),
                             mgl[level].schwarz.iblock,
                             mgl[level].schwarz.jblock,
                             mgl[level].schwarz.mask,
                             mgl[level].schwarz.maxa,
                             mgl[level].schwarz.au,
                             mgl[level].schwarz.al,
                             mgl[level].schwarz.rhsloc,
                             &(mgl[level].schwarz.memt));
					break;
			}
		}
        
        else {
            fasp_dcsr_postsmoothing(smoother,A_level0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,ndeg,smooth_order,ordering);
        }
        
    }
    
    else // coarsest level solver
    {
#if   WITH_MUMPS
        /* use MUMPS direct solver on the coarsest level */
        fasp_solver_mumps(A_level0, b0, e0, 0);
#elif WITH_UMFPACK
        /* use UMFPACK direct solver on the coarsest level */
        fasp_solver_umfpack(A_level0, b0, e0, 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(A_level0, b0, e0, 0);
#else
        /* use iterative solver on the coarest level */
        fasp_coarse_itsolver(A_level0, b0, e0, tol, prtlvl);
#endif
    }
    
    if (prtlvl>=PRINT_MOST)
        printf("Nonlinear AMLI level %d, post-smoother %d.\n", level, smoother);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_solver_nl_amli_bsr (AMG_data_bsr *mgl, AMG_param *param,
 *                                   INT level, INT num_levels)
 *
 * \brief Solve Ax=b with recursive nonlinear AMLI-cycle
 *
 * \param mgl         Pointer to AMG data: AMG_data
 * \param param       Pointer to AMG parameters: AMG_param
 * \param level       Current level
 * \param num_levels  Total numebr of levels
 *
 * \author Xiaozhe Hu
 * \date   04/06/2010
 *
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 */
void fasp_solver_nl_amli_bsr (AMG_data_bsr *mgl,
                              AMG_param *param,
                              INT level,
                              INT num_levels)
{
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol;
    INT i;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dBSRmat *A_level0 = &mgl[level].A; // fine level matrix
    dBSRmat *A_level1 = &mgl[level+1].A; // coarse level matrix
    const INT m0 = A_level0->ROW*A_level0->nb, m1 = A_level1->ROW*A_level1->nb;
    
    ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
    REAL *r = mgl[level].w.val; // for residual
    
    dvector uH, bH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;
    bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.ROW, mgl[0].A.COL, mgl[0].A.NNZ);
#endif
    
    if (prtlvl>=PRINT_MOST)
        printf("Nonlinear AMLI level %d, pre-smoother %d.\n", level, smoother);
    
    if (level < num_levels-1) {
        
        // pre smoothing
        if (level<param->ILU_levels) {
            fasp_smoother_dbsr_ilu(A_level0, b0, e0, LU_level);
        }
        else {
            SHORT steps = param->presmooth_iter;
            
            if (steps > 0){
                switch (smoother) {
                    case SMOOTHER_JACOBI:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_jacobi (A_level0, b0, e0);
                        break;
                    case SMOOTHER_GS:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_gs (A_level0, b0, e0, ASCEND, NULL);
                        break;
                    case SMOOTHER_SOR:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_sor (A_level0, b0, e0, ASCEND, NULL,relax);
                        break;
                    default:
                        printf("### ERROR: Wrong smoother type %d!\n", smoother);
                        exit(ERROR_INPUT_PAR);
                }
            }
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dbsr_aAxpy(-1.0,A_level0,e0->val,r);
        
        fasp_blas_dbsr_mxv(&mgl[level].R, r, b1->val);
        
        // call nonlinear AMLI-cycle recursively
        {
            fasp_dvec_set(m1,e1,0.0);
            
            // The coarsest problem is solved exactly.
            // No need to call krylov method on second coarest level
            if (level == num_levels-2) {
                fasp_solver_nl_amli_bsr(&mgl[level+1], param, 0, num_levels-1);
            }
            else {  // recursively call preconditioned Krylov method on coarse grid
                precond_data_bsr pcdata;
                
                fasp_param_amg_to_prec_bsr (&pcdata, param);
                pcdata.maxit = 1;
                pcdata.max_levels = num_levels-1;
                pcdata.mgl_data = &mgl[level+1];
                
                precond pc;
                pc.data = &pcdata;
                pc.fct = fasp_precond_dbsr_nl_amli;
                
                fasp_array_cp (m1, b1->val, bH.val);
                fasp_array_cp (m1, e1->val, uH.val);
                
                const INT maxit = param->amli_degree+1;
                
                // switch (param->nl_amli_krylov_type) {
                // default: // only FGMRES is supported so far --Chensong
                fasp_solver_dbsr_pvfgmres(A_level1,&bH,&uH,&pc,param->tol,
                                          maxit,MIN(maxit,30),1,PRINT_NONE);
                //     break;
                // }
                
                fasp_array_cp (m1, bH.val, b1->val);
                fasp_array_cp (m1, uH.val, e1->val);
            }
            
        }
        
        fasp_blas_dbsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
        
        // post smoothing
        if (level < param->ILU_levels) {
            fasp_smoother_dbsr_ilu(A_level0, b0, e0, LU_level);
        }
        else {
            SHORT steps = param->postsmooth_iter;
            
            if (steps > 0){
                switch (smoother) {
                    case SMOOTHER_JACOBI:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_jacobi (A_level0, b0, e0);
                        break;
                    case SMOOTHER_GS:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_gs(A_level0, b0, e0, ASCEND, NULL);
                        break;
                    case SMOOTHER_SOR:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_sor(A_level0, b0, e0, ASCEND, NULL,relax);
                        break;
                    default:
                        printf("### ERROR: Wrong smoother type!\n");
                        exit(ERROR_INPUT_PAR);
                }
            }
        }
        
    }
    
    else // coarsest level solver
    {
#if   WITH_MUMPS
        /* use MUMPS direct solver on the coarsest level */
        fasp_solver_mumps(&mgl[level].Ac, b0, e0, 0);
#elif WITH_UMFPACK
        /* use UMFPACK direct solver on the coarsest level */
        fasp_solver_umfpack(&mgl[level].Ac, b0, e0, 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(&mgl[level].Ac, b0, e0, 0);
#else
        /* use iterative solver on the coarest level */
        fasp_coarse_itsolver(&mgl[level].Ac, b0, e0, tol, prtlvl);
#endif
    }
    
    if (prtlvl>=PRINT_MOST)
        printf("Nonlinear AMLI level %d, post-smoother %d.\n", level, smoother);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_nl_amli ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_amg_amli_coef (const REAL lambda_max, const REAL lambda_min,
 *                              const INT degree, REAL *coef)
 *
 * \brief Compute the coefficients of the polynomial used by AMLI-cycle
 *
 * \param lambda_max  Maximal lambda
 * \param lambda_min  Minimal lambda
 * \param degree      Degree of polynomial approximation
 * \param coef        Coefficient of AMLI (output)
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 */
void fasp_amg_amli_coef (const REAL lambda_max,
                         const REAL lambda_min,
                         const INT degree,
                         REAL *coef)
{
    const REAL mu0 = 1.0/lambda_max, mu1 = 1.0/lambda_min;
    const REAL c = (sqrt(mu0)+sqrt(mu1))*(sqrt(mu0)+sqrt(mu1));
    const REAL a = (4*mu0*mu1)/(c);
    
    const REAL kappa = lambda_max/lambda_min; // condition number
    const REAL delta = (sqrt(kappa) - 1.0)/(sqrt(kappa)+1.0);
    const REAL b = delta*delta;
    
    if (degree == 0) {
        coef[0] = 0.5*(mu0+mu1);
    }
    
    else if (degree == 1) {
        coef[0] = 0.5*c;
        coef[1] = -1.0*mu0*mu1;
    }
    
    else if (degree > 1) {
        INT i;
        
        // allocate memory
        REAL *work = (REAL *)fasp_mem_calloc(2*degree-1, sizeof(REAL));
        REAL *coef_k, *coef_km1;
        coef_k = work; coef_km1 = work+degree;
        
        // get q_k
        fasp_amg_amli_coef(lambda_max, lambda_min, degree-1, coef_k);
        // get q_km1
        fasp_amg_amli_coef(lambda_max, lambda_min, degree-2, coef_km1);
        
        // get coef
        coef[0] = a - b*coef_km1[0] + (1+b)*coef_k[0];
        
        for (i=1; i<degree-1; i++) {
            coef[i] = -b*coef_km1[i] + (1+b)*coef_k[i] - a*coef_k[i-1];
        }
        
        coef[degree-1] = (1+b)*coef_k[degree-1] - a*coef_k[degree-2];
        
        coef[degree] = -a*coef_k[degree-1];
        
        // clean memory
        if (work) fasp_mem_free(work);
    }
    
    else {
        printf("### ERROR: Wrong degree number %d for AMLI polynomial!\n", degree);
        exit(ERROR_INPUT_PAR);
    }
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
