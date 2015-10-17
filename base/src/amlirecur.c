/*! \file amlirecur.c
 *
 *  \brief Abstract AMLI multilevel iteration -- recursive version
 *
 *  \note  AMLI and non-linear AMLI cycles
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
//       #include "forts_ns.h"

#include "mg_util.inl"

static SHORT fasp_krylov_cycle_dcsr_pgcg(dCSRmat *, dvector *, dvector *, precond *);
static SHORT fasp_krylov_cycle_dcsr_pgcr(dCSRmat *, dvector *, dvector *, precond *);

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
 * \note AMLI polynomial computed by the best approximation of 1/x.
 *       Refer to Johannes K. Kraus, Panayot S. Vassilevski, Ludmil T. Zikatanov,
 *       "Polynomial of best uniform approximation to $x^{-1}$ and smoothing in
 *        two-level methods", 2013.
 *
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 * Modified by Zheng Li on 11/10/2014: update direct solvers.
 */
void fasp_solver_amli (AMG_data *mgl,
                       AMG_param *param,
                       INT level)
{
    const SHORT  amg_type=param->AMG_type;
    const SHORT  prtlvl = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const SHORT  coarse_solver = param->coarse_solver;
    const SHORT  degree= param->amli_degree;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    // local variables
    REAL   alpha  = 1.0;
    REAL * coef   = param->amli_coef;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A0->row, m1 = A1->row;
    
    INT      *ordering = mgl[level].cfmark.val; // smoother ordering
    ILU_data *LU_level = &mgl[level].LU;        // fine level ILU decomposition
    REAL     *r        = mgl[level].w.val;      // work array for residual
    REAL     *r1       = mgl[level+1].w.val+m1; // work array for residual
    
    // Schwarz parameters
    Schwarz_param swzparam;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: n=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.nnz);
#endif
    
    if ( prtlvl >= PRINT_MOST )
        printf("AMLI level %d, smoother %d.\n", level, smoother);
    
    if ( level < mgl[level].num_levels-1 ) {
        
        // pre smoothing
        if ( level < mgl[level].ILU_levels ) {
            
            fasp_smoother_dcsr_ilu(A0, b0, e0, LU_level);
            
        }
        
        else if ( level < mgl->Schwarz_levels ) {
            
            switch (mgl[level].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[level].Schwarz, &swzparam,
                                                       &mgl[level].x, &mgl[level].b);
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[level].Schwarz, &swzparam,
                                                        &mgl[level].x, &mgl[level].b);
                    break;
                default:
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[level].Schwarz, &swzparam,
                                                       &mgl[level].x, &mgl[level].b);
                    break;
            }
        }
        
        else {
            fasp_dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,ndeg,smooth_order,ordering);
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dcsr_aAxpy(-1.0,A0,e0->val,r);
        
        // restriction r1 = R*r0
        switch (amg_type) {
            case UA_AMG:
                fasp_blas_dcsr_mxv_agg(&mgl[level].R, r, b1->val); break;
            default:
                fasp_blas_dcsr_mxv(&mgl[level].R, r, b1->val); break;
        }
        
        // coarse grid correction
        {
            INT i;
            
            fasp_array_cp(m1,b1->val,r1);
            
            for ( i=1; i<=degree; i++ ) {
                fasp_dvec_set(m1,e1,0.0);
                fasp_solver_amli(mgl, param, level+1);
                
                // b1 = (coef[degree-i]/coef[degree])*r1 + A1*e1;
                // First, compute b1 = A1*e1
                fasp_blas_dcsr_mxv(A1, e1->val, b1->val);
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
            / fasp_blas_dcsr_vmv(A1, e1->val, e1->val);
            alpha = MIN(alpha, 1.0);
        }
        
        // prolongation e0 = e0 + alpha * P * e1
        switch (amg_type) {
            case UA_AMG:
                fasp_blas_dcsr_aAxpy_agg(alpha, &mgl[level].P, e1->val, e0->val);
                break;
            default:
                fasp_blas_dcsr_aAxpy(alpha, &mgl[level].P, e1->val, e0->val);
                break;
        }
        
        // post smoothing
        if ( level < mgl[level].ILU_levels ) {
            
            fasp_smoother_dcsr_ilu(A0, b0, e0, LU_level);
            
        }
        
        else if (level<mgl->Schwarz_levels) {
        
            switch (mgl[level].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[level].Schwarz, &swzparam,
                                                        &mgl[level].x, &mgl[level].b);
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[level].Schwarz, &swzparam,
                                                       &mgl[level].x, &mgl[level].b);
                    break;
                default:
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[level].Schwarz, &swzparam,
                                                        &mgl[level].x, &mgl[level].b);
                    break;
            }
        }
        
        else {
            fasp_dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,ndeg,smooth_order,ordering);
        }
        
    }
    
    else { // coarsest level solver
        
        switch (coarse_solver) {
                
#if WITH_SuperLU
            case SOLVER_SUPERLU:
                /* use SuperLU direct solver on the coarsest level */
                fasp_solver_superlu(A0, b0, e0, 0);
                break;
#endif
                
#if WITH_UMFPACK
            case SOLVER_UMFPACK:
                /* use UMFPACK direct solver on the coarsest level */
                fasp_umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif
                
#if WITH_MUMPS
            case SOLVER_MUMPS:
                /* use MUMPS direct solver on the coarsest level */
                mgl[level].mumps.job = 2;
                fasp_solver_mumps_steps(A0, b0, e0, &mgl[level].mumps);
                break;
#endif
                
            default:
                /* use iterative solver on the coarsest level */
                fasp_coarse_itsolver(A0, b0, e0, tol, prtlvl);
                
        }
        
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
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
 * \param num_levels  Total number of levels
 *
 * \author Xiaozhe Hu
 * \date   04/06/2010
 *
 * \note Refer to Xiazhe Hu, Panayot S. Vassilevski, Jinchao Xu
 *       "Comparative Convergence Analysis of Nonlinear AMLI-cycle Multigrid", 2013.
 *
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 * Modified by Zheng Li on 11/10/2014: update direct solvers.
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
    const SHORT  coarse_solver = param->coarse_solver;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol*1e-4;
    const SHORT  ndeg = param->polynomial_degree;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x;   // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dCSRmat *A0 = &mgl[level].A;   // fine level matrix
    dCSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    
    const INT m0 = A0->row, m1 = A1->row;
    
    INT      *ordering = mgl[level].cfmark.val; // smoother ordering
    ILU_data *LU_level = &mgl[level].LU;        // fine level ILU decomposition
    REAL     *r        = mgl[level].w.val;      // work array for residual
    
    dvector uH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;
    
    // Schwarz parameters
    Schwarz_param swzparam;
    if ( param->Schwarz_levels > 0 ) {
        swzparam.Schwarz_blksolver = param->Schwarz_blksolver;
    }

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: n=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.nnz);
#endif
    
    if ( prtlvl >= PRINT_MOST )
        printf("Nonlinear AMLI level %d, smoother %d.\n", num_levels, smoother);
    
    if ( level < num_levels-1 ) {
        
        // pre smoothing
        if ( level < mgl[level].ILU_levels ) {

            fasp_smoother_dcsr_ilu(A0, b0, e0, LU_level);
            
        }
        
        else if ( level < mgl->Schwarz_levels ) {
            
            switch (mgl[level].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[level].Schwarz, &swzparam,
                                                       &mgl[level].x, &mgl[level].b);
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[level].Schwarz, &swzparam,
                                                        &mgl[level].x, &mgl[level].b);
                    break;
                default:
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[level].Schwarz, &swzparam,
                                                       &mgl[level].x, &mgl[level].b);
                    break;
            }
        }
        
        else {
            
            fasp_dcsr_presmoothing(smoother,A0,b0,e0,param->presmooth_iter,
                                   0,m0-1,1,relax,ndeg,smooth_order,ordering);
            
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dcsr_aAxpy(-1.0,A0,e0->val,r);
        
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
            
            // V-cycle will be enforced when needed !!!
            if ( mgl[level+1].cycle_type <= 1 ) {
            
                fasp_solver_nl_amli(&mgl[level+1], param, 0, num_levels-1);
            
            }
            
            else { // recursively call preconditioned Krylov method on coarse grid
            
                precond_data pcdata;
                
                fasp_param_amg_to_prec(&pcdata, param);
                pcdata.maxit = 1;
                pcdata.max_levels = num_levels-1;
                pcdata.mgl_data = &mgl[level+1];

                precond pc;
                pc.data = &pcdata;
                pc.fct = fasp_precond_nl_amli;

                fasp_array_cp (m1, e1->val, uH.val);
                
                switch (param->nl_amli_krylov_type) {
                        
                    case SOLVER_GCG: // Use GCG
                        fasp_krylov_cycle_dcsr_pgcg(A1,b1,&uH,&pc);
                        break;
                        
                    default: // Use GCR
                        fasp_krylov_cycle_dcsr_pgcr(A1,b1,&uH,&pc);
                        break;
                }

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
        if ( level < mgl[level].ILU_levels ) {
        
            fasp_smoother_dcsr_ilu(A0, b0, e0, LU_level);
        
        }
        else if ( level < mgl->Schwarz_levels ) {
        
            switch (mgl[level].Schwarz.Schwarz_type) {
                case SCHWARZ_SYMMETRIC:
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[level].Schwarz, &swzparam,
                                                        &mgl[level].x, &mgl[level].b);
                    fasp_dcsr_Schwarz_forward_smoother(&mgl[level].Schwarz, &swzparam,
                                                       &mgl[level].x, &mgl[level].b);
                    break;
                default:
                    fasp_dcsr_Schwarz_backward_smoother(&mgl[level].Schwarz, &swzparam,
                                                        &mgl[level].x, &mgl[level].b);
                    break;
            }
        
        }
        
        else {
            fasp_dcsr_postsmoothing(smoother,A0,b0,e0,param->postsmooth_iter,
                                    0,m0-1,-1,relax,ndeg,smooth_order,ordering);
        }
        
    }
    
    else { // coarsest level solver
        
        switch (coarse_solver) {
                
#if WITH_SuperLU
            case SOLVER_SUPERLU:
                /* use SuperLU direct solver on the coarsest level */
                fasp_solver_superlu(A0, b0, e0, 0);
                break;
#endif
                
#if WITH_UMFPACK
            case SOLVER_UMFPACK:
                /* use UMFPACK direct solver on the coarsest level */
                fasp_umfpack_solve(A0, b0, e0, mgl[level].Numeric, 0);
                break;
#endif
                
#if WITH_MUMPS
            case SOLVER_MUMPS:
                /* use MUMPS direct solver on the coarsest level */
                mgl[level].mumps.job = 2;
                fasp_solver_mumps_steps(A0, b0, e0, &mgl[level].mumps);
                break;
#endif
                
            default:
                /* use iterative solver on the coarsest level */
                fasp_coarse_itsolver(A0, b0, e0, tol, prtlvl);
                
        }
        
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
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
 * \param num_levels  Total number of levels
 *
 * \author Xiaozhe Hu
 * \date   04/06/2010
 *
 * \note Nonlinear AMLI-cycle.
 *       Refer to Xiazhe Hu, Panayot S. Vassilevski, Jinchao Xu
 *       "Comparative Convergence Analysis of Nonlinear AMLI-cycle Multigrid", 2013.
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
    const SHORT  coarse_solver = param->coarse_solver;
    const REAL   relax = param->relaxation;
    const REAL   tol = param->tol;
    INT i;
    
    dvector *b0 = &mgl[level].b,   *e0 = &mgl[level].x; // fine level b and x
    dvector *b1 = &mgl[level+1].b, *e1 = &mgl[level+1].x; // coarse level b and x
    
    dBSRmat *A0 = &mgl[level].A; // fine level matrix
    dBSRmat *A1 = &mgl[level+1].A; // coarse level matrix
    const INT m0 = A0->ROW*A0->nb, m1 = A1->ROW*A1->nb;
    
    ILU_data *LU_level = &mgl[level].LU; // fine level ILU decomposition
    REAL *r = mgl[level].w.val; // for residual
    
    dvector uH, bH;  // for coarse level correction
    uH.row = m1; uH.val = mgl[level+1].w.val + m1;
    bH.row = m1; bH.val = mgl[level+1].w.val + 2*m1;
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
    printf("### DEBUG: n=%d, nnz=%d\n", mgl[0].A.ROW, mgl[0].A.NNZ);
#endif
    
    if (prtlvl>=PRINT_MOST)
        printf("Nonlinear AMLI level %d, pre-smoother %d.\n", level, smoother);
    
    if (level < num_levels-1) {
        
        // pre smoothing
        if (level<param->ILU_levels) {
            fasp_smoother_dbsr_ilu(A0, b0, e0, LU_level);
        }
        else {
            SHORT steps = param->presmooth_iter;
            
            if (steps > 0) {
                switch (smoother) {
                    case SMOOTHER_JACOBI:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_jacobi (A0, b0, e0);
                        break;
                    case SMOOTHER_GS:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_gs (A0, b0, e0, ASCEND, NULL);
                        break;
                    case SMOOTHER_SOR:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_sor (A0, b0, e0, ASCEND, NULL,relax);
                        break;
                    default:
                        printf("### ERROR: Wrong smoother type %d!\n", smoother);
                        fasp_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
                }
            }
        }
        
        // form residual r = b - A x
        fasp_array_cp(m0,b0->val,r);
        fasp_blas_dbsr_aAxpy(-1.0,A0,e0->val,r);
        
        fasp_blas_dbsr_mxv(&mgl[level].R, r, b1->val);
        
        // call nonlinear AMLI-cycle recursively
        {
            fasp_dvec_set(m1,e1,0.0);
            
            // The coarsest problem is solved exactly.
            // No need to call Krylov method on second coarsest level
            if (level == num_levels-2) {
                fasp_solver_nl_amli_bsr(&mgl[level+1], param, 0, num_levels-1);
            }
            else { // recursively call preconditioned Krylov method on coarse grid
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
                
                fasp_solver_dbsr_pvfgmres(A1,&bH,&uH,&pc,param->tol,
                                          maxit,MIN(maxit,30),1,PRINT_NONE);
                
                fasp_array_cp (m1, bH.val, b1->val);
                fasp_array_cp (m1, uH.val, e1->val);
            }
            
        }
        
        fasp_blas_dbsr_aAxpy(1.0, &mgl[level].P, e1->val, e0->val);
        
        // post smoothing
        if (level < param->ILU_levels) {
            fasp_smoother_dbsr_ilu(A0, b0, e0, LU_level);
        }
        else {
            SHORT steps = param->postsmooth_iter;
            
            if (steps > 0) {
                switch (smoother) {
                    case SMOOTHER_JACOBI:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_jacobi (A0, b0, e0);
                        break;
                    case SMOOTHER_GS:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_gs(A0, b0, e0, ASCEND, NULL);
                        break;
                    case SMOOTHER_SOR:
                        for (i=0; i<steps; i++)
                            fasp_smoother_dbsr_sor(A0, b0, e0, ASCEND, NULL,relax);
                        break;
                    default:
                        printf("### ERROR: Wrong smoother type %d!\n", smoother);
                        fasp_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
                }
            }
        }
        
    }
    
    else { // coarsest level solver
        
        switch (coarse_solver) {
                
#if WITH_SuperLU
            case SOLVER_SUPERLU:
                /* use SuperLU direct solver on the coarsest level */
                fasp_solver_superlu(&mgl[level].Ac, b0, e0, 0);
                break;
#endif
                
#if WITH_UMFPACK
            case SOLVER_UMFPACK:
                /* use UMFPACK direct solver on the coarsest level */
                fasp_umfpack_solve(&mgl[level].Ac, b0, e0, mgl[level].Numeric, 0);
                break;
#endif
                
#if WITH_MUMPS
            case SOLVER_MUMPS:
                /* use MUMPS direct solver on the coarsest level */
                mgl[level].mumps.job = 2;
                fasp_solver_mumps_steps(&mgl[level].Ac, b0, e0, &mgl[level].mumps);
                break;
#endif
                
            default:
                /* use iterative solver on the coarsest level */
                fasp_coarse_itsolver(&mgl[level].Ac, b0, e0, tol, prtlvl);
                
        }
        
    }
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
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
        printf("### ERROR: Wrong AMLI degree %d!\n", degree);
        fasp_chkerr(ERROR_INPUT_PAR, __FUNCTION__);
    }
    
    return;
}

/*---------------------------------*/
/*--     Private Functions       --*/
/*---------------------------------*/

/**
 * \fn static SHORT fasp_krylov_cycle_dcsr_pgcg (dCSRmat *A, dvector *b,
 *                                               dvector *u, precond *pc)
 *
 * \brief A preconditioned GCR method for solving Au=b
 *
 * \param *A    Pointer to the coefficient matrix
 * \param *b    Pointer to the dvector of right hand side
 * \param *u    Pointer to the dvector of DOFs
 * \param *pre  Pointer to the structure of precondition (precond)
 *
 * \author Zheng Li, Chensong Zhang
 * \date   11/09/2014
 *
 * \note   Specified for unsmoothed aggregation cycle
 */
static SHORT fasp_krylov_cycle_dcsr_pgcg (dCSRmat *A,
                                          dvector *b,
                                          dvector *u,
                                          precond *pc)
{
    REAL   absres, relres, normb;
    REAL   alpha1, alpha2, gamma1, gamma2, rho1, rho2, beta1, beta2, beta3, beta4;
    REAL   *work, *r, *x1, *v1, *v2;
    
    INT    m=A->row;
    REAL   *x = u->val;
    
    // allocate temp memory
    work = (REAL *)fasp_mem_calloc(4*m,sizeof(REAL));
    r = work; x1 = r + m; v1 = r + 2*m; v2 = r + 3*m;
    
    normb = fasp_blas_array_norm2(m, b->val);
    
    fasp_array_cp(m, b->val, r);
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x, pc->data);
    else
        fasp_array_cp(m, r, x);
    
    // v1 = A*p
    fasp_blas_dcsr_mxv(A, x, v1);
    
    // rho1 = (p,v1)
    rho1 = fasp_blas_array_dotprod (m, x, v1);
    
    // alpha1 = (p, r)
    alpha1 = fasp_blas_array_dotprod (m, x, r);
    
    beta1 = alpha1/rho1;
    
    // r = r - beta1 *v1
    fasp_blas_array_axpy(m, -beta1, v1, r);
    
    // norm(r)
    absres = fasp_blas_array_norm2(m, r);
    
    // compute relative residual
    relres = absres/normb;
    
    // if relres reaches tol(0.2), pgcg will stop,
    // otherwise, another one pgcg iteration will do.
    if (relres < 0.2) {
        fasp_blas_array_ax(m, beta1, x);
        fasp_mem_free(work);
        return FASP_SUCCESS;
    }
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x1, pc->data);
    else
        fasp_array_cp(m, r, x1);
    
    //v2 = A*p
    fasp_blas_dcsr_mxv(A, x1, v2);
    
    //gamma0 = (x1,v1)
    gamma1 = fasp_blas_array_dotprod (m, x1, v1);
    
    //alpha2 = (x1,r)
    alpha2  = fasp_blas_array_dotprod(m, x1, r);
    
    //rho2 = (x1,v2)
    rho2 = fasp_blas_array_dotprod(m, x1, v2);
    
    gamma2 = gamma1;
    
    beta2 = rho2 - gamma1*gamma2/rho1;
    beta3 = (alpha1 - gamma2*alpha2/beta2)/rho1;
    beta4 = alpha2/beta2;
    
    fasp_blas_array_ax(m, beta3, x);
    
    fasp_blas_array_axpy(m, beta4, x1, x);
    
    // free
    fasp_mem_free(work);
    return FASP_SUCCESS;
}

/**
 * \fn static SHORT fasp_krylov_cycle_dcsr_pgcr (dCSRmat *A, dvector *b,
 *                                               dvector *u, precond *pc)
 *
 * \brief A preconditioned GCR method for solving Au=b
 *
 * \param *A    Pointer to the coefficient matrix
 * \param *b    Pointer to the dvector of right hand side
 * \param *u    Pointer to the dvector of DOFs
 * \param *pre  Pointer to the structure of precondition (precond)
 *
 * \author zheng Li, Chensong Zhang
 * \date   11/09/2014
 *
 * \note   Specified for unsmoothed aggregation cycle.
 */
static SHORT fasp_krylov_cycle_dcsr_pgcr (dCSRmat *A,
                                          dvector *b,
                                          dvector *u,
                                          precond *pc)
{
    REAL   absres = BIGREAL;
    REAL   relres  = BIGREAL, normb  = BIGREAL;
    REAL   alpha, alpha1, alpha2, alpha3, alpha4, beta, gamma, rho1, rho2;
    
    INT    m=A->row;
    REAL   *x = u->val;
    
    // allocate temp memory
    REAL *work, *r, *x1, *v1, *v2;
    work = (REAL *)fasp_mem_calloc(4*m,sizeof(REAL));
    r = work; x1 = r + m; v1 = r + 2*m; v2 = r + 3*m;
    
    normb=fasp_blas_array_norm2(m, b->val);
    fasp_array_cp(m, b->val, r);
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x, pc->data);
    else
        fasp_array_cp(m, r, x);
    
    // v1 = A*x
    fasp_blas_dcsr_mxv(A, x, v1);
    // rho1 = (v1,v1)
    rho1 = fasp_blas_array_dotprod (m, v1, v1);
    // alpha1 = (r, v1)
    alpha1 = fasp_blas_array_dotprod (m, v1, r);
    
    alpha = alpha1/rho1;
    
    // r = r - alpha *v1
    fasp_blas_array_axpy(m, -alpha, v1, r);
    
    // norm(r)
    absres = fasp_blas_array_norm2(m, r);
    
    // compute relative residual
    relres = absres/normb;
    
    // if relres reaches tol(0.2), pgcr will stop,
    // otherwise, another one pgcr iteration will do.
    if (relres < 0.2) {
        fasp_blas_array_ax(m, alpha, x);
        fasp_mem_free(work);
        return FASP_SUCCESS;
    }
    
    // Preconditioning
    if (pc != NULL)
        pc->fct(r, x1, pc->data);
    else
        fasp_array_cp(m, r, x1);
    
    //v2 = A*x1
    fasp_blas_dcsr_mxv(A, x1, v2);
    
    //gamma = (v1,v2)
    gamma = fasp_blas_array_dotprod (m, v1, v2);
    //beta = (v2,v2)
    beta  = fasp_blas_array_dotprod(m, v2, v2);
    //alpha2 = (r,v2)
    alpha2 = fasp_blas_array_dotprod(m, r, v2);
    
    rho2 = beta - gamma*gamma/rho1;
    
    alpha3 = alpha1/rho1 - gamma*alpha2/(rho1*rho2);
    
    alpha4 = alpha2/rho2;
    
    // x = alpha3*x + alpha4*x1
    fasp_blas_array_ax(m, alpha3, x);
    fasp_blas_array_axpy(m, alpha4, x1, x);
    
    // free
    fasp_mem_free(work);
    return FASP_SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
