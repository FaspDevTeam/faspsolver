/*! \file mgcycle.c
 *  \brief Abstract non-recursive multigrid cycle
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "forts_ns.h"

#include "mg_util.inl"

#if FASP_GSRB
INT  nx_rb=1 ;  /**< Red Black Gs Smoother Nx */
INT  ny_rb=1 ;  /**< Red Black Gs Smoother Ny */
INT  nz_rb=1 ;  /**< Red Black Gs Smoother Nz */
INT *IMAP=NULL; /**< Red Black Gs Smoother imap */
INT  MAXIMAP=1; /**< Red Black Gs Smoother max dofs of reservoir */
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_solver_mgcycle (AMG_data *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle 
 *
 * \param mgl     Pointer to AMG_data data
 * \param param   Pointer to AMG parameters
 *
 * \author Chensong Zhang
 * \date   10/06/2010
 *
 * Modified by Chensong Zhang on 12/13/2011
 * Modified by Chensong Zhang on 02/27/2013: update direct solvers.
 */
void fasp_solver_mgcycle (AMG_data *mgl, 
                          AMG_param *param)
{    
    const SHORT  amg_type=param->AMG_type;
    const SHORT  print_level = param->print_level;
    const SHORT  smoother = param->smoother;
    const SHORT  smooth_order = param->smooth_order;
    const SHORT  cycle_type = param->cycle_type;
    const SHORT  nl = mgl[0].num_levels;
    const REAL   relax = param->relaxation;
    const SHORT  ndeg = param->polynomial_degree;
    
    // local variables
    REAL alpha = 1.0;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_mgcycle ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    if (print_level >= PRINT_MOST) 
        printf("AMG_level = %d, ILU_level = %d\n", nl, param->ILU_levels);
    
    INT num_lvl[MAX_AMG_LVL] = {0}, l = 0;
    
 ForwardSweep:
    while (l<nl-1) { 
        num_lvl[l]++;
    
        // pre smoothing
        if (l<param->ILU_levels) {
            fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
        }
        else if (l<mgl->schwarz_levels) {
			switch (mgl[l].schwarz.schwarz_type) {
				case 3:
					fbgs2ns_(&(mgl[l].schwarz.A.row),
                             mgl[l].schwarz.A.IA,
                             mgl[l].schwarz.A.JA,
                             mgl[l].schwarz.A.val,
                             mgl[l].x.val,
                             mgl[l].b.val,
                             &(mgl[l].schwarz.nblk),
                             mgl[l].schwarz.iblock,
                             mgl[l].schwarz.jblock,
                             mgl[l].schwarz.mask,
                             mgl[l].schwarz.maxa,
                             mgl[l].schwarz.au,mgl[l].schwarz.al,
                             mgl[l].schwarz.rhsloc,
                             &(mgl[l].schwarz.memt));
					bbgs2ns_(&(mgl[l].schwarz.A.row),
                             mgl[l].schwarz.A.IA,
                             mgl[l].schwarz.A.JA,
                             mgl[l].schwarz.A.val,
                             mgl[l].x.val,
                             mgl[l].b.val,
                             &(mgl[l].schwarz.nblk),
                             mgl[l].schwarz.iblock,
                             mgl[l].schwarz.jblock,
                             mgl[l].schwarz.mask,
                             mgl[l].schwarz.maxa,
                             mgl[l].schwarz.au,
                             mgl[l].schwarz.al,
                             mgl[l].schwarz.rhsloc,
                             &(mgl[l].schwarz.memt));
					break;
				default:
					fbgs2ns_(&(mgl[l].schwarz.A.row),
                             mgl[l].schwarz.A.IA,
                             mgl[l].schwarz.A.JA,
                             mgl[l].schwarz.A.val,
                             mgl[l].x.val,
                             mgl[l].b.val,
                             &(mgl[l].schwarz.nblk),
                             mgl[l].schwarz.iblock,
                             mgl[l].schwarz.jblock,
                             mgl[l].schwarz.mask,
                             mgl[l].schwarz.maxa,
                             mgl[l].schwarz.au,
                             mgl[l].schwarz.al,
                             mgl[l].schwarz.rhsloc,
                             &(mgl[l].schwarz.memt));
					break;
			}
		}

        else {
#if FASP_GSRB
	        if (( l==0 )&&(nx_rb>1))
				fasp_smoother_dcsr_gs_rb3d(&mgl[l].x, &mgl[l].A, &mgl[l].b, param->presmooth_iter,
				                           1, IMAP, MAXIMAP, nx_rb, ny_rb, nz_rb);
			else
				fasp_dcsr_presmoothing(smoother,&mgl[l].A,&mgl[l].b,&mgl[l].x,param->presmooth_iter,
									   0,mgl[l].A.row-1,1,relax,ndeg,smooth_order,mgl[l].cfmark.val);
#else
            fasp_dcsr_presmoothing(smoother,&mgl[l].A,&mgl[l].b,&mgl[l].x,param->presmooth_iter,
                                   0,mgl[l].A.row-1,1,relax,ndeg,smooth_order,mgl[l].cfmark.val);
#endif
        }
    
        // form residual r = b - A x
        fasp_array_cp(mgl[l].A.row, mgl[l].b.val, mgl[l].w.val); 
        fasp_blas_dcsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);
    
        // restriction r1 = R*r0
        switch (amg_type) {    
        case UA_AMG: 
            fasp_blas_dcsr_mxv_agg(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
            break;
        default:
            fasp_blas_dcsr_mxv(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
            break;
        }
    
        // prepare for the next level
        ++l; fasp_dvec_set(mgl[l].A.row, &mgl[l].x, 0.0);    
    
    }    
    
    // CoarseSpaceSolver:    
    {
#if   WITH_MUMPS
        /* use MUMPS direct solver on the coarsest level */
        fasp_solver_mumps(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif WITH_UMFPACK
        /* use UMFPACK direct solver on the coarsest level */
        fasp_solver_umfpack(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#else    
        /* use iterative solver on the coarest level */
        fasp_coarse_itsolver(&mgl[nl-1].A, &mgl[nl-1].b, &mgl[nl-1].x, param->tol, print_level);
#endif
    }
    
    // BackwardSweep: 
    while (l>0) { 
    
        --l;
    
        // find the optimal scaling factor alpha
        if ( param->coarse_scaling == ON ) {
            alpha = fasp_blas_array_dotprod(mgl[l+1].A.row, mgl[l+1].x.val, mgl[l+1].b.val)
                / fasp_blas_dcsr_vmv(&mgl[l+1].A, mgl[l+1].x.val, mgl[l+1].x.val);
        }
    
        // prolongation u = u + alpha*P*e1
        switch (amg_type) {
        case UA_AMG:
            fasp_blas_dcsr_aAxpy_agg(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
            break;
        default:
            fasp_blas_dcsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
            break;
        }
    
        // post-smoothing
        if (l<param->ILU_levels) {
            fasp_smoother_dcsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
        }
        else if (l<mgl->schwarz_levels) {
			switch (mgl[l].schwarz.schwarz_type) {
                case 3:
                    bbgs2ns_(&(mgl[l].schwarz.A.row),
                             mgl[l].schwarz.A.IA,
                             mgl[l].schwarz.A.JA,
                             mgl[l].schwarz.A.val,
                             mgl[l].x.val,
                             mgl[l].b.val,
                             &(mgl[l].schwarz.nblk),
                             mgl[l].schwarz.iblock,
                             mgl[l].schwarz.jblock,
                             mgl[l].schwarz.mask,
                             mgl[l].schwarz.maxa,
                             mgl[l].schwarz.au,
                             mgl[l].schwarz.al,
                             mgl[l].schwarz.rhsloc,
                             &(mgl[l].schwarz.memt));
                    fbgs2ns_(&(mgl[l].schwarz.A.row),
                             mgl[l].schwarz.A.IA,
                             mgl[l].schwarz.A.JA,
                             mgl[l].schwarz.A.val,
                             mgl[l].x.val,
                             mgl[l].b.val,
                             &(mgl[l].schwarz.nblk),
                             mgl[l].schwarz.iblock,
                             mgl[l].schwarz.jblock,
                             mgl[l].schwarz.mask,
                             mgl[l].schwarz.maxa,
                             mgl[l].schwarz.au,
                             mgl[l].schwarz.al,
                             mgl[l].schwarz.rhsloc,
                             &(mgl[l].schwarz.memt));
                    break;
                default:
                    bbgs2ns_(&(mgl[l].schwarz.A.row),
                             mgl[l].schwarz.A.IA,
                             mgl[l].schwarz.A.JA,
                             mgl[l].schwarz.A.val,
                             mgl[l].x.val,
                             mgl[l].b.val,
                             &(mgl[l].schwarz.nblk),
                             mgl[l].schwarz.iblock,
                             mgl[l].schwarz.jblock,
                             mgl[l].schwarz.mask,
                             mgl[l].schwarz.maxa,
                             mgl[l].schwarz.au,
                             mgl[l].schwarz.al,
                             mgl[l].schwarz.rhsloc,
                             &(mgl[l].schwarz.memt));
                    break;
			}
		}

        else {

#if FASP_GSRB
	        if (( l==0 )&&(nx_rb>1))
				fasp_smoother_dcsr_gs_rb3d(&mgl[l].x, &mgl[l].A, &mgl[l].b, param->presmooth_iter,
				                           -1,IMAP,MAXIMAP,nx_rb,ny_rb,nz_rb);
			else
            	fasp_dcsr_postsmoothing(smoother,&mgl[l].A,&mgl[l].b,&mgl[l].x,param->postsmooth_iter,
                                        0,mgl[l].A.row-1,-1,relax,ndeg,smooth_order,mgl[l].cfmark.val);
#else
            fasp_dcsr_postsmoothing(smoother,&mgl[l].A,&mgl[l].b,&mgl[l].x,param->postsmooth_iter,
                                    0,mgl[l].A.row-1,-1,relax,ndeg,smooth_order,mgl[l].cfmark.val);
#endif
        }
    
        if (num_lvl[l]<cycle_type) break;
        else num_lvl[l] = 0;
    }
    
    if (l>0) goto ForwardSweep;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_solver_mgcycle ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_solver_mgcycle_bsr (AMG_data_bsr *mgl, AMG_param *param)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle 
 *
 * \param mgl     Pointer to AMG_data_bsr data
 * \param param   Pointer to AMG parameters
 *
 * \author Xiaozhe Hu 
 * \date   08/07/2011
 */
void fasp_solver_mgcycle_bsr (AMG_data_bsr *mgl, 
                              AMG_param *param)
{    
    const SHORT print_level = param->print_level;
    const SHORT nl          = mgl[0].num_levels;
    const SHORT smoother    = param->smoother;
    const SHORT cycle_type  = param->cycle_type;    
    //const INT smooth_order = param->smooth_order;
    const REAL relax = param->relaxation;
    
    // local variables
    INT nu_l[MAX_AMG_LVL] = {0}, l = 0;
    REAL alpha = 1.0;
    INT i;

    if (print_level >= PRINT_MOST) 
        printf("AMG_level = %d, ILU_level = %d\n", nl, param->ILU_levels);

 ForwardSweep:
    while (l<nl-1) { 
        nu_l[l]++;
        // pre smoothing
        if (l<param->ILU_levels) {
            fasp_smoother_dbsr_ilu (&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
        }
        else {
            INT steps = param->presmooth_iter;
    
            if (steps > 0) {
                switch (smoother) {
                case SMOOTHER_JACOBI:
                    for (i=0; i<steps; i++) 
                        fasp_smoother_dbsr_jacobi (&mgl[l].A, &mgl[l].b, &mgl[l].x);
                    break;
                case SMOOTHER_GS:
                    for (i=0; i<steps; i++) 
                        fasp_smoother_dbsr_gs (&mgl[l].A, &mgl[l].b, &mgl[l].x, ASCEND, NULL);
                    break;
                case SMOOTHER_SOR:
                    for (i=0; i<steps; i++) 
                        fasp_smoother_dbsr_sor (&mgl[l].A, &mgl[l].b, &mgl[l].x, ASCEND, NULL, relax);
                    break;
                default:
                    printf("### ERROR: Wrong smoother type %d!\n", smoother);
                    exit(ERROR_INPUT_PAR);                }
            }
        }
    
        // form residual r = b - A x
        fasp_array_cp(mgl[l].A.ROW*mgl[l].A.nb, mgl[l].b.val, mgl[l].w.val); 
        fasp_blas_dbsr_aAxpy(-1.0,&mgl[l].A, mgl[l].x.val, mgl[l].w.val);
    
        // restriction r1 = R*r0
        fasp_blas_dbsr_mxv(&mgl[l].R, mgl[l].w.val, mgl[l+1].b.val);
    
        // prepare for the next level
        ++l; fasp_dvec_set(mgl[l].A.ROW*mgl[l].A.nb, &mgl[l].x, 0.0);    
    
    }    
    
    // CoarseSpaceSolver:    
    {
#if   WITH_MUMPS
        /* use MUMPS direct solver on the coarsest level */
        fasp_solver_mumps(&mgl[nl-1].Ac, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif WITH_UMFPACK
        /* use UMFPACK direct solver on the coarsest level */
        fasp_solver_umfpack(&mgl[nl-1].Ac, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#elif WITH_SuperLU
        /* use SuperLU direct solver on the coarsest level */
        fasp_solver_superlu(&mgl[nl-1].Ac, &mgl[nl-1].b, &mgl[nl-1].x, 0);
#else    
        /* use iterative solver on the coarest level */
        const INT  csize = mgl[nl-1].A.ROW*mgl[nl-1].A.nb;
        const INT  cmaxit = MIN(csize*csize, 1000); // coarse level iteration number
        const REAL ctol = param->tol; // coarse level tolerance
        if ( fasp_solver_dbsr_pvgmres(&mgl[nl-1].A,&mgl[nl-1].b,&mgl[nl-1].x,NULL,ctol,cmaxit,25,1,0)<0 ) {
            printf("### WARNING: Coarse level solver does not converge in %d iterations!", cmaxit);
        }
#endif 
    }
    
    // BackwardSweep: 
    while (l>0) { 
        --l;
    
        // prolongation u = u + alpha*P*e1
        if ( param->coarse_scaling == ON ) {
            dvector PeH, Aeh;
            PeH.row = Aeh.row = mgl[l].b.row; 
            PeH.val = mgl[l].w.val + mgl[l].b.row;
            Aeh.val = PeH.val + mgl[l].b.row;
    
            fasp_blas_dbsr_mxv (&mgl[l].P, mgl[l+1].x.val,  PeH.val);
            fasp_blas_dbsr_mxv (&mgl[l].A, PeH.val, Aeh.val);
    
            alpha = (fasp_blas_array_dotprod (mgl[l].b.row, Aeh.val, mgl[l].w.val))
                / (fasp_blas_array_dotprod (mgl[l].b.row, Aeh.val, Aeh.val));
            fasp_blas_array_axpy (mgl[l].b.row, alpha, PeH.val, mgl[l].x.val);
        }
        else {
            fasp_blas_dbsr_aAxpy(alpha, &mgl[l].P, mgl[l+1].x.val, mgl[l].x.val);
        }    
        
        // post-smoothing
        if (l<param->ILU_levels) {
            fasp_smoother_dbsr_ilu(&mgl[l].A, &mgl[l].b, &mgl[l].x, &mgl[l].LU);
        }
        else {
            INT steps = param->postsmooth_iter;
    
            if (steps > 0) {
                switch (smoother) {
                case SMOOTHER_JACOBI:
                    for (i=0; i<steps; i++) 
                        fasp_smoother_dbsr_jacobi(&mgl[l].A, &mgl[l].b, &mgl[l].x);
                    break;
                case SMOOTHER_GS:
                    for (i=0; i<steps; i++) 
                        fasp_smoother_dbsr_gs(&mgl[l].A, &mgl[l].b, &mgl[l].x, ASCEND, NULL);
                    break;
                case SMOOTHER_SOR:
                    for (i=0; i<steps; i++) 
                        fasp_smoother_dbsr_sor(&mgl[l].A, &mgl[l].b, &mgl[l].x, ASCEND, NULL, relax);
                    break;
                default:
                    printf("### ERROR: Wrong smoother type %d!\n", smoother);
                    exit(ERROR_INPUT_PAR);
                }
            }
        }
    
        if (nu_l[l]<cycle_type) break;
        else nu_l[l] = 0;
    }
    
    if (l>0) goto ForwardSweep;
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
