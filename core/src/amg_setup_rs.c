/*! \file amg_setup_rs.c
 *  \brief Ruge-Stuben AMG: SETUP phase.
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_setup_rs(AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of Ruge and Stuben's classic AMG
 *
 * \param *mgl    pointer to AMG_data data
 * \param *param  pointer to AMG parameters
 *
 * \note Setup A, P, R, levels using classic AMG!
 *       Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller. 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Chensong Zhang
 * \date 05/09/2010 
 *
 *  Modified by Chensong Zhang on 04/04/2009.
 *  Modified by Chensong Zhang on 04/06/2010.
 *  Modified by Chensong Zhang on 05/09/2010. 
 *  Modified by Zhiyang Zhou on 11/17/2010.
 *  Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle.
 *  Modified by Chensong zhang on 09/09/2011: add min dof.
 */
INT fasp_amg_setup_rs (AMG_data *mgl, 
                       AMG_param *param)
{
	const INT print_level=param->print_level;
	const INT m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;	
    const INT cycle_type = param->cycle_type;
	
    // local variables
	INT     mm, size;
    SHORT   level=0, status=SUCCESS;
	SHORT   max_levels=param->max_levels;
    	
	clock_t setup_start=clock();

    ivector vertices=fasp_ivec_create(m); // stores level info (fine: 0; coarse: 1)

#if DEBUG_MODE
	printf("### DEBUG: fasp_amg_setup_rs ...... [Start]\n");
	printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    	
	param->tentative_smooth = 1.0; 
    // Xiaozhe 02/23/2011: make sure classical AMG will not call fasp_blas_dcsr_mxv_agg
	
    // setup AMLI coefficients
	if (cycle_type == AMLI_CYCLE) {
		param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
		REAL lambda_max = 2.0;
		REAL lambda_min = lambda_max/4;
		fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
	}
	
	// initialize ILU parameters
	mgl->ILU_levels = param->ILU_levels;
	ILU_param iluparam;
	if (param->ILU_levels>0) {
		iluparam.print_level = param->print_level;
		iluparam.ILU_lfil    = param->ILU_lfil;
		iluparam.ILU_droptol = param->ILU_droptol;
		iluparam.ILU_relax   = param->ILU_relax;
		iluparam.ILU_type    = param->ILU_type;
	}

#if DIAGONAL_PREF
    fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
	
    // main AMG setup loop
	while ((mgl[level].A.row>MAX(50,param->coarse_dof)) && (level<max_levels-1))
	{
#if DEBUG_MODE
		printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",
               level,mgl[level].A.row,mgl[level].A.nnz);
#endif
		
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) {
			status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
			if (status < 0) goto FINISHED;
		}
		
		/*-- Coarseing and form the structure of interpolation --*/
		status = fasp_amg_coarsening_rs(&mgl[level].A, &vertices, &mgl[level].P, param);
		
		/*-- Store the C/F marker: Zhiyang Zhou 2010/11/12 --*/
		size = mgl[level].A.row;
		mgl[level].cfmark = fasp_ivec_create(size);
		memcpy(mgl[level].cfmark.val, vertices.val, size*sizeof(int));
		
		if (mgl[level].P.col == 0) break;
		if (status < 0) goto FINISHED;
		
		/*-- Form interpolation --*/
		status = fasp_amg_interp(&mgl[level].A, &vertices, &mgl[level].P, param);
		if (status < 0) goto FINISHED;
        
		/*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/	
		fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
		
        // fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
		
		// TODO: Make a new ptap using (A,P) only. R is not needed as an input! --Chensong

        ++level;

#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif                
	}
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);	
	
	for (level=1; level<max_levels; ++level) {
		mm = mgl[level].A.row;
		mgl[level].num_levels = max_levels; 		
		mgl[level].b = fasp_dvec_create(mm);
		mgl[level].x = fasp_dvec_create(mm);
		
        if (cycle_type == NL_AMLI_CYCLE)  mgl[level].w = fasp_dvec_create(3*mm);	
        else mgl[level].w = fasp_dvec_create(2*mm);
	}
	
#if With_UMFPACK	
	// Need to sort the matrix A for UMFPACK to work
	dCSRmat Ac_tran;
	fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
	fasp_dcsr_sort(&Ac_tran);
	fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
	fasp_dcsr_free(&Ac_tran);
#endif
	
	if (print_level>PRINT_NONE) {
		clock_t setup_end=clock();
		REAL setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);

        print_amgcomplexity(mgl,print_level);
		printf("Ruge-Stuben AMG setup costs %f seconds.\n\n", setupduration);	
	}
	
	status = SUCCESS;
	
FINISHED:	
	fasp_ivec_free(&vertices);

#if DEBUG_MODE
	printf("### DEBUG: fasp_amg_setup_rs ...... [Finish]\n");
#endif

	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
