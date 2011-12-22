/*! \file amg_setup_cr.c
 *  \brief Brannick-Falgout compatible relaxation based AMG: SETUP phase
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_setup_cr(AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of Brannick Falgout CR coarsening for classic AMG
 *
 * \param *mgl     pointer to AMG_data data
 * \param *param   pointer to AMG parameters
 *
 * \note Setup A, P, R, levels using CR coarsening for 
 *       classic AMG interpolation
 *       Concrete algorithm see Brannick and Falgout 
 *          "Compatible relaxation and coarsening in AMG"  
 *
 * \author James Brannick
 * \date 04/21/2010
 */
INT fasp_amg_setup_cr (AMG_data *mgl, 
                       AMG_param *param)
{
	dCSRmat *A=&mgl[0].A;
	const INT m=A->row, n=A->col, nnz=A->nnz;
    const INT print_level=param->print_level;
	
    INT i_0=0,i_n,max_levels=param->max_levels;
	INT status = SUCCESS;
	
	// The variable vertices stores level info (fine: 0; coarse: 1)
	ivector vertices=fasp_ivec_create(m); // add by Fengchunsheng /Mar/10/2011
	
	if (param->print_level>=PRINT_MOST) printf("amg_setup_cr: %d, %d, %d\n",m,n,nnz);
	
	clock_t setup_start=clock();
	
	unsigned INT level=0;	
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1)) {
		
		/*-- Coarsen and form the structure of interpolation --*/
		i_n = mgl[level].A.row-1;
		
		fasp_amg_coarsening_cr(i_0,i_n,&mgl[level].A, &vertices, param);
		
		/*-- Form interpolation --*/
		/* 1. SPARSITY -- Form ip and jp */ 
		/* First a symbolic one
		 then gather the list */
		/* 2. COEFFICIENTS -- Form P */ 
		// energymin(mgl[level].A, &vertices[level], mgl[level].P, param);
		// fasp_mem_free(vertices[level].val);
		
		/*-- Form coarse level stiffness matrix --*/
		// fasp_dcsr_trans(mgl[level].P, mgl[level].R);
		
		/*-- Form coarse level stiffness matrix --*/	
		//fasp_blas_dcsr_rap(mgl[level].R, mgl[level].A, mgl[level].P, mgl[level+1].A);
		
		++level;
	}
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);	
	
	for (level=1; level<max_levels; ++level) {
		INT m = mgl[level].A.row;
		mgl[level].num_levels = max_levels; 		
		mgl[level].b = fasp_dvec_create(m);
		mgl[level].x = fasp_dvec_create(m);
		mgl[level].w = fasp_dvec_create(m);	
	}
	
	if (print_level>=PRINT_SOME) {
		REAL gridcom=0.0, opcom=0.0;
		
		printf("-----------------------------------------------\n");
		printf("  Level     Num of rows      Num of nonzeros\n");
		printf("-----------------------------------------------\n");
		for (level=0;level<max_levels;++level) {
			printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
			gridcom += mgl[level].A.row;
			opcom += mgl[level].A.nnz;
		}
		printf("-----------------------------------------------\n");
		
		gridcom /= (REAL)mgl[0].A.row;
		opcom /= (REAL)mgl[0].A.nnz;
		printf("Compatible Relaxation AMG grid complexity = %f\n", gridcom);	
		printf("Compatible Relaxation AMG operator complexity = %f\n", opcom);	
	}
	
	if (print_level>PRINT_NONE) {
		clock_t setup_end=clock();
		REAL setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
		printf("Compatible Relaxation AMG setup costs %f seconds.\n", setupduration);	
	}
	
	fasp_ivec_free(&vertices);	//add by Fengchunsheng /Mar/10/2011
	
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
