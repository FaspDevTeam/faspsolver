/*! \file amg_setup_rs_omp.c
 *  \brief Ruge-Stuben AMG: SETUP phase.
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
/*---------------------------------omp----------------------------------------*/

/**
 * \fn int fasp_amg_setup_rs2_omp (AMG_data *mgl, AMG_param *param, 
 *                          int nthreads, int openmp_holds)
 * \brief Setup phase of Ruge and Stuben's classic AMG
 * \param mgl    pointer to AMG_data data
 * \param param  pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note Setup A, P, R, levels using classic AMG!
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller. 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_setup_rs2_omp (AMG_data *mgl, 
													 AMG_param *param, 
													 int nthreads, 
													 int openmp_holds)
{
	int status=SUCCESS;
	
#if FASP_USE_OPENMP
	const int print_level=param->print_level;
	const int m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;	
	
	int max_levels=param->max_levels;
	int mm, level=0;
	int size; // Zhiyang Zhou 2010/11/12	
	ivector vertices=fasp_ivec_create(m); // stores level info (fine: 0; coarse: 1)
	
#if DEBUG_MODE
	printf("fasp_amg_setup_rs ...... [Start]\n");
	printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8) printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	double setup_start=omp_get_wtime();
	double setup_start_ysk = omp_get_wtime();
	double setup_end_ysk = setup_start_ysk;
	//double ol_setup_start_omp, ol_setup_end_omp;
	
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
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
		//ol_setup_start_omp = omp_get_wtime();
#if DEBUG_MODE
		printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
#endif
		
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) {
			status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
			if (status < 0) goto FINISHED;
		}
		
		/*-- Coarseing and form the structure of interpolation --*/
		status = fasp_amg_coarsening_rs_omp(&mgl[level].A, &vertices, &mgl[level].P, param, nthreads,openmp_holds);
		
		/*-- Store the C/F marker: Zhiyang Zhou 2010/11/12 --*/
		size = mgl[level].A.row;
		mgl[level].cfmark = fasp_ivec_create(size);
		fasp_iarray_cp_omp(size, vertices.val, mgl[level].cfmark.val, nthreads,openmp_holds);
		
		if (mgl[level].P.col == 0) break;
		if (status < 0) goto FINISHED;
		
		/*-- Form interpolation --*/
		status = fasp_amg_interp_omp(&mgl[level].A, &vertices, &mgl[level].P, param,nthreads,openmp_holds);
		if (status < 0) goto FINISHED;
		
		/*-- Form coarse level stiffness matrix --*/
		fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
		/*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/	
#if 0
		fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);  //commented by Feng Chunsheng
#else
		fasp_blas_dcsr_rap1_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, nthreads, openmp_holds);
#endif
		// fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
		//ol_setup_end_omp = omp_get_wtime();
		//printf(" level = %d, Used Time = %lf\n", level, ol_setup_end_omp-ol_setup_start_omp);
		++level;
	}
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);
	
	for (level=1; level<max_levels; ++level) {
		mm = mgl[level].A.row;
		mgl[level].num_levels = max_levels; 		
		mgl[level].b = fasp_dvec_create(mm);
		mgl[level].x = fasp_dvec_create(mm);
		mgl[level].w = fasp_dvec_create(mm);	
	}
	
#if With_UMFPACK	
	// Need to sort the matrix A for UMFPACK to work
	dCSRmat Ac_tran;
	fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
	fasp_dcsr_sort(&Ac_tran);
	fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
	fasp_dcsr_free(&Ac_tran);
#endif
	
	if (print_level>2) {
		double gridcom=0.0, opcom=0.0;
		
		printf("-----------------------------------------------\n");
		printf("  Level     Num of rows     Num of nonzeros\n");
		printf("-----------------------------------------------\n");
		for (level=0;level<max_levels;++level) {
			printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
			gridcom += mgl[level].A.row;
			opcom += mgl[level].A.nnz;
		}
		printf("-----------------------------------------------\n");
		
		gridcom /= mgl[0].A.row;
		opcom /= mgl[0].A.nnz;
		printf("Ruge-Stuben AMG grid complexity = %f\n", gridcom);
		printf("Ruge-Stuben AMG operator complexity = %f\n", opcom);
	}
	setup_end_ysk = omp_get_wtime();
	total_setup_time += setup_end_ysk - setup_start_ysk;
	printf(" Total Setup Time is %lf\n", total_setup_time);
	
	if (print_level>0) {
		double setup_end=omp_get_wtime();
		double setupduration = setup_end - setup_start;
		print_cputime("Ruge-Stuben AMG setup",setupduration);
	}
	
	status = SUCCESS;
	
FINISHED:
	fasp_ivec_free(&vertices);
#endif
	return status;
}

/**
 * \fn int fasp_amg_setup_rs1_omp (AMG_data *mgl, AMG_param *param, 
 *                          int nthreads, int openmp_holds)
 * \brief Setup phase of Ruge and Stuben's classic AMG
 * \param mgl    pointer to AMG_data data
 * \param param  pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note Setup A, P, R, levels using classic AMG!
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller. 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_setup_rs1_omp (AMG_data *mgl, 
													 AMG_param *param, 
													 int nthreads, 
													 int openmp_holds)
{
	int status=SUCCESS;
	
#if FASP_USE_OPENMP
	const int print_level=param->print_level;
	const int m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;	
	
	int max_levels=param->max_levels;
	int mm, level=0;
	int size=0; // Zhiyang Zhou 2010/11/12	
	ivector vertices=fasp_ivec_create(m); // stores level info (fine: 0; coarse: 1)
	int *icor_ysk = (int *)fasp_mem_calloc(5*nthreads+2, sizeof(int));
	
#if DEBUG_MODE
	printf("fasp_amg_setup_rs ...... [Start]\n");
	printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8) printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	double setup_start=omp_get_wtime();
	double ol_setup_start_omp, ol_setup_end_omp;
	
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
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
		ol_setup_start_omp = omp_get_wtime();
#if DEBUG_MODE
		printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
#endif
		
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) {
			status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
			if (status < 0) goto FINISHED;
		}
		
		/*-- Coarseing and form the structure of interpolation --*/
		status = fasp_amg_coarsening_rs_omp(&mgl[level].A, &vertices, &mgl[level].P, param, nthreads,openmp_holds);
		
		/*-- Store the C/F marker: Zhiyang Zhou 2010/11/12 --*/
		size = mgl[level].A.row;
		mgl[level].cfmark = fasp_ivec_create(size);
		fasp_iarray_cp_omp(size, vertices.val, mgl[level].cfmark.val, nthreads,openmp_holds);
		
		if (mgl[level].P.col == 0) break;
		if (status < 0) goto FINISHED;
		
		/*-- Form interpolation --*/
		status = fasp_amg_interp1_omp(&mgl[level].A, &vertices, &mgl[level].P, param, icor_ysk, nthreads, openmp_holds);
		if (status < 0) goto FINISHED;
		
		/*-- Form coarse level stiffness matrix --*/
		fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
		/*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/
#if 0
		fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);  //commented by Feng Chunsheng
#else
		fasp_blas_dcsr_rap3_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, icor_ysk, nthreads, openmp_holds);
#endif
		// fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
		ol_setup_end_omp = omp_get_wtime();
		//printf(" level = %d, Used Time = %lf\n", level, ol_setup_end_omp-ol_setup_start_omp);
		++level;
	}
	fasp_mem_free(icor_ysk);
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);
	
	for (level=1; level<max_levels; ++level) {
		mm = mgl[level].A.row;
		mgl[level].num_levels = max_levels;
		mgl[level].b = fasp_dvec_create(mm);
		mgl[level].x = fasp_dvec_create(mm);
		mgl[level].w = fasp_dvec_create(mm);
	}
	
#if With_UMFPACK	
	// Need to sort the matrix A for UMFPACK to work
	dCSRmat Ac_tran;
	fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
	fasp_dcsr_sort(&Ac_tran);
	fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
	fasp_dcsr_free(&Ac_tran);
#endif
	
	if (print_level>2) {
		double gridcom=0.0, opcom=0.0;
		
		printf("-----------------------------------------------\n");
		printf("  Level     Num of rows     Num of nonzeros\n");
		printf("-----------------------------------------------\n");
		for (level=0;level<max_levels;++level) {
			printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
			gridcom += mgl[level].A.row;
			opcom += mgl[level].A.nnz;
		}
		printf("-----------------------------------------------\n");
		
		gridcom /= mgl[0].A.row;
		opcom /= mgl[0].A.nnz;
		printf("Ruge-Stuben AMG grid complexity = %f\n", gridcom);
		printf("Ruge-Stuben AMG operator complexity = %f\n", opcom);
	}
	
	if (print_level>0) {
		double setup_end=omp_get_wtime();
		double setupduration = setup_end - setup_start;
		printf("Ruge-Stuben AMG setup costs %f seconds.\n", setupduration);	
	}
	
	status = SUCCESS;
	
FINISHED:	
	fasp_ivec_free(&vertices);
#endif
	return status;
}

/**
 * \fn int fasp_amg_setup_rs_omp (AMG_data *mgl, AMG_param *param, 
 *                          int nthreads, int openmp_holds)
 * \brief Setup phase of Ruge and Stuben's classic AMG
 * \param mgl    pointer to AMG_data data
 * \param param  pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note Setup A, P, R, levels using classic AMG!
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller. 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012   Modified by Feng Chunsheng
 */
int fasp_amg_setup_rs_omp (AMG_data *mgl, 
													 AMG_param *param, 
													 int nthreads, 
													 int openmp_holds)
{
	int status=SUCCESS;
	
#if FASP_USE_OPENMP
	const int print_level=param->print_level;
	const int m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;	
	
	int max_levels=param->max_levels;
	int mm, level=0;
	int size; // Zhiyang Zhou 2010/11/12	
	ivector vertices=fasp_ivec_create(m); // stores level info (fine: 0; coarse: 1)
	int *icor_ysk = (int *)fasp_mem_calloc(5*nthreads+2, sizeof(int));
	
#if DEBUG_MODE
	printf("fasp_amg_setup_rs ...... [Start]\n");
	printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8) printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	double setup_start = omp_get_wtime();
	double setup_start_ysk = omp_get_wtime();
	double setup_end_ysk = setup_start_ysk;
	//double ol_setup_start_omp, ol_setup_end_omp;
	
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
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
		//ol_setup_start_omp = omp_get_wtime();
#if DEBUG_MODE
		printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
#endif
		
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) {
			status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
			if (status < 0) goto FINISHED;
		}
		
		/*-- Coarseing and form the structure of interpolation --*/
		status = fasp_amg_coarsening_rs_omp(&mgl[level].A, &vertices, &mgl[level].P, param, nthreads,openmp_holds);
		
		/*-- Store the C/F marker: Zhiyang Zhou 2010/11/12 --*/
		size = mgl[level].A.row;
		mgl[level].cfmark = fasp_ivec_create(size);
		fasp_iarray_cp_omp(size, vertices.val, mgl[level].cfmark.val, nthreads,openmp_holds);
		
		if (mgl[level].P.col == 0) break;
		if (status < 0) goto FINISHED;
		
		/*-- Form interpolation --*/
		status = fasp_amg_interp1_omp(&mgl[level].A, &vertices, &mgl[level].P, param, icor_ysk, nthreads, openmp_holds);
//		status = fasp_amg_interp(&mgl[level].A, &vertices, &mgl[level].P, param); //, icor_ysk, nthreads, openmp_holds);
		if (status < 0) goto FINISHED;
		
		/*-- Form coarse level stiffness matrix --*/
		fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
		/*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/
#if 1 
		fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);  //commented by Feng Chunsheng
#else
		fasp_blas_dcsr_rap4_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, icor_ysk, nthreads, openmp_holds);
#endif
		// fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
		//ol_setup_end_omp = omp_get_wtime();
		//printf(" level = %d, Used Time = %lf\n", level, ol_setup_end_omp-ol_setup_start_omp);
		++level;
	}
	fasp_mem_free(icor_ysk);
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);
	
	for (level=1; level<max_levels; ++level) {
		mm = mgl[level].A.row;
		mgl[level].num_levels = max_levels;
		mgl[level].b = fasp_dvec_create(mm);
		mgl[level].x = fasp_dvec_create(mm);
		mgl[level].w = fasp_dvec_create(mm);
	}
	
#if With_UMFPACK	
	// Need to sort the matrix A for UMFPACK to work
	dCSRmat Ac_tran;
	fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
	fasp_dcsr_sort(&Ac_tran);
	fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
	fasp_dcsr_free(&Ac_tran);
#endif
	
	if (print_level>2) {
		double gridcom=0.0, opcom=0.0;
		
		printf("-----------------------------------------------\n");
		printf("  Level     Num of rows     Num of nonzeros\n");
		printf("-----------------------------------------------\n");
		for (level=0;level<max_levels;++level) {
			printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
			gridcom += mgl[level].A.row;
			opcom += mgl[level].A.nnz;
		}
		printf("-----------------------------------------------\n");
		
		gridcom /= mgl[0].A.row;
		opcom /= mgl[0].A.nnz;
		printf("Ruge-Stuben AMG grid complexity = %f\n", gridcom);
		printf("Ruge-Stuben AMG operator complexity = %f\n", opcom);
	}
	setup_end_ysk = omp_get_wtime();
	total_setup_time += setup_end_ysk - setup_start_ysk;
	printf(" Total Setup Time is %lf\n", total_setup_time);
	
	if (print_level>0) {
		double setup_end = omp_get_wtime();
		double setupduration = setup_end - setup_start;
		printf("Ruge-Stuben AMG setup costs %f seconds.\n", setupduration);	
	}
	
	status = SUCCESS;
	
FINISHED:	
	fasp_ivec_free(&vertices);
#endif
	return status;
}

/**
 * \fn int fasp_amg_setup_rs3_omp (AMG_data *mgl, AMG_param *param, 
 *                          int nthreads, int openmp_holds)
 * \brief Setup phase of Ruge and Stuben's classic AMG
 * \param mgl    pointer to AMG_data data
 * \param param  pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \note Setup A, P, R, levels using classic AMG!
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller. 
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *           Academic Press Inc., San Diego, CA, 2001. 
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_setup_rs3_omp (AMG_data *mgl, 
													 AMG_param *param, 
													 int nthreads, 
													 int openmp_holds)
{
	int status=SUCCESS;
	
#if FASP_USE_OPENMP
	const int print_level=param->print_level;
	const int m=mgl[0].A.row, n=mgl[0].A.col, nnz=mgl[0].A.nnz;	
	
	int max_levels=param->max_levels;
	int mm, level=0;
	int size; // Zhiyang Zhou 2010/11/12
	double rate_sparsity;
	ivector vertices=fasp_ivec_create(m); // stores level info (fine: 0; coarse: 1)
	int *icor_ysk = (int *)fasp_mem_calloc(5*nthreads+2, sizeof(int));
	
#if DEBUG_MODE
	printf("fasp_amg_setup_rs ...... [Start]\n");
	printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
	
	if (print_level>8) printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
	
	double setup_start = omp_get_wtime();
	double setup_start_ysk = omp_get_wtime();
	double setup_end_ysk = setup_start_ysk;
	//double ol_setup_start_omp, ol_setup_end_omp;
	
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
	
	while ((mgl[level].A.row>param->coarse_dof) && (level<max_levels-1))
	{
		//ol_setup_start_omp = omp_get_wtime();
#if DEBUG_MODE
		printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
#endif
		
		/*-- setup ILU decomposition if necessary */
		if (level<param->ILU_levels) {
			status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
			if (status < 0) goto FINISHED;
		}
		
		/*-- Coarseing and form the structure of interpolation --*/
		status = fasp_amg_coarsening_rs_omp(&mgl[level].A, &vertices, &mgl[level].P, param, nthreads,openmp_holds);
		
		/*-- Store the C/F marker: Zhiyang Zhou 2010/11/12 --*/
		size = mgl[level].A.row;
		mgl[level].cfmark = fasp_ivec_create(size);
		fasp_iarray_cp_omp(size, vertices.val, mgl[level].cfmark.val, nthreads,openmp_holds);
		
		if (mgl[level].P.col == 0) break;
		if (status < 0) goto FINISHED;
		
		/*-- Form interpolation --*/
		//Compute the rate of sparsity of the matrix
		rate_sparsity = (double)mgl[level].A.nnz / size;
		rate_sparsity /= size;
		if ((size > 50000) && (rate_sparsity < 0.5))
		{
			status = fasp_amg_interp1_omp(&mgl[level].A, &vertices, &mgl[level].P, param, icor_ysk, nthreads, openmp_holds);
			if (status < 0) goto FINISHED;
			
			/*-- Form coarse level stiffness matrix --*/
			fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
			/*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/
#if 0
			fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);  //commented by Feng Chunsheng
#else
			fasp_blas_dcsr_rap4_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, icor_ysk, nthreads, openmp_holds);
#endif
			// fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
			//ol_setup_end_omp = omp_get_wtime();
			//printf(" level = %d, Used Time = %lf\n", level, ol_setup_end_omp-ol_setup_start_omp);
		}
		else
		{
			/*-- Form interpolation --*/
			status = fasp_amg_interp2_omp(&mgl[level].A, &vertices, &mgl[level].P, param, nthreads, openmp_holds);
			if (status < 0) goto FINISHED;
			
			/*-- Form coarse level stiffness matrix --*/
			fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
			/*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/	
#if 0
			fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);  //commented by Feng Chunsheng
#else
			fasp_blas_dcsr_rap1_omp(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, nthreads, openmp_holds);
#endif
			// fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
			//ol_setup_end_omp = omp_get_wtime();
			//printf(" level = %d, Used Time = %lf\n", level, ol_setup_end_omp-ol_setup_start_omp);
		}
		++level;
	}
	fasp_mem_free(icor_ysk);
	
	// setup total level number and current level
	mgl[0].num_levels = max_levels = level+1;
	mgl[0].w = fasp_dvec_create(m);
	
	for (level=1; level<max_levels; ++level) {
		mm = mgl[level].A.row;
		mgl[level].num_levels = max_levels;
		mgl[level].b = fasp_dvec_create(mm);
		mgl[level].x = fasp_dvec_create(mm);
		mgl[level].w = fasp_dvec_create(mm);
	}
	
#if With_UMFPACK	
	// Need to sort the matrix A for UMFPACK to work
	dCSRmat Ac_tran;
	fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
	fasp_dcsr_sort(&Ac_tran);
	fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
	fasp_dcsr_free(&Ac_tran);
#endif
	
	if (print_level>2) {
		double gridcom=0.0, opcom=0.0;
		
		printf("-----------------------------------------------\n");
		printf("  Level     Num of rows     Num of nonzeros\n");
		printf("-----------------------------------------------\n");
		for (level=0;level<max_levels;++level) {
			printf("%5d  %14d  %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
			gridcom += mgl[level].A.row;
			opcom += mgl[level].A.nnz;
		}
		printf("-----------------------------------------------\n");
		
		gridcom /= mgl[0].A.row;
		opcom /= mgl[0].A.nnz;
		printf("Ruge-Stuben AMG grid complexity = %f\n", gridcom);
		printf("Ruge-Stuben AMG operator complexity = %f\n", opcom);
	}
	setup_end_ysk = omp_get_wtime();
	total_setup_time += setup_end_ysk - setup_start_ysk;
	printf(" Total Setup Time is %lf\n", total_setup_time);
	
	if (print_level>0) {
		double setup_end = omp_get_wtime();
		double setupduration = setup_end - setup_start;
		printf("Ruge-Stuben AMG setup costs %f seconds.\n", setupduration);	
	}
	
	status = SUCCESS;
	
FINISHED:	
	fasp_ivec_free(&vertices);
#endif
	return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
