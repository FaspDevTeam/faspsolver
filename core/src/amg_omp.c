/*! \file amg_omp.c 
 *  \brief AMG method as an iterative solver (main file)
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/
 
/** 
 * \fn int fasp_solver_amg_omp(dCSRmat *A, dvector *b, dvector *x, AMG_param *param,  
 *                             int nthreads, int openmp_holds) 
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG. 
 * 
 * \param *A      pointer to the coefficient matrix 
 * \param *b      pointer to the dvector of right hand side 
 * \param *x      pointer to the dvector of dofs 
 * \param *param  pointer to AMG parameters 
 * \param nthreads number of threads 
 * \param openmp_holds threshold of parallelization 
 * 
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller  
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben) 
 *           Academic Press Inc., San Diego, CA, 2001. 
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang 
 * \date 03/01/2011 
 */ 
int fasp_solver_amg_omp (dCSRmat *A,  
                         dvector *b,  
                         dvector *x,  
                         AMG_param *param,  
                         int nthreads,  
                         int openmp_holds) 
{ 
  int status = SUCCESS; 
	 
#if FASP_USE_OPENMP 
  const int nnz=A->nnz, m=A->row, n=A->col; 
  const int max_levels=param->max_levels; 
  const int print_level=param->print_level; 
  const int amg_type=param->AMG_type; 
  const int cycle_type=param->cycle_type; 
	 
  double AMG_start, AMG_end; 
  double AMG_duration; 
	 
#if DEBUG_MODE 
  printf("amg ...... [Start]\n"); 
  printf("amg: nr=%d, nc=%d, nnz=%d\n", m, n, nnz); 
#endif 
	 
  AMG_start=omp_get_wtime(); 
	 
  // initialize A, b, x for mgl[0]		 
  AMG_data *mgl=fasp_amg_data_create(max_levels); 
  mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp_omp(A,&mgl[0].A, nthreads, openmp_holds); 
  mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp_omp(b,&mgl[0].b, nthreads, openmp_holds); 
  mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp_omp(x,&mgl[0].x, nthreads, openmp_holds);  
	 
  // AMG setup	 
  switch (amg_type) { 
  case CLASSIC_AMG: 
    if ( (status=fasp_amg_setup_rs_omp(mgl, param, nthreads, openmp_holds))<0 ) goto FINISHED; 
    break; 
  default: 
    printf("Error: Wrong AMG type %d!\n", amg_type); goto FINISHED; 
  } 
	 
  // AMG solve 
  switch (cycle_type) 
  { 
    default: 
      if ( (status=fasp_amg_solve_omp(mgl, param, nthreads, openmp_holds)) < 0 ) goto FINISHED; 
      break; 
  } 
	 
  fasp_dvec_cp_omp(&mgl[0].x, x, nthreads, openmp_holds); // save solution vector 
	 
  if (print_level>0) { 
    AMG_end=omp_get_wtime();		 
    AMG_duration = AMG_end - AMG_start; 
    printf("AMG totally costs %f seconds.\n", AMG_duration); 
  }	 
	 
 FINISHED:	 
  fasp_amg_data_free(mgl);	// clean-up memory	 
  if (status<0) printf("Error: AMG solver quit unexpectedly!\n"); 
	 
#if DEBUG_MODE 
  printf("amg ...... [Finish]\n"); 
#endif 
	 
#endif 
	 
  return status; 
} 
 
/** 
 * \fn int fasp_solver_amg1_omp(dCSRmat *A, dvector *b, dvector *x, AMG_param *param,  
 *                             int nthreads, int openmp_holds) 
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG. 
 * 
 * \param *A      pointer to the coefficient matrix 
 * \param *b      pointer to the dvector of right hand side 
 * \param *x      pointer to the dvector of dofs 
 * \param *param  pointer to AMG parameters 
 * \param nthreads number of threads 
 * \param openmp_holds threshold of parallelization 
 * 
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller  
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben) 
 *           Academic Press Inc., San Diego, CA, 2001. 
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang 
 * \date 03/01/2011 
 */ 
int fasp_solver_amg1_omp (dCSRmat *A,  
                         dvector *b,  
                         dvector *x,  
                         AMG_param *param,  
                         int nthreads,  
                         int openmp_holds) 
{ 
  int status = SUCCESS; 
	 
#if FASP_USE_OPENMP 
  const int nnz=A->nnz, m=A->row, n=A->col; 
  const int max_levels=param->max_levels; 
  const int print_level=param->print_level; 
  const int amg_type=param->AMG_type; 
  const int cycle_type=param->cycle_type; 
	 
  double AMG_start, AMG_end; 
  double AMG_duration; 
	 
#if DEBUG_MODE 
  printf("amg ...... [Start]\n"); 
  printf("amg: nr=%d, nc=%d, nnz=%d\n", m, n, nnz); 
#endif 
	 
  AMG_start=omp_get_wtime(); 
	 
  // initialize A, b, x for mgl[0]		 
  AMG_data *mgl=fasp_amg_data_create(max_levels); 
  mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp_omp(A,&mgl[0].A, nthreads, openmp_holds); 
  mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp_omp(b,&mgl[0].b, nthreads, openmp_holds); 
  mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp_omp(x,&mgl[0].x, nthreads, openmp_holds);  
	 
  // AMG setup	 
  switch (amg_type) { 
  case CLASSIC_AMG: 
    if ( (status=fasp_amg_setup_rs1_omp(mgl, param, nthreads, openmp_holds))<0 ) goto FINISHED; 
    break; 
  default: 
    printf("Error: Wrong AMG type %d!\n", amg_type); goto FINISHED; 
  } 
	 
  // AMG solve 
  switch (cycle_type) 
  { 
    default: 
      if ( (status=fasp_amg_solve_omp(mgl, param, nthreads, openmp_holds)) < 0 ) goto FINISHED; 
      break; 
  } 
	 
  fasp_dvec_cp_omp(&mgl[0].x, x, nthreads, openmp_holds); // save solution vector 
	 
  if (print_level>0) { 
    AMG_end=omp_get_wtime();		 
    AMG_duration = AMG_end - AMG_start; 
    printf("AMG totally costs %f seconds.\n", AMG_duration); 
  }	 
	 
 FINISHED:	 
  fasp_amg_data_free(mgl);	// clean-up memory	 
  if (status<0) printf("Error: AMG solver quit unexpectedly!\n"); 
	 
#if DEBUG_MODE 
  printf("amg ...... [Finish]\n"); 
#endif 
	 
#endif 
	 
  return status; 
} 
 
/** 
 * \fn int fasp_solver_amg2_omp(dCSRmat *A, dvector *b, dvector *x, AMG_param *param,  
 *                             int nthreads, int openmp_holds) 
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG. 
 * 
 * \param *A      pointer to the coefficient matrix 
 * \param *b      pointer to the dvector of right hand side 
 * \param *x      pointer to the dvector of dofs 
 * \param *param  pointer to AMG parameters 
 * \param nthreads number of threads 
 * \param openmp_holds threshold of parallelization 
 * 
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller  
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben) 
 *           Academic Press Inc., San Diego, CA, 2001. 
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang 
 * \date 03/01/2011 
 */ 
int fasp_solver_amg2_omp (dCSRmat *A,  
                         dvector *b,  
                         dvector *x,  
                         AMG_param *param,  
                         int nthreads,  
                         int openmp_holds) 
{ 
  int status = SUCCESS; 
	 
#if FASP_USE_OPENMP 
  const int nnz=A->nnz, m=A->row, n=A->col; 
  const int max_levels=param->max_levels; 
  const int print_level=param->print_level; 
  const int amg_type=param->AMG_type; 
  const int cycle_type=param->cycle_type; 
	 
  double AMG_start, AMG_end; 
  double AMG_duration; 
	 
#if DEBUG_MODE 
  printf("amg ...... [Start]\n"); 
  printf("amg: nr=%d, nc=%d, nnz=%d\n", m, n, nnz); 
#endif 
	 
  AMG_start=omp_get_wtime(); 
	 
  // initialize A, b, x for mgl[0]		 
  AMG_data *mgl=fasp_amg_data_create(max_levels); 
  mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp_omp(A,&mgl[0].A, nthreads, openmp_holds); 
  mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp_omp(b,&mgl[0].b, nthreads, openmp_holds); 
  mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp_omp(x,&mgl[0].x, nthreads, openmp_holds);  
	 
  // AMG setup	 
  switch (amg_type) { 
  case CLASSIC_AMG: 
    if ( (status=fasp_amg_setup_rs2_omp(mgl, param, nthreads, openmp_holds))<0 ) goto FINISHED; 
    break; 
  default: 
    printf("Error: Wrong AMG type %d!\n", amg_type); goto FINISHED; 
  } 
	 
  // AMG solve 
  switch (cycle_type) 
  { 
    default: 
      if ( (status=fasp_amg_solve_omp(mgl, param, nthreads, openmp_holds)) < 0 ) goto FINISHED; 
      break; 
  } 
	 
  fasp_dvec_cp_omp(&mgl[0].x, x, nthreads, openmp_holds); // save solution vector 
	 
  if (print_level>0) { 
    AMG_end=omp_get_wtime();		 
    AMG_duration = AMG_end - AMG_start; 
    printf("AMG totally costs %f seconds.\n", AMG_duration); 
  }	 
	 
 FINISHED:	 
  fasp_amg_data_free(mgl);	// clean-up memory	 
  if (status<0) printf("Error: AMG solver quit unexpectedly!\n"); 
	 
#if DEBUG_MODE 
  printf("amg ...... [Finish]\n"); 
#endif 
	 
#endif 
	 
  return status; 
} 
 
/** 
 * \fn int fasp_solver_amg3_omp(dCSRmat *A, dvector *b, dvector *x, AMG_param *param,  
 *                             int nthreads, int openmp_holds) 
 * \brief Solve Ax=b by Ruge and Stuben's classic AMG. 
 * 
 * \param *A      pointer to the coefficient matrix 
 * \param *b      pointer to the dvector of right hand side 
 * \param *x      pointer to the dvector of dofs 
 * \param *param  pointer to AMG parameters 
 * \param nthreads number of threads 
 * \param openmp_holds threshold of parallelization 
 * 
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller  
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben) 
 *           Academic Press Inc., San Diego, CA, 2001. 
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang 
 * \date 03/01/2011 
 */ 
int fasp_solver_amg3_omp (dCSRmat *A,  
                         dvector *b,  
                         dvector *x,  
                         AMG_param *param,  
                         int nthreads,  
                         int openmp_holds) 
{ 
  int status = SUCCESS; 
	 
#if FASP_USE_OPENMP 
  const int nnz=A->nnz, m=A->row, n=A->col; 
  const int max_levels=param->max_levels; 
  const int print_level=param->print_level; 
  const int amg_type=param->AMG_type; 
  const int cycle_type=param->cycle_type; 
	 
  double AMG_start, AMG_end; 
  double AMG_duration; 
	 
#if DEBUG_MODE 
  printf("amg ...... [Start]\n"); 
  printf("amg: nr=%d, nc=%d, nnz=%d\n", m, n, nnz); 
#endif 
	 
  AMG_start=omp_get_wtime(); 
	 
  // initialize A, b, x for mgl[0]		 
  AMG_data *mgl=fasp_amg_data_create(max_levels); 
  mgl[0].A = fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp_omp(A,&mgl[0].A, nthreads, openmp_holds); 
  mgl[0].b = fasp_dvec_create(n);       fasp_dvec_cp_omp(b,&mgl[0].b, nthreads, openmp_holds); 
  mgl[0].x = fasp_dvec_create(n);       fasp_dvec_cp_omp(x,&mgl[0].x, nthreads, openmp_holds);  
	 
  // AMG setup	 
  switch (amg_type) { 
  case CLASSIC_AMG: 
    if ( (status=fasp_amg_setup_rs3_omp(mgl, param, nthreads, openmp_holds))<0 ) goto FINISHED; 
    break; 
  default: 
    printf("Error: Wrong AMG type %d!\n", amg_type); goto FINISHED; 
  } 
	 
  // AMG solve 
  switch (cycle_type) 
  { 
    default: 
      if ( (status=fasp_amg_solve_omp(mgl, param, nthreads, openmp_holds)) < 0 ) goto FINISHED; 
      break; 
  } 
	 
  fasp_dvec_cp_omp(&mgl[0].x, x, nthreads, openmp_holds); // save solution vector 
	 
  if (print_level>0) { 
    AMG_end=omp_get_wtime();		 
    AMG_duration = AMG_end - AMG_start; 
    printf("AMG totally costs %f seconds.\n", AMG_duration); 
  }	 
	 
 FINISHED:	 
  fasp_amg_data_free(mgl);	// clean-up memory	 
  if (status<0) printf("Error: AMG solver quit unexpectedly!\n"); 
	 
#if DEBUG_MODE 
  printf("amg ...... [Finish]\n"); 
#endif 
	 
#endif 
	 
  return status; 
} 
 
/*---------------------------------*/ 
/*--        End of File          --*/ 
/*---------------------------------*/ 
