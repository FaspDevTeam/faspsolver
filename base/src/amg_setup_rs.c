/*! \file amg_setup_rs.c
 *  \brief Ruge-Stuben AMG: SETUP phase.
 */

#include <time.h>
#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void dCSRmat_Division_Groups(dCSRmat A, INT *result, INT *groups)
 *
 * \brief Use the algebra method to get matrix's colored classes
 *
 * \param A    Input dCSRmat
 * \param result  Pointer to Return group flags
 * \param groups  Return group numbers
 *
 * \author Chunsheng Feng
 * \date   09/15/2012
 */
void dCSRmat_Division_Groups(dCSRmat A,
                             INT *result,
                             INT *groups)
{
    
    INT k,i,j,pre,group;
    INT front,rear;
    INT n=A.row;
    
    INT *cq = (INT *)malloc(sizeof(INT)*n);
    INT *newr = (INT *)malloc(sizeof(INT)*n);
    
#ifdef _OPENMP
#pragma omp parallel for if(n>OPENMP_HOLDS)
#endif
    for(k=0;k<n;k++) cq[k]=k+1;
    front=n-1;
    rear=n-1;
    memset(newr, 0, sizeof(INT)*n);
    
    group=1;
    pre=0;
    
    do{
        front=(front+1)%n;
        i=cq[front];
        if(i<pre)
        {   group++;
            result[i-1]=group;
            for(j= A.IA[i-1]; j< A.IA[i]; j++)
                if (A.JA[j] != i-1)  newr[ A.JA[j] ] = group;
        }
        else if(newr[i-1]==group)
        {   rear=(rear+1)%n;
            cq[rear]=i;
        }
        else
        {    result[i-1]=group;
            for(j= A.IA[i-1]; j< A.IA[i]; j++)
                if (A.JA[j] != i-1)  newr[ A.JA[j] ] = group;
            
            
        }
        pre=i;
    }while(rear!=front);
    
    *groups = group;
    free(cq);
    free(newr);
}

/**
 * \fn SHORT fasp_amg_setup_rs (AMG_data *mgl, AMG_param *param)
 *
 * \brief Setup phase of Ruge and Stuben's classic AMG
 *
 * \param mgl    Pointer to AMG_data data
 * \param param  Pointer to AMG parameters
 *
 * \author Chensong Zhang
 * \date   05/09/2010
 *
 * \note Setup A, P, R, levels using classic AMG!
 *       Refter to "Multigrid"
 *       by U. Trottenberg, C. W. Oosterlee and A. Schuller.
 *       Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *       Academic Press Inc., San Diego, CA, 2001.
 *
 * Modified by Chensong Zhang on 04/04/2009.
 * Modified by Chensong Zhang on 04/06/2010.
 * Modified by Chensong Zhang on 05/09/2010.
 * Modified by Zhiyang Zhou on 11/17/2010.
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle.
 * Modified by Chensong zhang on 09/09/2011: add min dof.
 */
SHORT fasp_amg_setup_rs (AMG_data *mgl,
                         AMG_param *param)
{
    const INT print_level = param->print_level;
    const INT m = mgl[0].A.row;
    const INT cycle_type = param->cycle_type;
    
    // local variables
    INT     mm, size;
    SHORT   level = 0, status = SUCCESS;
    SHORT   max_levels = param->max_levels;
    REAL    setup_start, setup_end;
    
    fasp_gettime(&setup_start);
    
    ivector vertices = fasp_ivec_create(m); // stores level info (fine: 0; coarse: 1)
    
    iCSRmat S; // strong n-couplings
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_rs ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", mgl[0].A.row, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    // Xiaozhe 02/23/2011: make sure classical AMG will not call fasp_blas_dcsr_mxv_agg
    param->tentative_smooth = 1.0;
    
    // setup AMLI coefficients
    if ( cycle_type == AMLI_CYCLE ) {
        param->amli_coef = (REAL *)fasp_mem_calloc(param->amli_degree+1,sizeof(REAL));
        REAL lambda_max = 2.0;
        REAL lambda_min = lambda_max/4;
        fasp_amg_amli_coef(lambda_max, lambda_min, param->amli_degree, param->amli_coef);
    }
    
    // initialize ILU parameters
    mgl->ILU_levels = param->ILU_levels;
    ILU_param iluparam;
    if ( param->ILU_levels > 0 ) {
        iluparam.print_level = param->print_level;
        iluparam.ILU_lfil    = param->ILU_lfil;
        iluparam.ILU_droptol = param->ILU_droptol;
        iluparam.ILU_relax   = param->ILU_relax;
        iluparam.ILU_type    = param->ILU_type;
    }
    
    // initialize Schwarz parameters
	mgl->schwarz_levels = param->schwarz_levels;
	INT schwarz_mmsize  = param->schwarz_mmsize;
	INT schwarz_maxlvl  = param->schwarz_maxlvl;
	INT schwarz_type    = param->schwarz_type;
    
#if DIAGONAL_PREF
    fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
    // main AMG setup loop
    while ((mgl[level].A.row>MAX(50,param->coarse_dof)) && (level<max_levels-1)) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",
               level,mgl[level].A.row,mgl[level].A.nnz);
#endif
        
        /*-- setup ILU decomposition if necessary */
        if ( level < param->ILU_levels ) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
            if (status < 0) goto FINISHED;
        }
        
        /* -- setup Schwarz smoother if necessary */
        if ( level < param->schwarz_levels ) {
            mgl[level].schwarz.A=fasp_dcsr_sympat(&mgl[level].A);
            fasp_dcsr_shift (&(mgl[level].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[level].schwarz, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
        }
        
        /*-- Coarseing and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[level].A, &vertices, &mgl[level].P, &S, param);
        
        /*-- Store the C/F marker --*/
        size = mgl[level].A.row;
        mgl[level].cfmark = fasp_ivec_create(size);
        memcpy(mgl[level].cfmark.val, vertices.val, size*sizeof(INT));
        
        if ( mgl[level].P.col == 0 ) break;
        if ( status < 0 ) goto FINISHED;
        
        /*-- Form interpolation --*/
        status = fasp_amg_interp(&mgl[level].A, &vertices, &mgl[level].P, &S, param);
        
        if ( status < 0 ) goto FINISHED;
        
        /*-- Form coarse level stiffness matrix: There are two RAP routines available! --*/
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        
        fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        
        // fasp_blas_dcsr_ptap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        
        // TODO: Make a new ptap using (A,P) only. R is not needed as an input! --Chensong
        
        /*-- clean up! --*/
        fasp_mem_free(S.IA);
        fasp_mem_free(S.JA);
        
        ++level;
        
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif
        
    }
    
    // setup total level number and current level
    mgl[0].num_levels = max_levels = level+1;
    mgl[0].w = fasp_dvec_create(m);
    
#if 0
    INT groups;
    INT *results; //=(INT *)malloc(sizeof(INT)*mgl[level].A.row);
    for (level=0; level<max_levels && mgl[level].A.row>2000; level++){
        results=(INT *)malloc(sizeof(INT)*mgl[level].A.row);
        dCSRmat_Division_Groups(mgl[level].A, results, &groups);
        printf("row = %d groups = %d \n", mgl[level].A.row, groups);
        free(results);
    }
#endif
    
    for (level=1; level<max_levels; ++level) {
        mm = mgl[level].A.row;
        mgl[level].num_levels = max_levels;
        mgl[level].b = fasp_dvec_create(mm);
        mgl[level].x = fasp_dvec_create(mm);
        
        if (cycle_type == NL_AMLI_CYCLE) mgl[level].w = fasp_dvec_create(3*mm);
        else mgl[level].w = fasp_dvec_create(2*mm);
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK to work
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran, &mgl[max_levels-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
    if ( print_level > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        REAL setupduration = setup_end - setup_start;
        print_amgcomplexity(mgl, print_level);
        print_cputime("Ruge-Stuben AMG setup", setupduration);
    }
    
    status = SUCCESS;
    
FINISHED:
    fasp_ivec_free(&vertices);
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_rs ...... [Finish]\n");
#endif
    
    return status;
}


/**
 * \fn INT fasp_amg_setup_rs_omp (AMG_data *mgl, AMG_param *param)
 *
 * \brief  Setup of AMG based on R-S coarsening
 *
 * \param mgl    Pointer to AMG_data data
 * \param param  Pointer to AMG parameters
 *
 * \return       SUCCESS if successed, otherwise, erro information.
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   03/11/2011
 *
 * \note Setup A, P, R, levels using classic AMG!
 *
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller.
 *           Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben).
 *           Academic Press Inc., San Diego, CA, 2001.
 */

INT fasp_amg_setup_rs_omp (AMG_data *mgl,
                           AMG_param *param)
{
    INT status=SUCCESS;
    
#ifdef _OPENMP
    const INT print_level=param->print_level;
    const INT m=mgl[0].A.row;
    const INT cycle_type = param->cycle_type;
    const INT interp_type = param->interpolation_type;
    
	//local variables
    INT mm, size, nthreads;
    INT level=0;
    INT max_levels=param->max_levels;
    REAL setup_start, setup_end, setup_duration;
    
    // set thread number
    nthreads = FASP_GET_NUM_THREADS();
    
    // stores level info (fine: 0; coarse: 1)
    ivector vertices=fasp_ivec_create(m);
    INT *icor_ysk = (INT *)fasp_mem_calloc(5*nthreads+2, sizeof(INT));
    memset(icor_ysk, 0x0, sizeof(INT)*(5*nthreads+2));
    
    // strong n-couplings
    iCSRmat S;
    
#if DEBUG_MODE
    printf("fasp_amg_setup_rs ...... [Start]\n");
    printf("fasp_amg_setup_rs: nr=%d, nc=%d, nnz=%d\n", m, mgl[0].A.col, mgl[0].A.nnz);
#endif
    
    param->tentative_smooth = 1.0;
    
    fasp_gettime(&setup_start);
    
    //setup AMLI coefficients
    if(cycle_type == AMLI_CYCLE) {
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
    
	// initialize Schwarz parameters
	mgl->schwarz_levels = param->schwarz_levels;
	INT schwarz_mmsize  = param->schwarz_mmsize;
    INT schwarz_maxlvl  = param->schwarz_maxlvl;
	INT schwarz_type    = param->schwarz_type;
    
#if DIAGONAL_PREF
	fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
	// main AMG setup loop
    while ((mgl[level].A.row>MAX(50, param->coarse_dof)) && (level<max_levels-1)) {
        
#if DEBUG_MODE
        printf("### DEBUG: level = %5d  row = %14d  nnz = %16d\n",level,mgl[level].A.row,mgl[level].A.nnz);
#endif
        
        /*-- setup ILU decomposition if necessary --*/
        if (level<param->ILU_levels) {
            status = fasp_ilu_dcsr_setup(&mgl[level].A,&mgl[level].LU,&iluparam);
            if (status < 0) goto FINISHED;
        }
        
        /*-- setup Schwarz smoother if necessary --*/
        if (level<param->schwarz_levels){
            mgl[level].schwarz.A=fasp_dcsr_sympat(&mgl[level].A);
            fasp_dcsr_shift (&(mgl[level].schwarz.A), 1);
            fasp_schwarz_setup(&mgl[level].schwarz, schwarz_mmsize, schwarz_maxlvl, schwarz_type);
        }
        /*-- Coarseing and form the structure of interpolation --*/
        status = fasp_amg_coarsening_rs(&mgl[level].A, &vertices, &mgl[level].P, &S, param);
        
        size = mgl[level].A.row;
        mgl[level].cfmark = fasp_ivec_create(size);
        fasp_iarray_cp(size, vertices.val, mgl[level].cfmark.val);
        
        if (mgl[level].P.col == 0) break;
        if (status < 0) goto FINISHED;
        
        status = fasp_amg_interp1(&mgl[level].A, &vertices, &mgl[level].P, param, &S, icor_ysk);
        
        if (status < 0) goto FINISHED;
        
        /*-- Form coarse level stiffness matrix --*/
        fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        
        if(interp_type==INTERP_REG) {
            fasp_blas_dcsr_rap4(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A, icor_ysk);
        }
        else {
            fasp_blas_dcsr_rap(&mgl[level].R, &mgl[level].A, &mgl[level].P, &mgl[level+1].A);
        }
        /*-- clean up! --*/
        fasp_mem_free(S.IA);
        fasp_mem_free(S.JA);
        
        ++level;
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif
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
        //        mgl[level].w = fasp_dvec_create(2*mm);
        
        if (cycle_type == NL_AMLI_CYCLE)  mgl[level].w = fasp_dvec_create(3*mm);
        else mgl[level].w = fasp_dvec_create(2*mm);
    }
    
#if WITH_UMFPACK
    // Need to sort the matrix A for UMFPACK to work
    dCSRmat Ac_tran;
    fasp_dcsr_trans(&mgl[max_levels-1].A, &Ac_tran);
    fasp_dcsr_sort(&Ac_tran);
    fasp_dcsr_cp(&Ac_tran,&mgl[max_levels-1].A);
    fasp_dcsr_free(&Ac_tran);
#endif
    
    if ( print_level > PRINT_NONE ) {
        fasp_gettime(&setup_end);
        setup_duration = setup_end - setup_start;
        print_amgcomplexity(mgl,print_level);
        print_cputime("Ruge-Stuben AMG setup",setup_duration);
    }
    
    status = SUCCESS;
    
FINISHED:
    fasp_ivec_free(&vertices);
    
#if DEBUG_MODE
	printf("### DEBUG: fasp_amg_setup_rs ...... [Finish]\n");
#endif
    
#endif  // end of OMP
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
