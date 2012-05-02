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
 * \fn SHORT fasp_amg_setup_cr (AMG_data *mgl, AMG_param *param)
 *
 * \brief Set up phase of Brannick Falgout CR coarsening for classic AMG
 *
 * \param mgl     Pointer to AMG_data data
 * \param param   Pointer to AMG parameters
 *
 * \author James Brannick
 * \date   04/21/2010
 *
 * \note Setup A, P, R, levels using CR coarsening for 
 *       classic AMG interpolation
 *       Concrete algorithm see Brannick and Falgout 
 *          "Compatible relaxation and coarsening in AMG"  
 */
SHORT fasp_amg_setup_cr (AMG_data *mgl, 
                         AMG_param *param)
{
    dCSRmat   *A=&mgl[0].A;
    const INT  m=A->row, n=A->col, nnz=A->nnz;
    const INT  print_level=param->print_level;
    
    INT     i_0=0,i_n;
    SHORT   level=0, status=SUCCESS;
    SHORT   max_levels=param->max_levels;
    
    clock_t setup_start=clock();
    
    // The variable vertices stores level info (fine: 0; coarse: 1)
    ivector vertices=fasp_ivec_create(m); // add by Fengchunsheng /Mar/10/2011
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_cr ...... [Start]\n");
    printf("### DEBUG: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
#endif
    
    if (print_level>=PRINT_MOST) printf("fasp_amg_setup_cr: nr=%d, nc=%d, nnz=%d\n", m, n, nnz);
    
#if DIAGONAL_PREF
    fasp_dcsr_diagpref(&mgl[0].A); // reorder each row to make diagonal appear first
#endif
    
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
        
#if DIAGONAL_PREF
        fasp_dcsr_diagpref(&mgl[level].A); // reorder each row to make diagonal appear first
#endif     
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
    
    if (print_level>PRINT_NONE) {
        clock_t setup_end=clock();
        REAL setupduration = (REAL)(setup_end - setup_start)/(REAL)(CLOCKS_PER_SEC);
        print_amgcomplexity(mgl,print_level);
        print_cputime("Compatible Relaxation AMG setup",setupduration);
    }
    
    fasp_ivec_free(&vertices);    //add by Fengchunsheng /Mar/10/2011
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_setup_cr ...... [Finish]\n");
#endif
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
