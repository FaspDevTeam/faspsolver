/*! \file  PreAMGAggregation.inl
 *
 *  \brief Utilities for aggregation methods
 *
 *  \note  This file contains Level-4 (Pre) functions, which are used in:
 *         PreAMGSetupSA.c, PreAMGSetupUA.c, PreAMGSetupSABSR.c, 
 *         and PreAMGSetupUABSR.c
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2012--2017 by the FASP team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *---------------------------------------------------------------------------------
 *
 *  \warning This file is also used in FASP4BLKOIL!!!
 *
 *  // TODO: Get rid of unused functions! --Chensong
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#define SYMMETRIC_PAIRWISE 1 // use symmetric pairwise aggregation

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static dCSRmat condenseBSR (const dBSRmat *A)
 *
 * \brief Form a dCSRmat matrix from a dBSRmat matrix: use the (1,1)-entry
 *
 * \param A    Pointer to the BSR format matrix
 *
 * \return     dCSRmat matrix if succeed, NULL if fail
 *
 * \author Xiaozhe Hu
 * \date   03/16/2012
 */
static dCSRmat condenseBSR (const dBSRmat *A)
{
    // information about A
    const INT   ROW = A->ROW;
    const INT   COL = A->COL;
    const INT   NNZ = A->NNZ;
    const SHORT  nc = A->nb;
    const SHORT nc2 = nc*nc;
    const REAL  TOL = 1e-8;
    
    const REAL *val = A->val;
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;
    
    // (1,1) block
    dCSRmat  P_csr = fasp_dcsr_create(ROW, COL, NNZ);
    REAL    *Pval  = P_csr.val;
    memcpy (P_csr.JA, JA, NNZ*sizeof(INT));
    memcpy (P_csr.IA, IA, (ROW+1)*sizeof(INT));
    
#ifdef _OPENMP
    INT i;
    
#pragma omp parallel for if(NNZ>OPENMP_HOLDS)
    for ( i=NNZ-1; i>=0; i-- ) Pval[i] = val[i*nc2];

#else
    INT i, j;
    
    for ( i=NNZ, j=NNZ*nc2-nc2 + (0*nc+0); i--; j-=nc2 ) Pval[i] = val[j];
#endif
    
    // compress CSR format
    fasp_dcsr_compress_inplace (&P_csr,TOL);
    
    // return P
    return P_csr;
}

/**
 * \fn static dCSRmat condenseBSRLinf (const dBSRmat *A)
 *
 * \brief Form a dCSRmat matrix from a dBSRmat matrix: use inf-norm of each block
 *
 * \param A    Pointer to the BSR format matrix
 *
 * \return     dCSRmat matrix if succeed, NULL if fail
 *
 * \author Xiaozhe Hu
 * \date   05/25/2014
 */
static dCSRmat condenseBSRLinf (const dBSRmat *A)
{
    // information about A
    const INT   ROW = A->ROW;
    const INT   COL = A->COL;
    const INT   NNZ = A->NNZ;
    const SHORT  nc = A->nb;
    const SHORT nc2 = nc*nc;
    const REAL  TOL = 1e-8;
    
    const REAL *val = A->val;
    const INT  *IA  = A->IA;
    const INT  *JA  = A->JA;
    
    // CSR matrix
    dCSRmat  Acsr = fasp_dcsr_create (ROW, COL, NNZ);
    REAL    *Aval = Acsr.val;
    
    // get structure
    memcpy (Acsr.JA, JA, NNZ*sizeof(INT));
    memcpy (Acsr.IA, IA, (ROW+1)*sizeof(INT));
    
    INT i, j, k;
    INT row_start, row_end;
    
    for ( i=0; i<ROW; i++ ) {
        
        row_start = A->IA[i]; row_end = A->IA[i+1];
        
        for ( k = row_start; k < row_end; k++ ) {
            j = A->JA[k];
            Aval[k] = fasp_smat_Linf (val+k*nc2, nc);
            if ( i != j ) Aval[k] = -Aval[k];
        }
        
    }
    
    // compress CSR format
    fasp_dcsr_compress_inplace(&Acsr,TOL);
    
    // return CSR matrix
    return Acsr;
}

/**
 * \fn static void form_boolean_p (const ivector *vertices, dCSRmat *tentp, 
 *                                 const INT NumLevels, const INT NumAggregates)
 *
 * \brief Form aggregation based on strong coupled neighbors
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param NumLevels          Level number
 * \param NumAggregates      Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 */
static void form_boolean_p (const ivector  *vertices,
                            dCSRmat        *tentp,
                            const INT       NumLevels,
                            const INT       NumAggregates)
{
    INT i, j;
    
    /* Form tentative prolongation */
    tentp->row = vertices->row;
    tentp->col = NumAggregates;
    tentp->nnz = vertices->row;
    tentp->IA  = (INT *)fasp_mem_calloc(tentp->row+1,sizeof(INT));
    
    // local variables
    INT       *IA   = tentp->IA;
    INT       *vval = vertices->val;
    const INT  row  = tentp->row;
    
    // first run
    for ( i = 0, j = 0; i < row; i++ ) {
        IA[i] = j;
        if (vval[i] > UNPT) j++;
    }
    IA[row] = j;
    
    // allocate memory for P
    tentp->nnz = j;
    tentp->JA  = (INT *)fasp_mem_calloc(tentp->nnz, sizeof(INT));
    tentp->val = (REAL *)fasp_mem_calloc(tentp->nnz, sizeof(REAL));
    
    INT  *JA = tentp->JA;
    REAL *val = tentp->val;
    
    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > UNPT) {
            JA[j] = vval[i];
            val[j] = 1.0;
            j ++;
        }
    }
}

/**
 * \fn static void form_pairwise (const dCSRmat *A, const INT pair,
 *                                const REAL k_tg, ivector *vertices,
 *                                INT *NumAggregates)
 *
 * \brief Form aggregation based on pairwise matching
 *
 * \param A                 Pointer to the coefficient matrices
 * \param pair              Number of pairs in matching
 * \param vertices          Pointer to the aggregation of vertices
 * \param NumAggregates     Pointer to number of aggregations
 *
 * \author Xiaoping Li, Zheng Li, Chensong Zhang
 * \date   04/21/2014
 *
 * \note Refer to Artem Napov and Yvan Notay "An algebraic multigrid
 *       method with guaranteed convergence rate" 2011.
 */
static void form_pairwise (const dCSRmat  *A,
                           const INT       pair,
                           const REAL      k_tg,
                           ivector        *vertices,
                           INT            *NumAggregates)
{
    const INT row  = A->row;
    
    const INT  *AIA  = A->IA;
    const INT  *AJA  = A->JA;
    const REAL *Aval = A->val;
    
    INT   i, j, row_start, row_end;
    REAL  sum;
    
    INT   col, index = 0;
    REAL  mu, min_mu, aii, ajj, aij;
    REAL  temp1, temp2, temp3, temp4;
    
    /*---------------------------------------------------------*/
    /* Step 1. select extremely strong diagonal dominate rows  */
    /*         and store in G0.                                */
    /*         G0        : vertices->val[i]=G0PT               */
    /*         Remaining : vertices->val[i]=UNPT               */
    /*---------------------------------------------------------*/
    
    fasp_ivec_alloc(row, vertices);
    
    if ( pair == 1 ) {
        for ( i = 0; i < row; i++ ) {
            sum = 0.0;
            row_start = AIA[i];
            row_end = AIA[i+1];
            
            for ( j = row_start+1; j < row_end; j++) sum += ABS(Aval[j]);
            
            if ( Aval[AIA[i]] >= ((k_tg+1.)/(k_tg-1.))*sum) {
                vertices->val[i] = G0PT;
            }
            else {
                vertices->val[i] = UNPT;
            }
        }
    }
    else {
        fasp_iarray_set(row, vertices->val, UNPT);
    }
    
    /*---------------------------------------------------------*/
    /* Step 2. compute row sum (off-diagonal) for each vertex  */
    /*---------------------------------------------------------*/
    
    REAL *s = (REAL *)fasp_mem_calloc(row, sizeof(REAL));
    
    for ( i = 0; i < row; i++ ) {
        s[i] = 0.0;
        
        if ( vertices->val[i] == G0PT ) continue;
        
        row_start = AIA[i]; row_end = AIA[i+1];
        for ( j = row_start + 1; j < row_end; j++ ) s[i] -= Aval[j];
    }
    
    /*---------------------------------------------------------*/
    /* Step 3. start the pairwise aggregation                  */
    /*---------------------------------------------------------*/
    
    *NumAggregates = 0;
    
    for ( i = 0; i < row; i++ ) {
        
        if ( vertices->val[i] != UNPT ) continue;
        
        aij = 0.0;
        min_mu = BIGREAL;
        
        row_start = AIA[i]; row_end = AIA[i+1];
        
        aii = Aval[row_start];
        
        for ( j= row_start + 1; j < row_end; j++ ) {
            col = AJA[j];
            if ( vertices->val[col] != UNPT ) continue;
            
            aij = Aval[j];
            ajj = Aval[AIA[col]];
            
            temp1 = aii+s[i]+2*aij;
            temp2 = ajj+s[col]+2*aij;
            temp2 = 1.0/temp1+1.0/temp2;
            
            temp3 = MAX(ABS(aii-s[i]), SMALLREAL); // avoid temp3 to be zero
            temp4 = MAX(ABS(ajj-s[col]), SMALLREAL); // avoid temp4 to be zero
            temp4 = -aij+1./(1.0/temp3+1.0/temp4);
            
            mu    = (-aij+1.0/temp2) / temp4;
            
            if ( min_mu > mu ) {
                min_mu = mu;
                index  = col;
            }
        }
        
        vertices->val[i] = *NumAggregates;
        
        if ( min_mu <= k_tg ) vertices->val[index] = *NumAggregates;
        
        *NumAggregates += 1;
    }
    
    fasp_mem_free(s);
}

/**
 * \fn static SHORT aggregation_pairwise (dCSRmat *A, AMG_param *param,
 *                                        const INT level, ivector *vertices,
 *                                        INT *NumAggregates)
 *
 * \brief AMG coarsening based on pairwise matching aggregation
 *
 * \param mgl               Pointer to multigrid data
 * \param param             Pointer to AMG parameters
 * \param level             Level number
 * \param vertices          Pointer to the aggregation of vertices
 * \param NumAggregates     Pointer to number of aggregations
 *
 * \author Xiaoping Li, Zheng Li, Chensong Zhang
 * \date   04/21/2014
 *
 * \note Setup A, P, PT and levels using the pairwise aggregation;
 *       Refer to A. Napov and Y. Notay
 *       "An algebraic multigrid method with guaranteed convergence rate", 2012
 *
 * Modified by Chensong Zhang, Zheng Li on 07/29/2014
 */
static SHORT aggregation_pairwise (AMG_data   *mgl,
                                   AMG_param  *param,
                                   const INT   level,
                                   ivector    *vertices,
                                   INT        *NumAggregates)
{
    const INT  pair_number = param->pair_number;
    dCSRmat  * ptrA = &mgl[level].A;
    REAL       quality_bound = param->quality_bound;
    
    INT        i, j, k, num_agg = 0, aggindex;
    INT        lvl = level;
    REAL       isorate;
    
    SHORT      dopass = 0, domin = 0;
    SHORT      status = FASP_SUCCESS;
    
#if SYMMETRIC_PAIRWISE == 0
    ivector  map1, map2;
    INT  *order;
    REAL *s = (REAL*)fasp_mem_calloc(ptrA->row, sizeof(REAL));
#endif
    
#if SYMMETRIC_PAIRWISE == 1
    INT bandwidth = fasp_dcsr_bandwidth(&mgl[level].A);
    if (bandwidth > 5.0)
        param->quality_bound = quality_bound = 1.0*bandwidth;
#endif
    
    for ( i = 1; i <= pair_number; ++i ) {
        
#if SYMMETRIC_PAIRWISE == 1
        /*-- generate aggregations by pairwise matching --*/
        form_pairwise(ptrA, i, quality_bound, &vertices[lvl], &num_agg);
#else
        if ( i == 1 ) {
            first_pairwise_unsymm(ptrA, quality_bound, order, &vertices[lvl], &map1, s, &num_agg);
        }
        else {
            second_pairwise_unsymm(&mgl[level].A, ptrA, quality_bound, i, order, &map1, &vertices[lvl-1],
                                   &vertices[lvl], &map2, s, &num_agg);
        }
#endif
        
        /*-- check number of aggregates in the first pass --*/
        if ( i == 1 && num_agg < MIN_CDOF ) {
            for ( domin=k=0; k<ptrA->row; k++ ) {
                if ( vertices[lvl].val[k] == G0PT ) domin ++;
            }
            isorate = (REAL)num_agg/domin;
            if ( isorate < 0.1 ) {
                status = ERROR_AMG_COARSEING; goto END;
            }
        }
        
        if ( i < pair_number ) {
            
            /*-- Form Prolongation --*/
            form_boolean_p(&vertices[lvl], &mgl[lvl].P, lvl+1, num_agg);
            
            /*-- Perform aggressive coarsening only up to the specified level --*/
            if ( mgl[lvl].P.col < MIN_CDOF ) break;
            
            /*-- Form restriction --*/
            fasp_dcsr_trans(&mgl[lvl].P, &mgl[lvl].R);
            
            /*-- Form coarse level stiffness matrix --*/
            fasp_blas_dcsr_rap_agg(&mgl[lvl].R, ptrA, &mgl[lvl].P, &mgl[lvl+1].A);
            
            ptrA = &mgl[lvl+1].A;
            
            fasp_dcsr_free(&mgl[lvl].P);
            fasp_dcsr_free(&mgl[lvl].R);
        }
        lvl ++; dopass ++;
    }
    
    // Form global aggregation indices
    if ( dopass > 1 ) {
        for ( i = 0; i < mgl[level].A.row; ++i ) {
            aggindex = vertices[level].val[i];
            if ( aggindex < 0 ) continue;
            for ( j = 1; j < dopass; ++j ) aggindex = vertices[level+j].val[aggindex];
            vertices[level].val[i] = aggindex;
        }
    }
    *NumAggregates = num_agg;
    
    /*-- clean memory --*/
    for ( i = 1; i < dopass; ++i ) {
        fasp_dcsr_free(&mgl[level+i].A);
        fasp_ivec_free(&vertices[level+i]);
    }
    
#if SYMMETRIC_PAIRWISE == 0
    fasp_ivec_free(&map1);
    fasp_ivec_free(&map2);
    fasp_mem_free(s);
#endif
    
END:
    return status;
}

/**
 * \fn static SHORT aggregation_vmb (dCSRmat *A, ivector *vertices, AMG_param *param,
 *                                   const INT NumLevels, dCSRmat *Neigh, 
 *                                   INT *NumAggregates)
 *
 * \brief Form aggregation based on strong coupled neighbors
 *
 * \param A                 Pointer to the coefficient matrices
 * \param vertices          Pointer to the aggregation of vertices
 * \param param             Pointer to AMG parameters
 * \param NumLevels          Level number
 * \param Neigh             Pointer to strongly coupled neighbors
 * \param NumAggregates  Pointer to number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * \note Setup A, P, PT and levels using the unsmoothed aggregation algorithm;
 *       Refer to P. Vanek, J. Madel and M. Brezina
 *       "Algebraic Multigrid on Unstructured Meshes", 1994
 *
 * Modified by Zheng Li, Chensong Zhang on 07/29/2014
 */
static SHORT aggregation_vmb (dCSRmat    *A,
                              ivector    *vertices,
                              AMG_param  *param,
                              const INT   NumLevels,
                              dCSRmat    *Neigh,
                              INT        *NumAggregates)
{
    const INT    row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
    const INT  * AIA = A->IA, * AJA = A->JA;
    const REAL * Aval = A->val;
    const INT    max_aggregation = param->max_aggregation;
    
    // return status
    SHORT  status = FASP_SUCCESS;
    
    // local variables
    INT    num_left = row;
    INT    subset, count;
    INT  * num_each_agg;
    
    REAL   strongly_coupled, strongly_coupled2;
    INT    i, j, index, row_start, row_end;
    INT  * NIA = NULL, * NJA = NULL;
    REAL * Nval = NULL;
    
    dvector diag;
    fasp_dcsr_getdiag(0, A, &diag);  // get the diagonal entries
    
    if ( GE(param->tentative_smooth, SMALLREAL) ) {
        strongly_coupled = param->strong_coupled * pow(0.5, NumLevels-1);
    }
    else {
        strongly_coupled = param->strong_coupled;
    }
    strongly_coupled2 = pow(strongly_coupled,2);
    
    /*------------------------------------------*/
    /*    Form strongly coupled neighborhood    */
    /*------------------------------------------*/
    fasp_dcsr_alloc(row, col, nnz, Neigh);
    
    NIA  = Neigh->IA;
    NJA  = Neigh->JA;
    Nval = Neigh->val;
    
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
    for ( i = row; i >= 0; i-- ) NIA[i] = AIA[i];
    
    for ( index = i = 0; i < row; ++i ) {
        NIA[i] = index;
        row_start = AIA[i]; row_end = AIA[i+1];
        for ( j = row_start; j < row_end; ++j ) {
            if ( (AJA[j] == i)
                || (pow(Aval[j],2) >= strongly_coupled2*ABS(diag.val[i]*diag.val[AJA[j]]))) {
                NJA[index] = AJA[j];
                Nval[index] = Aval[j];
                index++;
            }
        }
    }
    NIA[row] = index;
    
    Neigh->nnz = index;
    Neigh->JA  = (INT*) fasp_mem_realloc(Neigh->JA,  (Neigh->IA[row])*sizeof(INT));
    Neigh->val = (REAL*)fasp_mem_realloc(Neigh->val, (Neigh->IA[row])*sizeof(REAL));
    
    NIA  = Neigh->IA;
    NJA  = Neigh->JA;
    Nval = Neigh->val;
    
    fasp_dvec_free(&diag);
    
    /*------------------------------------------*/
    /*             Initialization               */
    /*------------------------------------------*/
    fasp_ivec_alloc(row, vertices);
    fasp_iarray_set(row, vertices->val, -2);
    *NumAggregates = 0;
    
    /*-------------*/
    /*   Step 1.   */
    /*-------------*/
    for ( i = 0; i < row; ++i ) {
        if ( (AIA[i+1] - AIA[i]) == 1 ) {
            vertices->val[i] = UNPT;
            num_left--;
        }
        else {
            subset = TRUE;
            row_start = NIA[i]; row_end = NIA[i+1];
            for ( j = row_start; j < row_end; ++j ) {
                if ( vertices->val[NJA[j]] >= UNPT ) {
                    subset = FALSE;
                    break;
                }
            }
            if ( subset ) {
                count = 0;
                vertices->val[i] = *NumAggregates;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for ( j = row_start; j < row_end; ++j ) {
                    if ( (NJA[j]!=i) && (count < max_aggregation) ) {
                        vertices->val[NJA[j]] = *NumAggregates;
                        num_left--;
                        count ++;
                    }
                }
                (*NumAggregates)++;
            }
        }
    }
    
    /*-------------*/
    /*   Step 2.   */
    /*-------------*/
    INT *temp_C = (INT*)fasp_mem_calloc(row,sizeof(INT));
    
    if ( *NumAggregates < MIN_CDOF ) {
        status = ERROR_AMG_COARSEING; goto END;
    }
    
    num_each_agg = (INT*)fasp_mem_calloc(*NumAggregates,sizeof(INT));
    
    //for ( i = 0; i < *NumAggregates; i++ ) num_each_agg[i] = 0; // initialize
    
    for ( i = row; i--; ) {
        temp_C[i] = vertices->val[i];
        if ( vertices->val[i] >= 0 ) num_each_agg[vertices->val[i]] ++;
    }
    
    for ( i = 0; i < row; ++i ) {
        if ( vertices->val[i] < UNPT ) {
            row_start = NIA[i]; row_end = NIA[i+1];
            for ( j = row_start; j < row_end; ++j ) {
                if ( temp_C[NJA[j]] > UNPT
                    && num_each_agg[temp_C[NJA[j]]] < max_aggregation ) {
                    vertices->val[i] = temp_C[NJA[j]];
                    num_left--;
                    num_each_agg[temp_C[NJA[j]]] ++ ;
                    break;
                }
            }
        }
    }
    
    /*-------------*/
    /*   Step 3.   */
    /*-------------*/
    while ( num_left > 0 ) {
        for ( i = 0; i < row; ++i ) {
            if ( vertices->val[i] < UNPT ) {
                count = 0;
                vertices->val[i] = *NumAggregates;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for ( j = row_start; j < row_end; ++j ) {
                    if ( (NJA[j]!=i) && (vertices->val[NJA[j]] < UNPT)
                        && (count<max_aggregation) ) {
                        vertices->val[NJA[j]] = *NumAggregates;
                        num_left--;
                        count++;
                    }
                }
                (*NumAggregates)++;
            }
        }
    }
    
    fasp_mem_free(num_each_agg);
    
END:
    fasp_mem_free(temp_C);
    
    return status;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
