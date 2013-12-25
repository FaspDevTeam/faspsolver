/*! \file interpolation.c
 *
 *  \brief Interpolation operators for AMG
 *
 *  \note Ref U. Trottenberg, C. W. Oosterlee, and A. Schuller
 *            Multigrid (Appendix A: An Intro to Algebraic Multigrid)
 *            Academic Press Inc., San Diego, CA, 2001
 *            With contributions by A. Brandt, P. Oswald and K. Stuben.
 */

#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"

static void interp_DIR (dCSRmat *, ivector *, dCSRmat *, AMG_param *);
static void interp_DIR1(dCSRmat *, ivector *, dCSRmat *, AMG_param *, INT *);
static void interp_STD (dCSRmat *, ivector *, dCSRmat *, iCSRmat *, AMG_param *);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_amg_interp (dCSRmat *A, ivector *vertices, dCSRmat *P, iCSRmat *S,
 *                           AMG_param *param)
 *
 * \brief Generate interpolation operator P
 *
 * \param A          Pointer to dCSRmat: the coefficient matrix (index starts from 0)
 * \param vertices   Indicator vector for the C/F splitting of the variables
 * \param P          Prolongation (input: nonzero pattern, output: prolongation)
 * \param S          Strong connection matrix
 * \param param      AMG parameters
 *
 * \author  Xuehai Huang, Chensong Zhang
 * \date    04/04/2010
 *
 * Modified by Xiaozhe Hu on 05/23/2012: add S as input
 * Modified by Chensong Zhang on 09/12/2012: clean up and debug interp_RS
 * Modified by Chensong Zhang on 05/14/2013: reconstruct the code
 */
void fasp_amg_interp (dCSRmat *A,
                      ivector *vertices,
                      dCSRmat *P,
                      iCSRmat *S,
                      AMG_param *param)
{
	const INT coarsening_type = param->coarsening_type;
    INT       interp_type     = param->interpolation_type;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_interp ...... [Start]\n");
#endif
    
    // make sure standard interpolaiton is used for aggressive coarsening
    if ( coarsening_type == COARSE_AC ) interp_type = INTERP_STD;
    
    switch ( interp_type ) {
            
        case INTERP_DIR: // Direct interpolation
            interp_DIR(A, vertices, P, param); break;
            
        case INTERP_STD: // Standard interpolation
            interp_STD(A, vertices, P, S, param); break;
        
        case INTERP_ENG: // Energy-min interpolation
            fasp_amg_interp_em(A, vertices, P, param); break;
        
        default:
            fasp_chkerr(ERROR_AMG_INTERP_TYPE, "fasp_amg_interp");
    
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_interp ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_amg_interp1 (dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param,
 *                            iCSRmat *S, INT *icor_ysk)
 *
 * \brief Generate interpolation operator P
 *
 * \param A          Pointer to dCSRmat: the coefficient matrix (index starts from 0)
 * \param vertices   Indicator vector for the C/F splitting of the variables
 * \param P          Prolongation (input: nonzero pattern, output: prolongation)
 * \param S          Strong connection matrix
 * \param param      AMG parameters
 * \param icor_ysk   Indices of coarse nodes in fine grid
 *
 * \return           SUCCESS or error message
 *
 * \author           Chunsheng Feng, Xiaoqiang Yue
 * \date             03/01/2011
 *
 * Modified by Chensong Zhang on 05/14/2013: reconstruct the code
 */
void fasp_amg_interp1 (dCSRmat *A,
                       ivector *vertices,
                       dCSRmat *P,
                       AMG_param *param,
                       iCSRmat *S,
                       INT *icor_ysk)
{
	const INT coarsening_type = param->coarsening_type;
    INT       interp_type     = param->interpolation_type;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_interp1 ...... [Start]\n");
#endif
    
    // make sure standard interpolaiton is used for aggressive coarsening
    if ( coarsening_type == COARSE_AC ) interp_type = INTERP_STD;
    
    switch ( interp_type ) {
            
        case INTERP_DIR: // Direct interpolation
            interp_DIR1(A, vertices, P, param, icor_ysk); break;
            
        case INTERP_STD: // Standard interpolation
            interp_STD(A, vertices, P, S, param); break;
            
        case INTERP_ENG: // Energy-min interpolation
            fasp_amg_interp_em(A, vertices, P, param); break;
            
        default:
            fasp_chkerr(ERROR_AMG_INTERP_TYPE, "fasp_amg_interp1");
            
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_interp1 ...... [Finish]\n");
#endif
}

/**
 * \fn void fasp_amg_interp_trunc (dCSRmat *P, AMG_param *param)
 *
 * \brief Trunction step for prolongation operators
 *
 * \param P        Prolongation (input: full, output: truncated)
 * \param param    Pointer to AMG_param: AMG parameters
 *
 * \author Chensong Zhang
 * \date   05/14/2013
 *
 * Originally by Xuehai Huang, Chensong Zhang on 01/31/2009
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012: add OMP support
 * Modified by Chensong Zhang on 05/14/2013: rewritten
 */
void fasp_amg_interp_trunc (dCSRmat *P,
                            AMG_param *param)
{
    const INT   row    = P->row;
    const INT   nnzold = P->nnz;
    const INT   prtlvl = param->print_level;
    const REAL  eps_tr = param->truncation_threshold;
    
    // local variables
    INT  num_nonzero = 0;    // number of non zeros after truncation
    REAL Min_neg, Max_pos;   // min negative and max positive entries
    REAL Fac_neg, Fac_pos;   // factors for negative and positive entries
    REAL Sum_neg, TSum_neg;  // sum and truncated sum of negative entries
    REAL Sum_pos, TSum_pos;  // sum and truncated sum of positive entries
    
    INT  index1 = 0, index2 = 0, begin, end;
    INT  i, j;

#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_interp_trunc ...... [Start]\n");
#endif

    for ( i = 0; i < row; ++i ) {
        
        begin = P->IA[i]; end = P->IA[i+1];

        P->IA[i] = num_nonzero;
        Min_neg  = Max_pos  = 0;
        Sum_neg  = Sum_pos  = 0;
        TSum_neg = TSum_pos = 0;
        
        // 1. Summations of positive and negative entries
        for ( j = begin; j < end; ++j ) {
            
            if ( P->val[j] > 0 ) {
                Sum_pos += P->val[j];
                Max_pos = MAX(Max_pos, P->val[j]);
            }

            else if ( P->val[j] < 0 ) {
                Sum_neg += P->val[j];
                Min_neg = MIN(Min_neg, P->val[j]);
            }

        }
        
        Max_pos *= eps_tr; Min_neg *= eps_tr;
        
        // 2. Set JA of truncated P
        for ( j = begin; j < end; ++j ) {
            
            if ( P->val[j] >= Max_pos ) {
                num_nonzero++;
                P->JA[index1++] = P->JA[j];
                TSum_pos += P->val[j];
            }

            else if ( P->val[j] <= Min_neg ) {
                num_nonzero++;
                P->JA[index1++] = P->JA[j];
                TSum_neg += P->val[j];
            }

        }
        
        // 3. Compute factors and set values of truncated P
        if ( TSum_pos > SMALLREAL ) {
            Fac_pos = Sum_pos / TSum_pos; // factor for positive entries
        }
        else {
            Fac_pos = 1.0;
        }
            
        if ( TSum_neg < -SMALLREAL ) {
            Fac_neg = Sum_neg / TSum_neg; // factor for negative entries
        }
        else {
            Fac_neg = 1.0;
        }

        for ( j = begin; j < end; ++j ) {
            
            if ( P->val[j] >= Max_pos )
                P->val[index2++] = P->val[j] * Fac_pos;

            else if ( P->val[j] <= Min_neg )
                P->val[index2++] = P->val[j] * Fac_neg;
        }
        
    }
    
    // resize the truncated prolongation P
    P->nnz = P->IA[row] = num_nonzero;
    P->JA  = (INT  *)fasp_mem_realloc(P->JA,  num_nonzero*sizeof(INT));
    P->val = (REAL *)fasp_mem_realloc(P->val, num_nonzero*sizeof(REAL));
    
    if ( prtlvl >= PRINT_MORE ) {
        printf("Truncate prolongation, #nz befor: %10d, after: %10d\n",
               nnzold, num_nonzero);
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_interp_trunc ...... [Finish]\n");
#endif

}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void interp_DIR (dCSRmat *A, ivector *vertices, dCSRmat *P,
 *                             AMG_param *param)
 *
 * \brief Direct interpolation
 *
 * \param A          Pointer to dCSRmat: the coefficient matrix (index starts from 0)
 * \param vertices   Indicator vector for the C/F splitting of the variables
 * \param P          Prolongation (input: nonzero pattern, output: prolongation)
 * \param param      Pointer to AMG_param: AMG parameters
 *
 * \author Xuehai Huang
 * \date   01/31/2009
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012: add OMP support
 * Modified by Chensong Zhang on 09/12/2012: compare with the old version
 * Modified by Chensong Zhang on 05/14/2013: reconstruct the code
 */
static void interp_DIR (dCSRmat *A,
                        ivector *vertices,
                        dCSRmat *P,
                        AMG_param *param )
{
    INT  row   = A->row;
    INT       *vec   = vertices->val;

    // local variables
    SHORT      IS_STRONG;   // is the variable strong coupled to i?
    INT        num_pcouple; // number of positive strong couplings
    INT        begin_row, end_row;
    INT        i, j, k, l, index = 0, idiag;
    
    // a_minus and a_plus for Neighbors and Prolongation support
    REAL       amN, amP, apN, apP;
    REAL       alpha, beta, aii;

    // indices of C-nodes
    INT      * cindex = (INT *)fasp_mem_calloc(row, sizeof(INT));

    INT        use_openmp = FALSE;
    
#ifdef _OPENMP
    INT myid, mybegin, myend, stride_i, nthreads;
    row = MIN(P->IA[P->row], row);
    if ( row > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    // Step 1. Fill in values for interpolation operator P
    if (use_openmp) {
#ifdef _OPENMP
        stride_i = row/nthreads;
#pragma omp parallel private(myid,mybegin,myend,i,begin_row,end_row,idiag,aii,amN,amP,apN,apP,num_pcouple,j,k,alpha,beta,l) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            mybegin = myid*stride_i;
            if(myid < nthreads-1) myend = mybegin+stride_i;
            else myend = row;
            for (i=mybegin; i<myend; ++i){
                begin_row=A->IA[i]; end_row=A->IA[i+1]-1;
                for(idiag=begin_row;idiag<=end_row;idiag++){
                    if (A->JA[idiag]==i) {
                        aii=A->val[idiag];
                        break;
                    }
                }
                if(vec[i]==0){  // if node i is on fine grid
                    amN=0, amP=0, apN=0, apP=0,  num_pcouple=0;
                    for(j=begin_row;j<=end_row;++j){
                        if(j==idiag) continue;
                        for(k=P->IA[i];k<P->IA[i+1];++k) {
                            if(P->JA[k]==A->JA[j]) break;
                        }
                        if(A->val[j]>0) {
                            apN+=A->val[j];
                            if(k<P->IA[i+1]) {
                                apP+=A->val[j];
                                num_pcouple++;
                            }
                        }
                        else {
                            amN+=A->val[j];
                            if(k<P->IA[i+1]) {
                                amP+=A->val[j];
                            }
                        }
                    } // j
                    
                    alpha=amN/amP;
                    if(num_pcouple>0) {
                        beta=apN/apP;
                    }
                    else {
                        beta=0;
                        aii+=apN;
                    }
                    for(j=P->IA[i];j<P->IA[i+1];++j){
                        k=P->JA[j];
                        for(l=A->IA[i];l<A->IA[i+1];l++){
                            if(A->JA[l]==k) break;
                        }
                        if(A->val[l]>0){
                            P->val[j]=-beta*A->val[l]/aii;
                        }
                        else {
                            P->val[j]=-alpha*A->val[l]/aii;
                        }
                    }
                }
                else if(vec[i]==2) // if node i is a special fine node
                {
                    
                }
                else {// if node i is on coarse grid
                    P->val[P->IA[i]]=1;
                }
            }
        }
#endif
    }
    
    else {
        
        for ( i = 0; i < row; ++i ) {
            
            begin_row = A->IA[i]; end_row = A->IA[i+1];
            
            // find diagonal entry first!!!
            for ( idiag = begin_row; idiag < end_row; idiag++ ) {
                if ( A->JA[idiag] == i ) {
                    aii = A->val[idiag]; break;
                }
            }
            
            if ( vec[i] == FGPT ) { // fine grid nodes
                
                amN = amP = apN = apP = 0.0;
                
                num_pcouple = 0;
                
                for ( j = begin_row; j < end_row; ++j ) {
                
                    if ( j == idiag ) continue; // skip diagonal
                    
                    // check a point strong-coupled to i or not
                    IS_STRONG = FALSE;
                    for ( k = P->IA[i]; k < P->IA[i+1]; ++k ) {
                        if ( P->JA[k] == A->JA[j] ) { IS_STRONG = TRUE; break; }
                    }
                    
                    if ( A->val[j] > 0 ) {
                        apN += A->val[j]; // sum up positive entries
                        if ( IS_STRONG ) { apP += A->val[j]; num_pcouple++; }
                    }
                    else {
                        amN += A->val[j]; // sum up negative entries
                        if ( IS_STRONG ) { amP += A->val[j]; }
                    }
                } // end for j
                
                // set weight factors
                alpha = amN / amP;
                if ( num_pcouple > 0 ) {
                    beta = apN / apP;
                }
                else {
                    beta = 0.0; aii += apN;
                }
                
                // keep aii inside the loop to avoid floating pt error! --Chensong
                for ( j = P->IA[i]; j < P->IA[i+1]; ++j ) {
                    k = P->JA[j];
                    for ( l = A->IA[i]; l < A->IA[i+1]; l++ ) {
                        if ( A->JA[l] == k ) break;
                    }
                    if ( A->val[l] > 0 ) {
                        P->val[j] = -beta  * A->val[l] / aii;
                    }
                    else {
                        P->val[j] = -alpha * A->val[l] / aii;
                    }
                }
                
            } // end if vec
            
            else if ( vec[i] == CGPT ) { // coarse grid nodes
                P->val[P->IA[i]] = 1.0;
            }
        }
    }
    
    // Step 2. Generate coarse level indices and set values of P.JA
    for ( index = i = 0; i < row; ++i ) {
        if ( vec[i] == CGPT ) cindex[i] = index++;
    }
    P->col = index;
    
    if (use_openmp) {
#ifdef _OPENMP
        stride_i = P->IA[P->row]/nthreads;
#pragma omp parallel private(myid,mybegin,myend,i,j) num_threads(nthreads)
        {
            myid = omp_get_thread_num();
            mybegin = myid*stride_i;
            if ( myid < nthreads-1 ) myend = mybegin+stride_i;
            else myend = P->IA[P->row];
            for ( i = mybegin; i < myend; ++i ) {
                j = P->JA[i];
                P->JA[i] = cindex[j];
            }
        }
#endif
    }
    else {
        for ( i = 0; i < P->nnz; ++i ) {
            j = P->JA[i];
            P->JA[i] = cindex[j];
        }
    }
    
    // clean up
    fasp_mem_free(cindex);
    
    // Step 3. Truncate the prolongation operator to reduce cost
    fasp_amg_interp_trunc(P, param);
}

/**
 * \fn static void interp_STD (dCSRmat *A, ivector *vertices, dCSRmat *P, 
 *                             iCSRmat *S, AMG_param *param)
 *
 * \brief Standard interpolation
 *
 * \param A          Pointer to dCSRmat: the coefficient matrix (index starts from 0)
 * \param vertices   Indicator vector for the C/F splitting of the variables
 * \param P          Interpolation matrix (input: nnz pattern, output: prolongation)
 * \param S          Strong connection matrix
 * \param param      Pointer to AMG_param: AMG parameters
 *
 * \author Kai Yang, Xiaozhe Hu
 * \date   05/21/2012
 *
 * Modified by Chunsheng Feng, Zheng Li on 10/17/2012: add OMP support
 * Modified by Chensong Zhang on 05/15/2013: reconstruct the code
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 12/25/2013: cfsplitting of RS coarsening check C1 Criterion
 */
static void interp_STD (dCSRmat *A,
                        ivector *vertices,
                        dCSRmat *P,
                        iCSRmat *S,
                        AMG_param *param)
{
    const INT   row    = A->row;
    INT        *vec    = vertices->val;
    
    // local variables
    INT    i, j, k, l, m, index;
    REAL   alpha, factor, alN, alP;
    REAL   akk, akl, aik, aki;
    
    // indices for coarse neighbor node for every node
    INT  * cindex = (INT *)fasp_mem_calloc(row, sizeof(INT));
    
    // indices from column number to index in nonzeros in i-th row
    INT  * rindi  = (INT *)fasp_mem_calloc(2*row, sizeof(INT));
    
    // indices from column number to index in nonzeros in k-th row
    INT  * rindk  = (INT *)fasp_mem_calloc(2*row, sizeof(INT));
    
    // sums of strongly connected C neighbors
    REAL * csum   = (REAL *)fasp_mem_calloc(row, sizeof(REAL));

#if RS_C1
    // sums of all neighbors except ISPT
    REAL * psum   = (REAL *)fasp_mem_calloc(row, sizeof(REAL));
#endif    
    // sums of all neighbors
    REAL * nsum   = (REAL *)fasp_mem_calloc(row, sizeof(REAL));
    
    // diagonal entries
    REAL * diag   = (REAL *)fasp_mem_calloc(row, sizeof(REAL));
    
    // coefficents hat a_ij for relevant CGPT of the i-th node
    REAL * Ahat   = (REAL *)fasp_mem_calloc(row, sizeof(REAL));
    
    // Step 0. Prepare diagonal, Cs-sum, and N-sum
    fasp_iarray_set(row, cindex, -1);
    fasp_array_set(row, csum, 0.0);
    fasp_array_set(row, nsum, 0.0);
    
    for ( i = 0; i < row; i++ ) {
        
        // set flags for strong-connected C nodes
        for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
            k = S->JA[j];
            if ( vec[k] == CGPT ) cindex[k] = i;
        }
        
        for ( j = A->IA[i]; j < A->IA[i+1]; j++ ) {
            k = A->JA[j];
            
            if ( cindex[k] == i ) csum[i] += A->val[j]; // strong C-couplings
            
            if ( k == i ) diag[i]  = A->val[j];
#if RS_C1
            else {
                nsum[i] += A->val[j];
                if ( vec[k] != ISPT ) {
                    psum[i] += A->val[j];
                }
#else
            else          nsum[i] += A->val[j];
#endif
            }
        }
        
    }
    
    // Step 1. Fill in values for interpolation operator P
    for ( i = 0; i < row; i++ ) {
        
        if ( vec[i] == FGPT ) {
#if RS_C1            
            alN = psum[i];
#else
            alN = nsum[i];
#endif
            alP = csum[i];
            
            // form the reverse indices for i-th row
            for ( j = A->IA[i]; j < A->IA[i+1]; j++ ) rindi[A->JA[j]] = j;
            
            // clean up Ahat for relevent nodes only
            for ( j = P->IA[i]; j < P->IA[i+1]; j++ ) Ahat[P->JA[j]] = 0.0;
            
            // set values of Ahat
            Ahat[i] = diag[i];
            
            for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
                
                k = S->JA[j]; aik = A->val[rindi[k]];
                
                if ( vec[k] == CGPT ) Ahat[k] += aik;
                
                else if ( vec[k] == FGPT ) {
                    
                    akk = diag[k];
                    
                    // form the reverse indices for k-th row
                    for ( m = A->IA[k]; m < A->IA[k+1]; m++ ) rindk[A->JA[m]] = m;
                    
                    factor = aik / akk;
                    
                    // visit the strong-connected C neighbors of k, compute
                    // Ahat in the i-th row, set aki if found
                    aki = 0.0;
#if 0               // modified by Xiaoqiang Yue 12/25/2013
                    for ( m = S->IA[k]; m < S->IA[k+1]; m++ ) {
                        l   = S->JA[m];
                        akl = A->val[rindk[l]];
                        if ( vec[l] == CGPT ) Ahat[l] -= factor * akl;
                        else if ( l == i ) {
                            aki = akl; Ahat[l] -= factor * aki;
                        }
                    } // end for m
#else
                    for ( m = A->IA[k]; m < A->IA[k+1]; m++ ) {
                        if ( A->JA[m] == i ) {
                            aki = A->val[m];
                            Ahat[i] -= factor * aki;
                        }
                    } // end for m
#endif
                    for ( m = S->IA[k]; m < S->IA[k+1]; m++ ) {
                        l   = S->JA[m];
                        akl = A->val[rindk[l]];
                        if ( vec[l] == CGPT ) Ahat[l] -= factor * akl;
                    } // end for m
                    
                    // compute Cs-sum and N-sum for Ahat
                    alN -= factor * (nsum[k]-aki+akk);
                    alP -= factor *  csum[k];
                    
                } // end if vec[k]
                
            } // end for j
            
            // How about positive entries? --Chensong
            alpha = alN/alP;
            for ( j = P->IA[i]; j < P->IA[i+1]; j++ ) {
                k = P->JA[j];
                P->val[j] = -alpha*Ahat[k]/Ahat[i];
            }
            
        }
        
        else if ( vec[i] == CGPT ) {
            P->val[P->IA[i]] = 1.0;
        }
        
    } // end for i
    
    // Step 2. Generate coarse level indices and set values of P.JA
    for ( index = i = 0; i < row; ++i ) {
        if ( vec[i] == CGPT ) cindex[i] = index++;
    }
    P->col = index;
    
#ifdef _OPENMP
#pragma omp parallel for private(i,j) if(P->IA[P->row]>OPENMP_HOLDS)
#endif
    for ( i = 0; i < P->IA[P->row]; ++i ) {
        j = P->JA[i];
        P->JA[i] = cindex[j];
    }
    
    // clean up
    fasp_mem_free(cindex);
    fasp_mem_free(rindi);
    fasp_mem_free(rindk);
    fasp_mem_free(nsum);

#if RS_C1
    fasp_mem_free(psum);
#endif

    fasp_mem_free(csum);
    fasp_mem_free(diag);
    fasp_mem_free(Ahat);
    
    // Step 3. Truncate the prolongation operator to reduce cost
    fasp_amg_interp_trunc(P, param);
}

/**
 * \fn static void get_nwidth (dCSRmat *A, INT *nbl_ptr, INT *nbr_ptr)
 *
 * \brief Get the bandwidth of the matrix
 *
 * \param A          Pointer to dCSRmat: the coefficient matrix (index starts from 0)
 * \param nbl_ptr    Left bandwidth
 * \param nbr_ptr    Right bandwidth
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   03/01/2011
 */
static void get_nwidth (dCSRmat *A,
                        INT *nbl_ptr,
                        INT *nbr_ptr )
{
#ifdef _OPENMP
    INT *IA = A->IA;
    INT *JA = A->JA;
    INT myid, mybegin, myend, nthreads;
    INT max_l, max_r;
    INT i, end_row_A, j;
    
    nthreads = FASP_GET_NUM_THREADS();
    
    INT *max_left_right = (INT *)fasp_mem_calloc(2*nthreads, sizeof(INT));
    
#pragma omp parallel for private(myid,mybegin,myend,max_l,max_r,i,end_row_A,j)
    for (myid = 0; myid < nthreads; myid ++)
    {
        FASP_GET_START_END(myid, nthreads, A->row, &mybegin, &myend);
        max_l = 0;
        max_r = 0;
        for (i = mybegin; i < myend; i ++) {
            end_row_A = IA[i+1];
            for (j = IA[i]; j < end_row_A; j ++) {
                max_l = MAX(i-JA[j], max_l);
                max_r = MAX(JA[j]-i, max_r);
            }
        }
        max_left_right[myid*2] = max_l;
        max_left_right[myid*2+1] = max_r;
    }
    max_l = max_left_right[0];
    max_r = max_left_right[1];
    for (i = 1; i < nthreads; i ++) {
        max_l = MAX(max_l, max_left_right[i*2]);
        max_r = MAX(max_r, max_left_right[i*2+1]);
    }
    fasp_mem_free(max_left_right);
    *nbl_ptr = max_l;
    *nbr_ptr = max_r;
#endif
}

/**
 * \fn static void mod_cindex (INT nrows, INT *cindex)
 *
 * \brief Modify coarse indices
 *
 * \param nrows      Length of coarse indices
 * \param cindex     Indices of coarse nodes
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   03/01/2011
 */
static void mod_cindex (INT nrows,
                        INT *cindex)
{
#ifdef _OPENMP
    INT myid, mybegin, myend, i, nthreads;
    
    nthreads = FASP_GET_NUM_THREADS();
    
#pragma omp parallel for private(myid,mybegin,myend,i)
    for (myid = 0; myid < nthreads; myid ++)
    {
        FASP_GET_START_END(myid, nthreads, nrows, &mybegin, &myend);
        if (myid == 0)
        {
            mybegin ++;
        }
        for (i = mybegin; i < myend; i ++)
        {
            if (cindex[i] < cindex[i-1])
            {
                cindex[i] = cindex[i-1];
            }
        }
    }
#endif
}

/**
 * \fn static void get_cindex (INT nrows, INT ncols, INT *cindex,
 *                             INT nbl_ysk, INT nbr_ysk, INT *CF_marker,
 *                             INT *icor_ysk)
 *
 * \brief Get indices of coarse nodes in fine grid.
 *
 * \param nrows      Length of cindex
 * \param ncols      Length of cindex
 * \param cindex     Indices of nodes in coarse grid
 * \param nbl_ysk    Left  bandwith
 * \param nbr_ysk    Right bandwith
 * \param CF_marker  C/F marker of the nodes, 1: coarse nodes, 0: fine nodes.
 * \param icor_ysk   Indices of coarse nodes in fine grid
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   03/01/2011
 */
static void get_cindex (INT nrows,
                        INT ncols,
                        INT *cindex,
                        INT nbl_ysk,
                        INT nbr_ysk,
                        INT *CF_marker,
                        INT *icor_ysk)
{
#ifdef _OPENMP
    INT myid, FiveMyid, mybegin, myend, min_A, max_A, i, first_f_node, min_P, max_P, myend_minus_one;
    INT lengthAA = 0, lengthPP = 0;
    INT nthreads;
    
    nthreads = FASP_GET_NUM_THREADS();
    
#pragma omp parallel for private(myid,FiveMyid,mybegin,myend,min_A,max_A,i,first_f_node,min_P,max_P,myend_minus_one) reduction(+: lengthAA,lengthPP)
    for (myid = 0; myid < nthreads; myid ++)
    {
        FiveMyid = myid * 5;
        FASP_GET_START_END(myid, nthreads, ncols, &mybegin, &myend);
        icor_ysk[FiveMyid] = mybegin;
        if (mybegin == myend) {
            lengthAA = 0;
            lengthPP = 0;
            icor_ysk[FiveMyid+1] = 0;
            icor_ysk[FiveMyid+3] = 0;
        }
        else {
            first_f_node = fasp_BinarySearch(cindex, mybegin, nrows);
            for (i = first_f_node; i > -1; i --) {
                if (cindex[i] != mybegin) {
                    break;
                }
            }
            min_A = i + 1;
            min_A = MAX(0, min_A-2*nbl_ysk);
            myend_minus_one = myend - 1;
            first_f_node = fasp_BinarySearch(cindex, myend_minus_one, nrows);
            for (i = first_f_node; i > -1; i --) {
                if (cindex[i] != myend_minus_one) {
                    max_A = i;
                    break;
                }
            }
            max_A = MIN(nrows, max_A+2*nbr_ysk+1);
            lengthAA = max_A - min_A + 2;
            icor_ysk[FiveMyid+1] = lengthAA;
            icor_ysk[FiveMyid+2] = min_A;
            for (i = min_A; i >= 0; i --) {
                if (CF_marker[i] == 0) {
                    first_f_node = i;
                    break;
                }
            }
            if (i == -1) {
                min_P = 0;
            }
            else {
                first_f_node -= nbl_ysk;
                if (first_f_node <= 0) {
                    min_P = 0;
                }
                else {
                    for (i = first_f_node; i >= 0; i --) {
                        if (CF_marker[i] == 1) {
                            min_P = cindex[i];
                            break;
                        }
                    }
                    if (i == -1) {
                        min_P = 0;
                    }
                }
            }
            for (i = max_A-1; i < nrows; i ++) {
                if (CF_marker[i] == 0) {
                    first_f_node = i;
                    break;
                }
            }
            if (i == nrows) {
                max_P = ncols;
            }
            else {
                first_f_node += nbr_ysk;
                if (first_f_node >= nrows) {
                    max_P = ncols;
                }
                else {
                    for (i = first_f_node; i < nrows; i ++) {
                        if (CF_marker[i] == 1) {
                            max_P = cindex[i] + 1;
                            break;
                        }
                    }
                    if (i == nrows) {
                        max_P = ncols;
                    }
                }
            }
            lengthPP = max_P - min_P + 2;
            icor_ysk[FiveMyid+3] = lengthPP;
            icor_ysk[FiveMyid+4] = min_P;
        }
    }
    icor_ysk[5*nthreads] = lengthAA;
    icor_ysk[5*nthreads+1] = lengthPP;
#endif
}

/**
 * \fn static void interp_DIR1 (dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param)
 *
 * \brief Direct interpolation
 *
 * \param A          Pointer to dCSRmat: the coefficient matrix (index starts from 0)
 * \param vertices   Indicator vector for the C/F splitting of the variables
 * \param Ptr        Interpolation matrix (input: nnz pattern, output: prolongation)
 * \param param      Pointer to AMG_param: AMG parameters
 *
 * \author Chunsheng Feng, Xiaoqiang Yue
 * \date   03/01/2011
 *
 * \note This function is redesigned for OpenMP efficiency. 
 * 
 * TODO: Need to be cleaned up! --Chensong
 */
static void interp_DIR1 (dCSRmat *A,
                         ivector *vertices,
                         dCSRmat *Ptr,
                         AMG_param *param,
                         INT *icor_ysk)
{
    REAL eps_tr = param->truncation_threshold;
    REAL amN, amP, apN, apP;
    REAL alpha, beta, aii=0;
    INT *vec = vertices->val;
    INT num_pcouple, idiag;
    
    INT i,j,k,l,index=0;
    INT begin_row, end_row;
    INT myid;
    INT mybegin;
    INT myend;
    INT *indexs = NULL;
    INT shift;
    INT nbl_ysk, nbr_ysk;
    
    /* Generate interpolation P */
    dCSRmat P=fasp_dcsr_create(Ptr->row,Ptr->col,Ptr->nnz);
    
    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP
    INT row = MIN(P.IA[P.row], A->row);
    if ( row > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    /* step 1: Find first the structure IA of P */
    fasp_iarray_cp(P.row+1, Ptr->IA, P.IA);
    
    /* step 2: Find the structure JA of P */
    fasp_iarray_cp(P.nnz, Ptr->JA, P.JA);
    
    /* step 3: Fill the data of P */
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid,mybegin,myend,i,begin_row,end_row,idiag,aii,amN,amP,apN,apP,num_pcouple,j,k,alpha,beta,l)
#endif
        for (myid = 0; myid < nthreads; myid++ )
        {
            FASP_GET_START_END(myid, nthreads, A->row, &mybegin, &myend);
            for (i=mybegin; i<myend; ++i)
            {
                begin_row=A->IA[i]; end_row=A->IA[i+1]-1;
                
                for(idiag=begin_row;idiag<=end_row;idiag++) {
                    if (A->JA[idiag]==i) {
                        aii=A->val[idiag];
                        break;
                    }
                }
                
                if(vec[i]==0)  // if node i is on fine grid
                {
                    amN=0, amP=0, apN=0, apP=0,  num_pcouple=0;
                    
                    for(j=begin_row;j<=end_row;++j)
                    {
                        if(j==idiag) continue;
                        
                        for(k=Ptr->IA[i];k<Ptr->IA[i+1];++k) {
                            if(Ptr->JA[k]==A->JA[j]) break;
                        }
                        
                        if(A->val[j]>0) {
                            apN+=A->val[j];
                            if(k<Ptr->IA[i+1]) {
                                apP+=A->val[j];
                                num_pcouple++;
                            }
                        }
                        else
                        {
                            amN+=A->val[j];
                            if(k<Ptr->IA[i+1]) {
                                amP+=A->val[j];
                            }
                        }
                    } // j
                    
                    alpha=amN/amP;
                    if(num_pcouple>0) {
                        beta=apN/apP;
                    }
                    else {
                        beta=0;
                        aii+=apN;
                    }
                    
                    for(j=P.IA[i];j<P.IA[i+1];++j)
                    {
                        k=P.JA[j];
                        for(l=A->IA[i];l<A->IA[i+1];l++)
                        {
                            if(A->JA[l]==k) break;
                        }
                        
                        if(A->val[l]>0)
                        {
                            P.val[j]=-beta*A->val[l]/aii;
                        }
                        else
                        {
                            P.val[j]=-alpha*A->val[l]/aii;
                        }
                    }
                }
                else if(vec[i]==2)
                {
                    // if node i is a special fine node
                }
                else
                {
                    // if node i is on coarse grid
                    P.val[P.IA[i]]=1;
                }
            }
        }
    }
    else {
        for(i=0;i<A->row;++i) {
            begin_row=A->IA[i]; end_row=A->IA[i+1]-1;
            
            for(idiag=begin_row;idiag<=end_row;idiag++) {
                if (A->JA[idiag]==i) {
                    aii=A->val[idiag];
                    break;
                }
            }
            
            if(vec[i]==0)  // if node i is on fine grid
            {
                amN=0, amP=0, apN=0, apP=0,  num_pcouple=0;
                
                for(j=begin_row;j<=end_row;++j)
                {
                    if(j==idiag) continue;
                    
                    for(k=Ptr->IA[i];k<Ptr->IA[i+1];++k) {
                        if(Ptr->JA[k]==A->JA[j]) break;
                    }
                    
                    if(A->val[j]>0) {
                        apN+=A->val[j];
                        if(k<Ptr->IA[i+1]) {
                            apP+=A->val[j];
                            num_pcouple++;
                        }
                    }
                    else
                    {
                        amN+=A->val[j];
                        if(k<Ptr->IA[i+1]) {
                            amP+=A->val[j];
                        }
                    }
                } // j
                
                alpha=amN/amP;
                if(num_pcouple>0) {
                    beta=apN/apP;
                }
                else {
                    beta=0;
                    aii+=apN;
                }
                
                for(j=P.IA[i];j<P.IA[i+1];++j)
                {
                    k=P.JA[j];
                    for(l=A->IA[i];l<A->IA[i+1];l++)
                    {
                        if(A->JA[l]==k) break;
                    }
                    
                    if(A->val[l]>0)
                    {
                        P.val[j]=-beta*A->val[l]/aii;
                    }
                    else
                    {
                        P.val[j]=-alpha*A->val[l]/aii;
                    }
                }
            }
            else if(vec[i]==2)  // if node i is a special fine node
            {
            }
            else // if node i is on coarse grid
            {
                P.val[P.IA[i]]=1;
            }
        }
    }
    fasp_mem_free(Ptr->IA);
    fasp_mem_free(Ptr->JA);
    fasp_mem_free(Ptr->val);
    
    INT *cindex;
    
    cindex=(INT*)fasp_mem_calloc(A->row, sizeof(INT));
    
#if CHMEM_MODE
    total_alloc_mem += (A->row)*sizeof(INT);
#endif
    
    // The following is one of OPTIMAL parts ...0802...
    // Generate cindex in parallel
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp master
        {
            indexs = (INT *)fasp_mem_calloc(nthreads, sizeof(INT));
        }
#pragma omp parallel for private(myid, mybegin, myend, index, i)
#endif
        for (myid = 0; myid < nthreads; myid ++)
        {
            FASP_GET_START_END(myid, nthreads, A->row, &mybegin, &myend);
            index = 0;
            for (i=mybegin;i<myend;++i) {
                if(vec[i]==1)
                {
                    cindex[i]=index;
                    index++;
                }
            }
            indexs[myid] = index;
        }
        for (i = 1; i < nthreads; i ++) {
            indexs[i] += indexs[i-1];
        }
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin, myend, shift, i)
#endif
        for (myid = 0; myid < nthreads; myid ++)
        {
            FASP_GET_START_END(myid, nthreads, A->row, &mybegin, &myend);
            shift = 0;
            if (myid > 0) {
                shift = indexs[myid-1];
            }
            for (i=mybegin;i<myend;++i) {
                if(vec[i]==1)
                {
                    cindex[i] += shift;
                }
            }
        }
        fasp_mem_free(indexs);
    }
    else {
        index=0;
        for(i=0;i<A->row;++i) {
            if(vec[i]==1)
            {
                cindex[i]=index;
                index++;
            }
        }
    }
    
    if (use_openmp) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin,myend,i,j)
#endif
        for (myid = 0; myid < nthreads; myid++ )
        {
            FASP_GET_START_END(myid, nthreads, P.IA[P.row], &mybegin, &myend);
            for (i=mybegin; i<myend; ++i)
            {
                j=P.JA[i];
                P.JA[i]=cindex[j];
            }
        }
    }
    else {
        for(i=0;i<P.IA[P.row];++i)
        {
            j=P.JA[i];
            P.JA[i]=cindex[j];
        }
    }
    //    if (Ptr->col > OPENMP_HOLDS) {
    // The following is another OPTIMAL part ...0802...
    ////////////////////////////////////////////////////////////////////////////////
    get_nwidth(A, &nbl_ysk, &nbr_ysk);
    mod_cindex(A->row, cindex);
    get_cindex(A->row, Ptr->col, cindex, nbl_ysk, nbr_ysk, vec, icor_ysk);
    ////////////////////////////////////////////////////////////////////////////////
    //    }
    fasp_mem_free(cindex);
    
    /* Truncation of interpolation */
    REAL Min_neg, Max_pos;
    REAL Sum_neg, Sum_pos;
    REAL TSum_neg, TSum_pos;
    INT num_neg, num_pos;
    INT num_nonzero=0;
    
    
    Ptr->val=(REAL*)fasp_mem_calloc(P.IA[Ptr->row],sizeof(REAL));
    
#if CHMEM_MODE
    total_alloc_mem += (P.IA[Ptr->row])*sizeof(REAL);
#endif
    Ptr->JA=(INT*)fasp_mem_calloc(P.IA[Ptr->row],sizeof(INT));
#if CHMEM_MODE
    total_alloc_mem += (P.IA[Ptr->row])*sizeof(INT);
#endif
    Ptr->IA=(INT*)fasp_mem_calloc(Ptr->row+1, sizeof(INT));
#if CHMEM_MODE
    total_alloc_mem += (Ptr->row+1)*sizeof(INT);
#endif
    
    INT index1=0, index2=0;
    for(i=0;i<P.row;++i)
    {
        Min_neg=0;
        Max_pos=0;
        Sum_neg=0;
        Sum_pos=0;
        TSum_neg=0;
        TSum_pos=0;
        num_neg=0;
        num_pos=0;
        
        Ptr->IA[i]-=num_nonzero;
        
        for(j=P.IA[i];j<P.IA[i+1];++j)
        {
            if(P.val[j]<0)
            {
                Sum_neg+=P.val[j];
                if(P.val[j]<Min_neg)
                {
                    Min_neg=P.val[j];
                }
            }
            
            if(P.val[j]>0)
            {
                Sum_pos+=P.val[j];
                if(P.val[j]>Max_pos)
                {
                    Max_pos=P.val[j];
                }
            }
        }
        
        for(j=P.IA[i];j<P.IA[i+1];++j)
        {
            if(P.val[j]<0)
            {
                if(P.val[j]>Min_neg*eps_tr)
                {
                    num_neg++;
                }
                else
                {
                    num_nonzero--;
                }
            }
            
            if(P.val[j]>0)
            {
                if(P.val[j]<Max_pos*eps_tr)
                {
                    num_pos++;
                }
                else
                {
                    num_nonzero--;
                }
            }
        }
        
        // step 2: Find the structure JA and fill the data A of Ptr
        for(j=P.IA[i];j<P.IA[i+1];++j)
        {
            if(P.val[j]<0)
            {
                if(!(P.val[j]>Min_neg*eps_tr))
                {
                    Ptr->JA[index1]=P.JA[j];
                    TSum_neg+=P.val[j];
                    index1++;
                }
            }
            
            if(P.val[j]>0)
            {
                if(!(P.val[j]<Max_pos*eps_tr))
                {
                    Ptr->JA[index1]=P.JA[j];
                    TSum_pos+=P.val[j];
                    index1++;
                }
            }
        }
        
        // step 3: Fill the data A of Ptr
        for(j=P.IA[i];j<P.IA[i+1];++j)
        {
            if(P.val[j]<0)
            {
                if(!(P.val[j]>Min_neg*eps_tr))
                {
                    Ptr->val[index2]=P.val[j]/TSum_neg*Sum_neg;
                    index2++;
                }
            }
            
            if(P.val[j]>0)
            {
                if(!(P.val[j]<Max_pos*eps_tr))
                {
                    Ptr->val[index2]=P.val[j]/TSum_pos*Sum_pos;
                    index2++;
                }
            }
        }
    }
    Ptr->IA[P.row]-=num_nonzero;
    Ptr->nnz=Ptr->IA[Ptr->row];
    
    Ptr->JA=(INT*)fasp_mem_realloc(Ptr->JA, Ptr->IA[Ptr->row]*sizeof(INT));
    Ptr->val=(REAL*)fasp_mem_realloc(Ptr->val, Ptr->IA[Ptr->row]*sizeof(REAL));
    
    fasp_mem_free(P.IA);
    fasp_mem_free(P.JA);
    fasp_mem_free(P.val);    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
