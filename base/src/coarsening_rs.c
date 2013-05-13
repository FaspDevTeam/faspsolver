/*! \file coarsening_rs.c
 *  \brief Coarsening with a modified Ruge-Stuben strategy.
 *
 *  \note Ref Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller
 *            Appendix P475 A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *            Academic Press Inc., San Diego, CA, 2001.
 */

#include <omp.h>

#include "fasp.h"
#include "fasp_functs.h"

// Private routines for RS coarsening
static void find_strong_couple    (dCSRmat *, iCSRmat *, AMG_param *);
static INT  form_coarse_level_std (dCSRmat *, iCSRmat *, ivector *, INT);
static INT  remove_ff_connections (iCSRmat *, ivector *, INT, INT);
static INT  form_coarse_level_agg (dCSRmat *, iCSRmat *, ivector *, INT, INT);
static void form_P_pattern_dir    (dCSRmat *, iCSRmat *, ivector *, INT, INT);
static void form_P_pattern_std    (dCSRmat *, iCSRmat *, ivector *, INT, INT);

#include "linklist.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_coarsening_rs (dCSRmat *A, ivector *vertices, dCSRmat *P,
 *                                 iCSRmat *S, AMG_param *param)
 *
 * \brief Standard and aggressive coarsening schemes
 *
 * \param A          Coefficient matrix, the index starts from zero
 * \param vertices   Indicator vector for the C/F splitting of the variables
 * \param P          Interpolation matrix (nonzero pattern only)
 * \param S          Strong connection matrix
 * \param param      AMG parameters
 *
 * \return           SUCCESS or error message
 *
 * \author Xuehai Huang, Chensong Zhang, Xiaozhe Hu, Ludmil Zikatanov
 * \date   09/06/2010
 *
 * \note vertices = 0: fine; 1: coarse; 2: isolated or special
 *
 * Modified by Xiaozhe Hu on 05/23/2011: add strength matrix as an argument
 * Modified by Chensong Zhang on 04/21/2013
 * Modified by Xiaozhe Hu on 04/24/2013: modfiy aggressive coarsening
 * Mofified by Chensong Zhang on 04/28/2013: remove linked list
 * Mofified by Chensong Zhang on 05/11/2013: restructure the code
 */
INT fasp_amg_coarsening_rs (dCSRmat *A,
                            ivector *vertices,
                            dCSRmat *P,
                            iCSRmat *S,
                            AMG_param *param)
{
    const SHORT coarse_type = param->coarsening_type;
    const INT   agg_path    = param->aggressive_path;
    const INT   row         = A->row;
    
    // local variables
    SHORT       interp_type = param->interpolation_type;
    INT         col         = 0;
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_coarsening_rs ...... [Start]\n");
#endif
    
#if DEBUG_MODE
    printf("### DEBUG: Step 1. Find strong connections ......\n");
#endif
    
    // make sure standard interp is used for aggressive coarsening
    if ( coarse_type == COARSE_AC ) {
        param->interpolation_type = interp_type = INTERP_STD;
    }
    
    // find stong couplings and return them in S
    find_strong_couple(A, S, param);
    if ( S->nnz == 0 ) {
        printf("### ERROR: Fail to find strong connections!\n");        
        return RUN_FAIL;
    }
    
#if DEBUG_MODE
    printf("### DEBUG: Step 2. C/F splitting ......\n");
#endif
    
    switch ( coarse_type ) {
            
        case COARSE_RS: // modified Ruge-Stuben
            col = form_coarse_level_std(A, S, vertices, row);
            break;
            
        case COARSE_AC: // aggressive coarsening
            col = form_coarse_level_agg(A, S, vertices, row, agg_path);
            break;
            
        case COARSE_CR: // compatible relaxation (Need to be modified --Chensong)
            col = fasp_amg_coarsening_cr(0, A->row-1, A, vertices, param);
            break;
            
        default:
            printf("### ERROR: Coarsening type %d is not recognized!\n", coarse_type);
            return ERROR_AMG_COARSE_TYPE;
            
    }
    
    if ( col <= 0 ) {
        printf("### ERROR: Fail to find C variables!\n");
        return RUN_FAIL;
    }

#if DEBUG_MODE
    printf("### DEBUG: Step 3. Find support of C points ......\n");
#endif
    
    switch ( interp_type ) {
            
        case INTERP_DIR: // direct interpolation
        case INTERP_ENG: // energy-min interpolation
            col = remove_ff_connections(S, vertices, row, col);
            form_P_pattern_dir(P, S, vertices, row, col);
            break;
            
        case INTERP_STD: // standard interpolaiton
            form_P_pattern_std(P, S, vertices, row, col);
            break;
            
        default:
            printf("### ERROR: Interpoltion type %d is not recognized!\n", interp_type);
            return ERROR_AMG_INTERP_TYPE;
            
    }
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_coarsening_rs ...... [Finish]\n");
#endif
    
    return SUCCESS;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void find_strong_couple (dCSRmat *A, iCSRmat *S, AMG_param *param)
 *
 * \brief Generate the set of all strong negative couplings
 *
 * \param A          Coefficient matrix, the index starts from zero
 * \param S          Strong connection matrix
 * \param param      AMG parameters
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/25/2012: add OMP support
 * Mofified by Chensong Zhang on 05/11/2013: restructure the code
 */
static void find_strong_couple (dCSRmat *A,
                                iCSRmat *S,
                                AMG_param *param )
{
    const REAL  max_row_sum = param->max_row_sum;
    const REAL  epsilon_str = param->strong_threshold;
    const INT   row = A->row, col = A->col, row1 = row+1;
    const INT   nnz = A->nnz;
    
    INT  *ia = A->IA, *ja = A->JA;
    REAL *aj = A->val;

    // local variables
    INT   index, i, j, begin_row, end_row;
    REAL  row_min, row_sum;

    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP
    if ( row > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    // get the diagnal entry of A
    dvector diag; fasp_dcsr_getdiag(0, A, &diag);
        
    // copy the structure of A to S
    S->row = row; S->col = col; S->nnz = nnz; S->val = NULL;
    S->IA = (INT *)fasp_mem_calloc(row1, sizeof(INT));
    S->JA = (INT *)fasp_mem_calloc(nnz,  sizeof(INT));
    fasp_iarray_cp(row1, ia, S->IA); // init, no necessary --Chensong
    fasp_iarray_cp(nnz,  ja, S->JA);
    
    if ( use_openmp ) {
        
        INT mybegin, myend, myid;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin,myend,i,row_min,row_sum,begin_row,end_row,j)
#endif
        for ( myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
            for ( i = mybegin; i < myend; i++) {
                
                // Compute most negative entry in each row and row sum
                row_min = row_sum = 0.0;
                begin_row = ia[i]; end_row = ia[i+1];
                for ( j = begin_row; j < end_row; j++ ) {
                    row_min  = MIN(row_min, aj[j]);
                    row_sum += aj[j];
                }
                
                // Find diagonal entries of S and remove them later
                for ( j = begin_row; j < end_row; j++ ) {
                    if ( ja[j] == i ) { S->JA[j] = -1; break; }
                }
                
                // If sum_{j=1}^n a_{ij} > max_row_sum * |a_{ii}, mark entile row as
                // weak dependencies
                if ( max_row_sum < 1.0 && ABS(row_sum) > max_row_sum*ABS(diag.val[i]) ) {
                    //(row_sum > max_row_sum) ) {
                    for ( j = begin_row; j < end_row; j++ ) S->JA[j] = -1;
                }
                else {
                    for ( j = begin_row; j < end_row; j++) {
                        // If a_{ij} >= \epsilon_{str} * \min a_{ij}, the connection
                        // j->i is set to be weak; positive entries result in weak
                        // connections
                        if ( A->val[j] >= epsilon_str*row_min ) S->JA[j] = -1;
                    }
                }
                
            } // end for i
        } // end for myid
        
    }
    
    else {
    
        for ( i = 0; i < row; ++i ) {
            
            // Compute most negative entry in each row and row sum
            row_min = row_sum = 0.0;            
            begin_row = ia[i]; end_row = ia[i+1];
            for ( j = begin_row; j < end_row; j++ ) {
                row_min  = MIN(row_min, aj[j]);
                row_sum += aj[j];
            }
            
            // Find diagonal entries of S and remove them later
            for ( j = begin_row; j < end_row; j++ ) {
                if ( ja[j] == i ) { S->JA[j] = -1; break; }
            }
            
            // If sum_{j=1}^n a_{ij} > max_row_sum * |a_{ii}, mark entile row as
            // weak dependencies
            if ( max_row_sum < 1.0 && ABS(row_sum) > max_row_sum*ABS(diag.val[i]) ) {
                for ( j = begin_row; j < end_row; j++ ) S->JA[j] = -1;
            }
            else {
                for ( j = begin_row; j < end_row; j++ ) {
                    // If a_{ij} >= \epsilon_{str} * \min a_{ij}, the connection
                    // j->i is set to be weak; positive entries result in weak
                    // connections
                    if ( A->val[j] >= epsilon_str*row_min ) S->JA[j] = -1;
                }
            }
        } // end for i
    }
    
    // compress S: remove weak connections and form strong coupling matrix
    index = 0;
    for ( i = 0; i < row; ++i ) {
        S->IA[i] = index;
        begin_row = ia[i]; end_row = ia[i+1]-1;
        for ( j = begin_row; j<= end_row; j++ ) {
            if ( S->JA[j] > -1 ) S->JA[index++] = S->JA[j];
        }
    }
    
    if ( index > 0 ) {
        S->nnz = S->IA[row] = index;
    }
    else { // No strong coupling!!!
        S->nnz = 0; S->JA = NULL;
    }
    
    fasp_dvec_free(&diag);
}

/**
 * \fn static void find_strong_couple_agg1 (dCSRmat *A, iCSRmat *S, iCSRmat *Sh, 
 *                                          ivector *vertices, ivector *CGPT_index,
 *                                          ivector *CGPT_rindex)
 *
 * \brief Generate the set of all strong negative or absolute couplings using 
 *        aggressive coarsening A1 or A2
 *
 * \param A            Coefficient matrix, the index starts from zero
 * \param S            Strong connection matrix
 * \param Sh           Strong couplings matrix between coarse grid points
 * \param vertices     Type of variables--C/F splitting
 * \param CGPT_index   Index of CGPT from CGPT to all points
 * \param CGPT_rindex  Index of CGPT from all points to CGPT
 *
 * \author Kai Yang, Xiaozhe Hu
 * \date   09/06/2010
 *
 * \note The difference between find_strong_couple_agg1 and find_strong_couple_agg2
 *       is that find_strong_couple_agg1 uses path 1 to detetermine strongly coupled
 *       C points while find_strong_couple_agg2 uses path 2 to determinestrongly
 *       coupled C points. Usually find_strong_couple_agg1 gives more aggresive
 *       coarsening!
 *
 * Mofified by Chensong Zhang on 05/11/2013: restructure the code
 */
static void find_strong_couple_agg1 (dCSRmat *A,
                                     iCSRmat *S,
                                     iCSRmat *Sh,
                                     ivector *vertices,
                                     ivector *CGPT_index,
                                     ivector *CGPT_rindex)
{
    
    INT   i,j,k;
    INT   num_c,count,ci,cj,ck,fj,cck;
    INT  *cp_index, *cp_rindex, *times_visited, *vec=vertices->val;
    
    CGPT_rindex->row=A->row;
    CGPT_rindex->val=(INT*)fasp_mem_calloc(vertices->row,sizeof(INT));// for the reverse indexing of coarse grid points
    cp_rindex=CGPT_rindex->val;
    
    //count the number of coarse grid points
    num_c=0;
    for(i=0;i<vertices->row;i++)
    {
        if(vec[i]==CGPT) num_c++;
    }
    
    CGPT_index->row=num_c;
    
    //generate coarse grid point index
    CGPT_index->val=(INT *)fasp_mem_calloc(num_c,sizeof(INT));
    cp_index=CGPT_index->val;
    j=0;
    for (i=0;i<vertices->row;i++) {
        if(vec[i]==CGPT) {
            cp_index[j]=i;
            cp_rindex[i]=j;
            j++;
        }
    }
    
    Sh->row=num_c;
    Sh->col=num_c;
    Sh->val=NULL;
    Sh->JA=NULL;
    Sh->IA=(INT*)fasp_mem_calloc(Sh->row+1,sizeof(INT));
    
    times_visited=(INT*)fasp_mem_calloc(num_c,sizeof(INT)); // record the number of times some coarse point is visited
    
    //for (i=0; i<num_c; i++) times_visited[i]=0;
    memset(times_visited, 0, sizeof(INT)*num_c);
    
    // step 1: Find first the structure IA of Sh
    Sh->IA[0]=0;
    
    for(ci=0;ci<Sh->row;ci++)
    {
        count=0; // count the the number of coarse point that i is strongly connected to w.r.t. (p,2)
        i=cp_index[ci];//find the index of the ci-th coarse grid point
        
        //visit all the fine neighbors that ci is strongly connected to
        for(j=S->IA[i];j<S->IA[i+1];j++)
        {
            
            fj=S->JA[j];
            
            if(vec[fj]==CGPT&&fj!=i)
            {
                cj=cp_rindex[fj];
                
                if(times_visited[cj]!=ci+1)
                {
                    //newly visited
                    times_visited[cj]=ci+1;//marked as strongly connected from ci
                    count++;
                }
                
            }
            else if(vec[fj]==FGPT) // it is a fine grid point,
            {
                
                //find all the coarse neighbors that fj is strongly connected to
                for(k=S->IA[fj];k<S->IA[fj+1];k++)
                {
                    ck=S->JA[k];
                    
                    if(vec[ck]==CGPT&&ck!=i)// it is a coarse grid point
                    {
                        if(cp_rindex[ck]>=num_c) printf("find_strong_couple_agg1: index exceed bound!\n");
                        
                        cck=cp_rindex[ck];
                        
                        if (times_visited[cck]!=ci+1) {
                            //newly visited
                            times_visited[cck]=ci+1;//marked as strongly connected from ci
                            count++;
                        }
                        
                        /*
                         if (times_visited[cck] == ci+1){
                         
                         }
                         else if (times_visited[cck] == -ci-1){
                         times_visited[cck]=ci+1;//marked as strongly connected from ci
                         count++;
                         }
                         else{
                         times_visited[cck]=-ci-1;//marked as visited
                         }
                         */
                        
                    }//end if
                }//end for k
                
            }//end if
        }//end for j
        
        Sh->IA[ci+1]=Sh->IA[ci]+count;
        
    }//end for i
    
    
    // step 2: Find JA of Sh
    
    for (i=0; i<num_c; i++) times_visited[i]=0; // clean up times_visited
    Sh->nnz=Sh->IA[Sh->row];
    Sh->JA=(INT*)fasp_mem_calloc(Sh->nnz,sizeof(INT));
    
    for(ci=0;ci<Sh->row;ci++)
    {
        i=cp_index[ci]; //find the index of the i-th coarse grid point
        count=Sh->IA[ci]; //count for coarse points
        
        //visit all the fine neighbors that ci is strongly connected to
        for(j=S->IA[i];j<S->IA[i+1];j++)
        {
            fj=S->JA[j];
            if(vec[fj]==CGPT&&fj!=i)
            {
                cj=cp_rindex[fj];
                if(times_visited[cj]!=ci+1)
                {
                    //newly visited
                    times_visited[cj]=ci+1;
                    Sh->JA[count]=cj;
                    count++;
                }
                
                
            }
            else if(vec[fj]==FGPT) // it is a fine grid point,
            {
                //find all the coarse neighbors that fj is strongly connected to
                for(k=S->IA[fj];k<S->IA[fj+1];k++)
                {
                    ck=S->JA[k];
                    if(vec[ck]==CGPT&&ck!=i)// it is a coarse grid point
                    {
                        cck=cp_rindex[ck];
                        
                        
                        if(times_visited[cck]!=ci+1)
                        {
                            //newly visited
                            times_visited[cck]=ci+1;
                            Sh->JA[count]=cck;
                            count++;
                        }
                        
                        /*
                         if (times_visited[cck] == ci+1){
                         
                         }
                         else if (times_visited[cck] == -ci-1){
                         times_visited[cck]=ci+1;
                         Sh->JA[count]=cck;
                         count++;
                         }
                         else {
                         times_visited[cck]=-ci-1;
                         }
                         */
                        
                    }//end if
                }//end for k
                
            }//end if
        }//end for j
        if(count!=Sh->IA[ci+1]) printf("find_strong_couple_agg1: inconsistency in number of nonzeros values\n ");
    }//end for ci
    fasp_mem_free(times_visited);
}

/**
 * \fn static void find_strong_couple_agg2 (dCSRmat *A, iCSRmat *S, iCSRmat *Sh,
 *                                          ivector *vertices, ivector *CGPT_index,
 *                                          ivector *CGPT_rindex)
 *
 * \brief Generate the set of all strong negative or absolute couplings using
 *        aggressive coarsening A1 or A2
 *
 * \param A            Coefficient matrix, the index starts from zero
 * \param S            Strong connection matrix
 * \param Sh           Strong couplings matrix between coarse grid points
 * \param vertices     Type of variables--C/F splitting
 * \param CGPT_index   Index of CGPT from CGPT to all points
 * \param CGPT_rindex  Index of CGPT from all points to CGPT
 *
 * \author Xiaozhe Hu
 * \date   04/24/2013
 *
 * \note The difference between find_strong_couple_agg1 and find_strong_couple_agg2 
 *       is that find_strong_couple_agg1 uses path 1 to detetermine strongly coupled 
 *       C points while find_strong_couple_agg2 uses path 2 to determinestrongly 
 *       coupled C points. Usually find_strong_couple_agg1 gives more aggresive 
 *       coarsening!
 *
 * Mofified by Chensong Zhang on 05/11/2013: restructure the code
 */
static void find_strong_couple_agg2 (dCSRmat *A,
                                     iCSRmat *S,
                                     iCSRmat *Sh,
                                     ivector *vertices,
                                     ivector *CGPT_index,
                                     ivector *CGPT_rindex)
{
    
    INT   i,j,k;
    INT   num_c,count,ci,cj,ck,fj,cck;
    INT  *cp_index, *cp_rindex, *times_visited, *vec=vertices->val;
    
    CGPT_rindex->row=A->row;
    CGPT_rindex->val=(INT*)fasp_mem_calloc(vertices->row,sizeof(INT));// for the reverse indexing of coarse grid points
    cp_rindex=CGPT_rindex->val;
    
    //count the number of coarse grid points
    num_c=0;
    for(i=0;i<vertices->row;i++)
    {
        if(vec[i]==CGPT) num_c++;
    }
    
    CGPT_index->row=num_c;
    
    //generate coarse grid point index
    CGPT_index->val=(INT *)fasp_mem_calloc(num_c,sizeof(INT));
    cp_index=CGPT_index->val;
    j=0;
    for (i=0;i<vertices->row;i++) {
        if(vec[i]==CGPT) {
            cp_index[j]=i;
            cp_rindex[i]=j;
            j++;
        }
    }
    
    Sh->row=num_c;
    Sh->col=num_c;
    Sh->val=NULL;
    Sh->JA=NULL;
    Sh->IA=(INT*)fasp_mem_calloc(Sh->row+1,sizeof(INT));
    
    times_visited=(INT*)fasp_mem_calloc(num_c,sizeof(INT)); // record the number of times some coarse point is visited
    
    //for (i=0; i<num_c; i++) times_visited[i]=0;
    memset(times_visited, 0, sizeof(INT)*num_c);
    
    // step 1: Find first the structure IA of Sh
    Sh->IA[0]=0;
    
    for(ci=0;ci<Sh->row;ci++)
    {
        count=0; // count the the number of coarse point that i is strongly connected to w.r.t. (p,2)
        i=cp_index[ci];//find the index of the ci-th coarse grid point
        
        //visit all the fine neighbors that ci is strongly connected to
        for(j=S->IA[i];j<S->IA[i+1];j++)
        {
            
            fj=S->JA[j];
            
            if(vec[fj]==CGPT&&fj!=i)
            {
                cj=cp_rindex[fj];
                
                if(times_visited[cj]!=ci+1)
                {
                    //newly visited
                    times_visited[cj]=ci+1;//marked as strongly connected from ci
                    count++;
                }
                
            }
            else if(vec[fj]==FGPT) // it is a fine grid point,
            {
                
                //find all the coarse neighbors that fj is strongly connected to
                for(k=S->IA[fj];k<S->IA[fj+1];k++)
                {
                    ck=S->JA[k];
                    
                    if(vec[ck]==CGPT&&ck!=i)// it is a coarse grid point
                    {
                        if(cp_rindex[ck]>=num_c) printf("find_strong_couple_agg2: index exceed bound!\n");
                        
                        cck=cp_rindex[ck];
                        
                        if (times_visited[cck] == ci+1){
                            
                        }
                        else if (times_visited[cck] == -ci-1){
                            times_visited[cck]=ci+1;//marked as strongly connected from ci
                            count++;
                        }
                        else{
                            times_visited[cck]=-ci-1;//marked as visited
                        }
                        
                    }//end if
                }//end for k
                
            }//end if
        }//end for j
        
        Sh->IA[ci+1]=Sh->IA[ci]+count;
        
    }//end for i
    
    
    // step 2: Find JA of Sh
    
    for (i=0; i<num_c; i++) times_visited[i]=0; // clean up times_visited
    Sh->nnz=Sh->IA[Sh->row];
    Sh->JA=(INT*)fasp_mem_calloc(Sh->nnz,sizeof(INT));
    
    for(ci=0;ci<Sh->row;ci++)
    {
        i=cp_index[ci]; //find the index of the i-th coarse grid point
        count=Sh->IA[ci]; //count for coarse points
        
        //visit all the fine neighbors that ci is strongly connected to
        for(j=S->IA[i];j<S->IA[i+1];j++)
        {
            fj=S->JA[j];
            if(vec[fj]==CGPT&&fj!=i)
            {
                cj=cp_rindex[fj];
                if(times_visited[cj]!=ci+1)
                {
                    //newly visited
                    times_visited[cj]=ci+1;
                    Sh->JA[count]=cj;
                    count++;
                }
                
                
            }
            else if(vec[fj]==FGPT) // it is a fine grid point,
            {
                //find all the coarse neighbors that fj is strongly connected to
                for(k=S->IA[fj];k<S->IA[fj+1];k++)
                {
                    ck=S->JA[k];
                    if(vec[ck]==CGPT&&ck!=i)// it is a coarse grid point
                    {
                        cck=cp_rindex[ck];
                        
                        if (times_visited[cck] == ci+1){
                            
                        }
                        else if (times_visited[cck] == -ci-1){
                            times_visited[cck]=ci+1;
                            Sh->JA[count]=cck;
                            count++;
                        }
                        else {
                            times_visited[cck]=-ci-1;
                        }
                        
                    }//end if
                }//end for k
                
            }//end if
        }//end for j
        if(count!=Sh->IA[ci+1]) printf("find_strong_couple_agg2: inconsistency in number of nonzeros values\n ");
    }//end for ci
    fasp_mem_free(times_visited);
}

/**
 * \fn static INT form_coarse_level_std (dCSRmat *A, iCSRmat *S, ivector *vertices,
 *                                       INT row)
 *
 * \brief Find coarse level variables (C/F splitting): standard
 *
 * \param A            Coefficient matrix, the index starts from zero
 * \param S            Strong connection matrix
 * \param vertices     Indicator vector for the C/F splitting of the variables
 * \param row          Number of rows of P
 *
 * \return Number of cols of P
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 *
 * \note Coarsening Phase ONE: find coarse level points
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/24/2012: add OMP support
 * Modified by Chensong Zhang on 07/06/2012: fix a data type bug
 * Mofified by Chensong Zhang on 05/11/2013: restructure the code
 */
static INT form_coarse_level_std (dCSRmat *A,
                                  iCSRmat *S,
                                  ivector *vertices,
                                  INT row)
{
    // local variables
    INT col = 0;
    INT maxlambda, maxnode, num_left = 0;
    INT measure, newmeas;
    INT *ia = A->IA, *vec = vertices->val;
    INT i, j, k, l;
    INT myid, mybegin, myend;
    
    INT *work = (INT*)fasp_mem_calloc(3*row,sizeof(INT));
    INT *lists = work, *where = lists+row, *lambda = where+row;
    
    LinkList LoL_head = NULL, LoL_tail = NULL, list_ptr = NULL;
        
    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP
    if ( row > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    iCSRmat ST; fasp_icsr_trans(S, &ST);

    // 1. Initialize lambda
    if ( use_openmp ) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin,myend,i)
#endif
        for ( myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
            for ( i = mybegin; i < myend; i++ ) lambda[i] = ST.IA[i+1] - ST.IA[i];
        }
    }
    else {
        for ( i = 0; i < row; ++i ) lambda[i] = ST.IA[i+1] - ST.IA[i];
    }
    
    // 2. Before C/F splitting algorithm starts, filter out the variables which
    //    have no connections at all and mark them as special F-variables.
    if ( use_openmp ) {
        
#ifdef _OPENMP
#pragma omp parallel for reduction(+:num_left) private(myid, mybegin, myend, i)
#endif
        for ( myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
            for ( i = mybegin; i < myend; i++ ) {
                if ( (ia[i+1]-ia[i]) <= 1 ) {
                    vec[i] = ISPT; // set i as an ISOLATED fine node
                    lambda[i] = 0;
                }
                else {
                    vec[i] = UNPT; // set i as a undecided node
                    num_left++;
                }
            }
        } // end for myid
        
    }
    else {
        
        for ( i = 0; i < row; ++i ) {
            if ( (ia[i+1]-ia[i]) <= 1 ) {
                vec[i] = ISPT; // set i as an ISOLATED fine node
                lambda[i] = 0;
            }
            else {
                vec[i] = UNPT; // set i as a undecided node
                num_left++;
            }
        } // end for i
        
    }
    
    // 3. Form linked list for lambda (max to min)
    for ( i = 0; i < row; ++i ) {
        
        if ( vec[i] == ISPT ) continue; // skip isolated variables
            
        measure = lambda[i];
        
        if ( measure > 0 ) {
            enter_list(&LoL_head, &LoL_tail, lambda[i], i, lists, where);
        }
        else {
            
            if ( measure < 0 ) printf("### WARNING: Negative lambda[%d]!\n", i);
            
            // Set variables with nonpositvie measure as F-variables
            vec[i] = FGPT; // no strong connections, set i as fine node
            --num_left;

            // Update lambda and linked list after i->F
            for ( k = S->IA[i]; k < S->IA[i+1]; ++k ) {
                j = S->JA[k];
                if ( vec[j] == ISPT ) continue; // skip isolate variables
                if ( j < i ) { // why <? --Chensong
                    newmeas = lambda[j];
                    if ( newmeas > 0 ) {
                        remove_node(&LoL_head, &LoL_tail, newmeas, j, lists, where);
                    }
                    newmeas = ++(lambda[j]);
                    enter_list(&LoL_head, &LoL_tail, newmeas, j, lists, where);
                }
                else {
                    newmeas = ++(lambda[j]);
                }
            }
            
        } // end if measure
        
    } // end for i
    
    // 4. Main loop
    while ( num_left > 0 ) {
        
        // pick $i\in U$ with $\max\lambda_i: C:=C\cup\{i\}, U:=U\\{i\}$
        maxnode   = LoL_head->head;
        maxlambda = lambda[maxnode];
        
        vec[maxnode] = CGPT; // set maxnode as coarse node
        lambda[maxnode] = 0;
        --num_left;
        remove_node(&LoL_head, &LoL_tail, maxlambda, maxnode, lists, where);
        col++;
        
        // for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$
        for ( i = ST.IA[maxnode]; i < ST.IA[maxnode+1]; ++i ) {
            
            j = ST.JA[i];
            
            if ( vec[j] != UNPT ) continue; // skip decided variables
            
            vec[j] = FGPT;  // set j as fine node
            remove_node(&LoL_head, &LoL_tail, lambda[j], j, lists, where);
            --num_left;
            
            // Update lambda and linked list after j->F
            for ( l = S->IA[j]; l < S->IA[j+1]; l++ ) {
                k = S->JA[l];
                if ( vec[k] == UNPT ) { // k is unkown
                    remove_node(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
                    newmeas = ++(lambda[k]);
                    enter_list(&LoL_head, &LoL_tail, newmeas, k, lists, where);
                }
            }
            
        } // end for i
        
        // Update lambda and linked list after maxnode->C
        for ( i = S->IA[maxnode]; i < S->IA[maxnode+1]; ++i ) {
            
            j = S->JA[i];
            
            if ( vec[j] != UNPT ) continue; // skip decided variables

            measure = lambda[j];
            remove_node(&LoL_head, &LoL_tail, measure, j, lists, where);
            lambda[j] = --measure;
            
            if ( measure > 0 ) {
                enter_list(&LoL_head, &LoL_tail, measure, j, lists, where);
            }
            else { // j is the only point left, set as fine variable
                vec[j] = FGPT;
                --num_left;
                
                // Update lambda and linked list after j->F
                for ( l = S->IA[j]; l < S->IA[j+1]; l++ ) {
                    k = S->JA[l];
                    if ( vec[k] == UNPT ) { // k is unkown
                        remove_node(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
                        newmeas = ++(lambda[k]);
                        enter_list(&LoL_head, &LoL_tail, newmeas, k, lists, where);
                    }
                } // end for l
            } // end if
            
        } // end for
        
    } // end while
    
    if ( LoL_head ) {
        list_ptr = LoL_head;
        LoL_head->prev_node = NULL;
        LoL_head->next_node = NULL;
        LoL_head = list_ptr->next_node;
        fasp_mem_free(list_ptr);
    }
    
    fasp_icsr_free(&ST);
    fasp_mem_free(work);
    
    return col;
}

/**
 * \fn static INT remove_ff_connections (iCSRmat *S, ivector *vertices,
 *                                       INT row, INT col)
 *
 * \brief Find coarse level variables (C/F splitting): standard
 *
 * \param S            Strong connection matrix
 * \param vertices     Indicator vector for the C/F splitting of the variables
 * \param row          Number of rows of P
 * \param col          Number of columns of P
 *
 * \return Number of cols of P
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 *
 * \note Coarsening Phase TWO: check fine points with F-F connections. Needed by
 *       direct and energy-min interpolations!
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/24/2012: add OMP support
 * Mofified by Chensong Zhang on 05/12/2013: restructure the code
 */
static INT remove_ff_connections (iCSRmat *S,
                                  ivector *vertices,
                                  INT row,
                                  INT col)
{
    // local variables
    INT *vec         = vertices->val;
    INT *graph_array = (INT*)fasp_mem_calloc(row,sizeof(INT));
    INT  set_empty   = TRUE,  C_i_nonempty  = FALSE;
    INT  ci_tilde    = -1,    ci_tilde_mark = -1;
    
    INT  ji, jj, i, j, index;
    INT  myid, mybegin, myend;

    INT nthreads = 1, use_openmp = FALSE;

    if ( use_openmp ) {
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin,myend,i)
#endif
        for ( myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
            for ( i = mybegin; i < myend; ++i ) graph_array[i] = -1;
        }
    }
    else {
        for ( i = 0; i < row; ++i ) graph_array[i] = -1;
    }
    
    for ( i = 0; i < row; ++i ) {
        
        if ( vec[i] != FGPT ) continue; // skip non F-variables
        
        for ( ji = S->IA[i]; ji < S->IA[i+1]; ++ji ) {
            j = S->JA[ji];
            if ( vec[j] == CGPT ) graph_array[j] = i; // mark C-neighbors
        } // end for ji
        
        if ( ci_tilde_mark |= i ) ci_tilde = -1;

        for ( ji = S->IA[i]; ji < S->IA[i+1]; ++ji ) {
            
            j = S->JA[ji];
            
            if ( vec[j] != FGPT ) continue; // skip non F-variables
            
            set_empty = TRUE;
            
            for ( jj = S->IA[j]; jj < S->IA[j+1]; ++jj ) {
                index = S->JA[jj];
                if ( graph_array[index] == i ) {
                    set_empty = FALSE; break; // if there is a C-connections, skip
                }
            } // end for jj
            
            if ( set_empty ) {
                if ( C_i_nonempty ) {
                    vec[i] = CGPT;
                    col++;
                    if ( ci_tilde > -1 ) {
                        vec[ci_tilde] = FGPT;
                        col--;
                        ci_tilde = -1;
                    }
                    C_i_nonempty = FALSE;
                    break;
                }
                else {
                    ci_tilde = j;
                    ci_tilde_mark = i;
                    vec[j] = CGPT;
                    col++;
                    C_i_nonempty = TRUE;
                    i--; // roll back to check i-point again
                    break;
                } // end if C_i_nonempty
            } // end if set_empty
            
        } // end for ji
        
    } // end for i
    
    fasp_mem_free(graph_array);

    return col;
}

/**
 * \fn static INT form_coarse_level_agg (dCSRmat *A, iCSRmat *S, ivector *vertices,
 *                                       INT row, INT aggressive_path)
 *
 * \brief Find coarse level variables (C/F splitting): aggressive
 *
 * \param A                Coefficient matrix, the index starts from zero
 * \param S                Strong connection matrix
 * \param vertices         Indicator vector for the C/F splitting of the variables
 * \param row              Number of rows of P
 * \param aggressive_path  Aggressive path
 *
 * \return Number of cols of P
 *
 * \author Kai Yang, Xiaozhe Hu
 * \date   09/06/2010
 *
 * Modified by Chensong Zhang on 07/05/2012: Fix a data type bug
 * Modified by Chunsheng Feng, Zheng Li on 10/13/2012
 * Modified by Xiaozhe Hu on 04/24/2013: modify aggresive coarsening
 */
static INT form_coarse_level_agg (dCSRmat *A,
                                  iCSRmat *S,
                                  ivector *vertices,
                                  INT row,
                                  INT aggressive_path)
{
    INT col = 0; // initialize col(P): returning output
    INT maxlambda, maxnode, num_left=0;
    INT measure, newmeas;
    INT *vec = vertices->val;
    INT i,j,k,l,m,flag,ci,cj,ck,cl,num_c;
    
    INT *work = (INT*)fasp_mem_calloc(4*row,sizeof(INT));
    INT *lists = work, *where = lists+row, *lambda = where+row;
    INT *cp_index;
    
    ivector CGPT_index, CGPT_rindex;
    LinkList LoL_head = NULL, LoL_tail = NULL, list_ptr = NULL;
    iCSRmat ST,Sh,ShT;
    // Sh is for the strong coupling matrix between temporary CGPTs
    // ShT is the transpose of Sh
    // Snew is for combining the information from S and Sh
        
    fasp_icsr_trans(S, &ST);
    
#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend;
    INT sub_col = 0;
    INT nthreads = FASP_GET_NUM_THREADS();
#endif
    
    vertices->row = A->row;
    
    /************************************************************/
    /* Coarsening Phase ONE: find temporary coarse level points */
    /************************************************************/
    
    num_c = form_coarse_level_std(A, S, vertices, row);
    
    /************************************************************/
    /* Coarsening Phase TWO: find real coarse level points      */
    /************************************************************/
    
    //find Sh, the strong coupling between coarse grid points w.r.t. (path,2)
    if ( aggressive_path < 2 )
        find_strong_couple_agg1(A, S, &Sh, vertices, &CGPT_index, &CGPT_rindex);
    else
        find_strong_couple_agg2(A, S, &Sh, vertices, &CGPT_index, &CGPT_rindex);
    
    fasp_icsr_trans(&Sh, &ShT);
    
    CGPT_index.row  = num_c;
    CGPT_rindex.row = row;
    cp_index        = CGPT_index.val;
    
    // 1. Initialize lambda
#ifdef _OPENMP
#pragma omp parallel for if(num_c>OPENMP_HOLDS)
#endif
    for ( ci = 0; ci < num_c; ++ci ) lambda[ci] = ShT.IA[ci+1]-ShT.IA[ci];
    
    // 2. Set the variables with nonpositvie measure as F-variables
    num_left=0; // number of CGPT in the list LoL
    for (ci=0;ci<num_c;++ci) {
        i=cp_index[ci];
        measure=lambda[ci];
        if (vec[i]!=ISPT) {
            if (measure>0) {
                enter_list(&LoL_head, &LoL_tail, lambda[ci], ci, lists, where);
                num_left++;
            }
            else {
                if ( measure < 0) printf("### WARNING: Negative lambda[%d]!\n", i);
                vec[i]=FGPT; // set i as fine node
                //update the lambda value in the CGPT neighbor of i
                for (ck=Sh.IA[ci];ck<Sh.IA[ci+1];++ck) {
                    cj=Sh.JA[ck];
                    j=cp_index[cj];
                    if (vec[j]!=ISPT) {
                        if (cj<ci) {
                            newmeas=lambda[cj];
                            if (newmeas>0) {
                                remove_node(&LoL_head, &LoL_tail, newmeas, cj, lists, where);
                                num_left--;
                            }
                            newmeas= ++(lambda[cj]);
                            enter_list(&LoL_head, &LoL_tail,  newmeas, cj, lists, where);
                            num_left++;
                        }
                        else{
                            newmeas= ++(lambda[cj]);
                        }//end if cj<ci
                    }//end if vec[j]!=ISPT
                }//end for ck
                
            }//end if
        }
    }
    
    
    // 4. Main loop
    while (num_left>0) {
        // pick $i\in U$ with $\max\lambda_i: C:=C\cup\{i\}, U:=U\\{i\}$
        maxnode=LoL_head->head;
        maxlambda=lambda[maxnode];
        if (maxlambda==0) { printf("### WARNING: Head of the list has measure 0\n");}
        vec[cp_index[maxnode]]=3; // set maxnode as real coarse node, labeled as num 3
        --num_left;
        remove_node(&LoL_head, &LoL_tail, maxlambda, maxnode, lists, where);
        lambda[maxnode]=0;
        col++;//count for the real coarse node after aggressive coarsening
        
        // for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$
        for (ci=ShT.IA[maxnode];ci<ShT.IA[maxnode+1];++ci) {
            cj=ShT.JA[ci];
            j=cp_index[cj];
            
            /* if j is temporary CGPT */
            if (vec[j]==CGPT) {
                vec[j]=4; // set j as 4--fake CGPT
                remove_node(&LoL_head, &LoL_tail, lambda[cj], cj, lists, where);
                --num_left;
                //update the measure for neighboring points
                for (cl=Sh.IA[cj];cl<Sh.IA[cj+1];cl++) {
                    ck=Sh.JA[cl];
                    k=cp_index[ck];
                    if (vec[k]==CGPT) {// k is temporary CGPT
                        remove_node(&LoL_head, &LoL_tail, lambda[ck], ck, lists, where);
                        newmeas= ++(lambda[ck]);
                        enter_list(&LoL_head, &LoL_tail,newmeas, ck, lists, where);
                    }
                }
            } // if
        } // ci
        
        for (ci=Sh.IA[maxnode]; ci<Sh.IA[maxnode+1];++ci) {
            cj=Sh.JA[ci];
            j=cp_index[cj];
            if (vec[j]==CGPT) {// j is temporary CGPT
                measure=lambda[cj];
                remove_node(&LoL_head, &LoL_tail, measure, cj, lists, where);
                measure--;
                lambda[cj]=measure;
                if (measure>0) {
                    enter_list(&LoL_head, &LoL_tail,measure, cj, lists, where);
                }
                else {
                    vec[j]=4; // set j as fake CGPT variable
                    --num_left;
                    for (cl=Sh.IA[cj];cl<Sh.IA[cj+1];cl++) {
                        ck=Sh.JA[cl];
                        k=cp_index[ck];
                        if (vec[k]==CGPT) {// k is temporary CGPT
                            remove_node(&LoL_head, &LoL_tail, lambda[ck], ck, lists, where);
                            newmeas= ++(lambda[ck]);
                            enter_list(&LoL_head, &LoL_tail,newmeas, ck, lists, where);
                        }
                    } // end for l
                } // end if
            } // end if
        } // end for
        
    } // while
    
    // organize the variable type
    // make temporary CGPT--1 and fake CGPT--4 become FGPT
    // make real CGPT--3 to be CGPT
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
    for(i=0;i<row;i++) {
        if(vec[i]==CGPT||vec[i]==4)//change all the temporary or fake CGPT into FGPT
            vec[i]=FGPT;
    }
    
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
    for(i=0;i<row;i++) {
        if(vec[i]==3)// if i is real CGPT
            vec[i]=CGPT;
    }
    
    /************************************************************/
    /* Coarsening Phase THREE: all the FGPTs which have no CGPT */
    /* neighbors within distance 2. Change them into CGPT such  */
    /* the standard interpolation works                         */
    /************************************************************/

#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS) private(myid,mybegin,myend,i,flag,j,k,l,m,sub_col)
    for (myid=0; myid<nthreads; myid++) {
        FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
        for (i=mybegin; i<myend; i++) {
#else
        for (i=0; i<row; i++) {
#endif
            if (vec[i]==FGPT) {
                flag=0; //flag for whether there is any CGPT neighbor within distance of 2
                for (j=S->IA[i]; j<S->IA[i+1]; j++) {
                    k=S->JA[j];
                    if (flag==1) break;
                    else if (vec[k]==CGPT) {
                        flag=1;
                        break;
                    }
                    else if(vec[k]==FGPT) {
                        for (l=S->IA[k]; l<S->IA[k+1]; l++) {
                            m=S->JA[l];
                            if (vec[m]==CGPT) {
                                flag=1;
                                break;
                            }
                        } //end for l
                    }
                } //end for j
                    
                if (flag==0) {
                    vec[i]=CGPT;
#ifdef _OPENMP
                    sub_col++;
#else
                    col++;
#endif
                }
            } //end if
                
        } //end for i
            
#ifdef _OPENMP
#pragma omp critical(col)
        col += sub_col;
    } // myid
#endif
        
    if (LoL_head) {
        list_ptr=LoL_head;
        LoL_head->prev_node=NULL;
        LoL_head->next_node=NULL;
        LoL_head = list_ptr->next_node;
        fasp_mem_free(list_ptr);
    }
        
    fasp_icsr_free(&Sh);
    fasp_icsr_free(&ST);
    fasp_icsr_free(&ShT);
    fasp_mem_free(work);
        
    return col;
}
    
/**
 * \fn static void form_P_pattern_dir (dCSRmat *P, iCSRmat *S, ivector *vertices,
 *                                     INT row, INT col)
 *
 * \brief Generate sparsity pattern of prolongation for direct interpolation
 *
 * \param P         Pointer to the prolongation matrix
 * \param S         Pointer to the set of all strong couplings matrix
 * \param vertices  Pointer to the type of variables
 * \param row       Number of rows of P
 * \param col       Number of cols of P
 *
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue on 05/23/2012: add OMP support
 * Mofified by Chensong Zhang on 05/13/2013: restructure the code
 */
static void form_P_pattern_dir (dCSRmat *P,
                                iCSRmat *S,
                                ivector *vertices,
                                INT row,
                                INT col )
{
    // local variables
    INT i, j, k, index;
    INT *vec = vertices->val;
    
    INT nthreads = 1, use_openmp = FALSE;
    
#ifdef _OPENMP
    if ( row > OPENMP_HOLDS ) {
        use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
    }
#endif
    
    // Initialize P matrix
    P->row = row; P->col = col;
    P->IA  = (INT*)fasp_mem_calloc(row+1, sizeof(INT));
    
    // step 1: Find the structure IA of P first: using P as a counter
    if ( use_openmp ) {
        
        INT mybegin,myend,myid;
#ifdef _OPENMP
#pragma omp parallel for private(myid, mybegin,myend,i,j,k)
#endif
        for ( myid = 0; myid < nthreads; myid++ ) {
            FASP_GET_START_END(myid, nthreads, row, &mybegin, &myend);
            for ( i = mybegin; i < myend; ++i ) {
                switch ( vec[i] ) {
                    case FGPT: // fine grid points
                        for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
                            k = S->JA[j];
                            if ( vec[k] == CGPT ) P->IA[i+1]++;
                        }
                        break;
                        
                    case CGPT: // coarse grid points
                        P->IA[i+1] = 1; break;
                        
                    default: // treat everything else as isolated
                        P->IA[i+1] = 0; break;
                }
            }
        }
        
    }
    
    else {
        
        for ( i = 0; i < row; ++i ) {
            switch ( vec[i] ) {
                case FGPT: // fine grid points
                    for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
                        k = S->JA[j];
                        if ( vec[k] == CGPT ) P->IA[i+1]++;
                    }
                    break;
                     
                case CGPT: // coarse grid points
                    P->IA[i+1] = 1; break;

                default: // treat everything else as isolated
                    P->IA[i+1] = 0; break;
            }
        } // end for i

    } // end if
    
    // Form P->IA from the counter P
    for ( i = 0; i < P->row; ++i ) P->IA[i+1] += P->IA[i];
    P->nnz = P->IA[P->row]-P->IA[0];
    
    // step 2: Find the structure JA of P
    P->JA  = (INT*)fasp_mem_calloc(P->nnz,sizeof(INT));
    P->val = (REAL*)fasp_mem_calloc(P->nnz,sizeof(REAL));
    
    for ( index = i = 0; i < row; ++i ) {
        if ( vec[i] == FGPT ) { // fine node
            for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
                k = S->JA[j];
                if ( vec[k] == CGPT ) P->JA[index++] = k;
            } // end for j
        } // end if
        else if ( vec[i] == CGPT ) { // if node i is on coarse grid
            P->JA[index++] = i;
        }
    }

}

/**
 * \fn static void form_P_pattern_std (dCSRmat *P, iCSRmat *S, ivector *vertices, 
 *                                     INT row, INT col)
 *
 * \brief Generate sparsity pattern of prolongation for standard interpolation
 *
 * \param P         Pointer to the prolongation matrix
 * \param S         Pointer to the set of all strong couplings matrix
 * \param vertices  Pointer to the type of variables
 * \param row       Number of rows of P
 * \param col       Number of cols of P
 *
 * \author Kai Yang, Xiaozhe Hu
 * \date   05/21/2012
 *
 * Modified by Chunsheng Feng, Zheng Li on 10/13/2012: add OMP support
 * Mofified by Chensong Zhang on 05/13/2013: restructure the code
 */
static void form_P_pattern_std (dCSRmat *P,
                                iCSRmat *S,
                                ivector *vertices,
                                INT row,
                                INT col)
{
    // local variables
    INT i, j, k, l, h, index;
    INT *vec = vertices->val;
    
    // number of times a C-point is visited
    INT *times_visited = (INT*)fasp_mem_calloc(row,sizeof(INT));

    P->row = row; P->col = col;
    P->IA  = (INT*)fasp_mem_calloc(row+1, sizeof(INT));
    
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
    for ( i = 0; i < row; i++ ) times_visited[i] = -1;  // initilize

    // step 1: Find the structure IA of P first: use P as a counter
    for ( i = 0; i < row; ++i ) {
        
        if ( vec[i] == FGPT ) { // if node i is a F point
            for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
                
                k = S->JA[j];
                
                // if neighbor of i is a C point, good
                if ( (vec[k] == CGPT) && (times_visited[k] != i) ) {
                    times_visited[k] = i;
                    P->IA[i+1]++;
                }
                
                // if neighbor of i is a F point and k is not i, look for indirect C
                else if ( (vec[k] == FGPT) && (k != i) ) {
                    for ( l = S->IA[k]; l < S->IA[k+1]; l++ ) { // neighbors of k
                        h = S->JA[l];
                        if ( (vec[h] == CGPT) && (times_visited[h] != i) ) {
                            times_visited[h] = i;
                            P->IA[i+1]++;
                        }
                    }  // end for(l=S->IA[k];l<S->IA[k+1];l++)
                }  // end if (vec[k]==CGPT)
                
            } // end for (j=S->IA[i];j<S->IA[i+1];j++)
        }
        
        else if ( vec[i] == CGPT ) { // if node i is a C point
            P->IA[i+1] = 1;
        }
        
        else { // treat everything else as isolated points
            P->IA[i+1] = 0;
        } // end if (vec[i]==FGPT)
        
    } // end for (i=0;i<row;++i)
    
    // Form P->IA from the counter P
    for ( i = 0; i < P->row; ++i ) P->IA[i+1] += P->IA[i];
    P->nnz = P->IA[P->row]-P->IA[0];
    
    // step 2: Find the structure JA of P
    P->JA  = (INT*)fasp_mem_calloc(P->nnz,sizeof(INT));
    P->val = (REAL*)fasp_mem_calloc(P->nnz,sizeof(REAL));
    
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
    for ( i = 0; i < row; ++i ) times_visited[i] = -1;  // reinitilize
    
    for ( i = 0; i < row; ++i ) {
        
        if ( vec[i] == FGPT ) { // if node i is a F point

            index = 0;

            for ( j = S->IA[i]; j < S->IA[i+1]; j++ ) {
                
                k = S->JA[j];
                
                // if neighbor k of i is a C point
                if ( (vec[k] == CGPT) && (times_visited[k] != i) ) {
                    times_visited[k] = i;
                    P->JA[P->IA[i]+index] = k;
                    index++;
                }
                
                // if neighbor k of i is a F point and k is not i
                else if ( (vec[k] == FGPT) && (k != i) ) {
                    for ( l = S->IA[k]; l < S->IA[k+1]; l++ ) { // neighbors of k
                        h = S->JA[l];
                        if ( (vec[h] == CGPT) && (times_visited[h] != i) ) {
                            times_visited[h] = i;
                            P->JA[P->IA[i]+index] = h;
                            index++;
                        }
                        
                    }  // end for (l=S->IA[k];l<S->IA[k+1];l++)
                    
                }  // end if (vec[k]==CGPT)
                
            } // end for (j=S->IA[i];j<S->IA[i+1];j++)
        }
        
        else if ( vec[i] == CGPT ) {
            P->JA[P->IA[i]] = i;
        }
    }
    
    // clea up
    fasp_mem_free(times_visited);
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
