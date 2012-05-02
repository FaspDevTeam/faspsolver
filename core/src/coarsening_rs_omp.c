/*! \file coarsening_rs_omp.c
 *  \brief Coarsening with a modified Ruge-Stuben strategy.
 */

#include "fasp.h"
#include "fasp_functs.h"

static void generate_S_omp (dCSRmat *A, iCSRmat *S, AMG_param *param, int nthreads, int openmp_holds);
static void generate_sparsity_P_omp (dCSRmat *P, iCSRmat *S, ivector *vertices, int row, int col, int nthreads, int openmp_holds);
static int form_coarse_level_omp (dCSRmat *A, iCSRmat *S, ivector *vertices, int row, int nthreads, int openmp_holds);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_amg_coarsening_rs_omp(dCSRmat *A, ivector *vertices, dCSRmat *P,
 *                             AMG_param *param, int nthreads, int openmp_holds)
 * \brief RS coarsening
 *
 * \param A          pointer to the coefficient matrix, the index starts from zero
 * \param vertices   pointer to the indicator ivector of the CF splitting of the vertices
 *                        0: fine gird points
 *                        1: coarse grid points
 *                        2: isolated grid points
 * \param P          pointer to the resulted interpolation matrix (nonzero pattern only)
 * \param param      pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 * \return            SUCCESS or Error message
 *
 * Refter to Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller 
 *           Appendix P475 A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *           Academic Press Inc., San Diego, CA, 2001. 
 * 
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
int fasp_amg_coarsening_rs_omp (dCSRmat *A, 
                                ivector *vertices, 
                                dCSRmat *P, 
                                AMG_param *param, 
                                int nthreads, 
                                int openmp_holds)
{
    int status = SUCCESS;
    
#if FASP_USE_OPENMP
    
    const int   coarsening_type = param->coarsening_type;
    const int   row = A->row;
    const double epsilon_str = param->strong_threshold;
    
    int col = 0;    
    iCSRmat S; // strong n-couplings
    
#if DEBUG_MODE
    printf("coarsening_rs ...... [Start]\n");
#endif
    
#if CHMEM_MODE
    fasp_mem_usage();
#endif
    
#if DEBUG_MODE
    printf("Step 1. form dependent sets ......\n");
#endif
    
    // Step 1: generate S and S transpose
    switch (coarsening_type) {
    case 1: // modified Ruge-Stuben
        generate_S_omp(A, &S, param, nthreads, openmp_holds); break;
    case 3: // compatible relaxation still uses S from modified Ruge-Stuben
        generate_S(A, &S, param); break;
    default: // classical Ruge-Stuben
        generate_S_rs(A, &S, epsilon_str, coarsening_type); break;
    }
    if (S.nnz == 0) return RUN_FAIL;
    
#if DEBUG_MODE
    printf("Step 2. choose C points ......\n");
#endif
    
    // Step 2: standard coarsening algorithm 
    switch (coarsening_type){
    default:
        col = form_coarse_level_omp(A, &S, vertices, row, nthreads,openmp_holds);
        break;
    }
    
#if DEBUG_MODE
    printf("Step 3. find support of C points ......\n");
#endif
    
    // Step 3: generate sparsity pattern of P
    generate_sparsity_P_omp(P, &S, vertices, row, col, nthreads,openmp_holds);
    
#if DEBUG_MODE
    printf("coarsening_rs ...... [Finish]\n");
#endif
    
    fasp_mem_free(S.IA);
    fasp_mem_free(S.JA);
    
#if CHMEM_MODE
    fasp_mem_usage();
#endif
    
#endif
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void generate_S_omp(dCSRmat *A, iCSRmat *S, AMG_param *param, int nthreads)
 * \brief generate the set of all strong couplings S
 * \param A pointer to the coefficient matrix
 * \param S pointer to the set of all strong couplings matrix
 * \param param pointer to AMG parameters
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
static void generate_S_omp(dCSRmat *A, iCSRmat *S, AMG_param *param, int nthreads, int openmp_holds)
{
#if FASP_USE_OPENMP
    double max_row_sum=param->max_row_sum;
    double epsilon_str=param->strong_threshold;
    const int row=A->row, col=A->col;
    const int row_plus_one = row+1;
    const int nnz=A->IA[row]-A->IA[0];
    
    int index, i, j, begin_row, end_row;
    int *ia=A->IA, *ja=A->JA;
    double *aj=A->val;
    
    // get the diagnal entry of A
    dvector diag; fasp_dcsr_getdiag_omp(0, A, &diag,nthreads,openmp_holds);
    
    /** Step 1: generate S */
    double row_scale, row_sum;
    
    // copy the structure of A to S
    S->row=row; S->col=col; S->nnz=nnz; S->val=NULL;
    
    S->IA=(int*)fasp_mem_calloc(row_plus_one, sizeof(INT));
    
#if CHMEM_MODE
    total_alloc_mem += (row+1)*sizeof(double);
#endif
    
    S->JA=(int*)fasp_mem_calloc(nnz, sizeof(INT));
    
#if CHMEM_MODE
    total_alloc_mem += (nnz)*sizeof(INT);
#endif
    
    fasp_iarray_cp_omp(row_plus_one, ia, S->IA, nthreads,openmp_holds);
    fasp_iarray_cp_omp(nnz, ja, S->JA, nthreads,openmp_holds);
    
    if (row > openmp_holds) {
        int mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin, myend, i, row_scale,row_sum,begin_row,end_row,j) ////num_threads(nthreads)
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
                for (i=mybegin; i<myend; i++)
                    {
                        /** compute scaling factor and row sum */
                        row_scale=0; row_sum=0;
    
                        begin_row=ia[i]; end_row=ia[i+1];
                        for (j=begin_row;j<end_row;j++) {
                            row_scale=MIN(row_scale, aj[j]);
                            row_sum+=aj[j];
                        }
                        row_sum=ABS(row_sum)/MAX(SMALLREAL,ABS(diag.val[i]));
    
                        /* compute row entries of S */
                        for (j=begin_row;j<end_row;j++) {
                            if (ja[j]==i) {S->JA[j]=-1; break;}
                        }
    
                        /** if $|\sum_{j=1}^n a_{ij}|> \theta_2 |a_{ii}|$ */
                        if ((row_sum>max_row_sum)&&(max_row_sum<1)) { 
                            /** make all dependencies weak */
                            for (j=begin_row;j<end_row;j++) S->JA[j]=-1;
                        }
                        /** otherwise */
                        else {
                            for (j=begin_row;j<end_row;j++) {
                                /** if $a_{ij}>=\epsilon_{str}*\min a_{ij}$, the connection $a_{ij}$ is set to be weak connection */
                                if (A->val[j]>=epsilon_str*row_scale) S->JA[j]=-1;
                            }
                        }
                    } // end for i
            }
    }
    else {
        for (i=0;i<row;++i) {
            /** compute scaling factor and row sum */
            row_scale=0; row_sum=0;
    
            begin_row=ia[i]; end_row=ia[i+1];
            for (j=begin_row;j<end_row;j++) {
                row_scale=MIN(row_scale, aj[j]);
                row_sum+=aj[j];
            }
            row_sum=ABS(row_sum)/MAX(SMALLREAL,ABS(diag.val[i]));
    
            /* compute row entries of S */
            for (j=begin_row;j<end_row;j++) {
                if (ja[j]==i) {S->JA[j]=-1; break;}
            }
    
            /** if $|\sum_{j=1}^n a_{ij}|> \theta_2 |a_{ii}|$ */
            if ((row_sum>max_row_sum)&&(max_row_sum<1)) { 
                /** make all dependencies weak */
                for (j=begin_row;j<end_row;j++) S->JA[j]=-1;
            }
            /** otherwise */
            else {
                for (j=begin_row;j<end_row;j++) {
                    /** if $a_{ij}>=\epsilon_{str}*\min a_{ij}$, the connection $a_{ij}$ is set to be weak connection */
                    if (A->val[j]>=epsilon_str*row_scale) S->JA[j]=-1; 
                }
            }
        } // end for i
    }
    
    /* Compress the strength matrix */
    index=0;
    for (i=0;i<row;++i) {
        S->IA[i]=index;
        begin_row=ia[i]; end_row=ia[i+1]-1;
        for (j=begin_row;j<=end_row;j++) {
            if (S->JA[j]>-1) {
                S->JA[index]=S->JA[j];
                index++;
            }
        }
    }
    
    if (index > 0) {
        S->IA[row]=index;
        S->nnz=index;
        S->JA=(int*)fasp_mem_realloc(S->JA,index*sizeof(INT));
    }
    else {
        S->nnz = 0;
        S->JA = NULL;
    }
    
    fasp_dvec_free(&diag);
#endif
}

/**
 * \fn static int form_coarse_level_omp(dCSRmat *A, iCSRmat *S, ivector *vertices, int row, int nthreads, int openmp_holds)
 * \brief find coarse level points
 * \param A pointer to the coefficient matrix
 * \param S pointer to the set of all strong couplings matrix
 * \param vertices pointer to the type of vertices (points)
 * \param row integer number of rows of P
 * \param openmp_holds threshold of parallelization
 * \return col integer number of cols of P
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012 Modified by  FENG Chunsheng
 */
static INT form_coarse_level_omp(dCSRmat *A, iCSRmat *S, ivector *vertices, INT row, INT nthreads, INT openmp_holds)
{
    int col = 0; 
#if FASP_USE_OPENMP
    unsigned int maxlambda, maxnode, num_left=0;
    double measure, new_meas;
    
    int *ia=A->IA, *vec=vertices->val;
    int ci_tilde = -1, ci_tilde_mark = -1;
    int set_empty = 1, C_i_nonempty = 0;
    int ji,jj,i,j,k,l,index;
    int myid;
    int mybegin;
    int myend;
    int stride_i;
    
    int *work = (int*)fasp_mem_calloc(4*row,sizeof(INT));
    int *lists = work, *where = lists+row, *lambda = where+row, *graph_array = lambda+row;
    
    LinkList LoL_head = NULL, LoL_tail = NULL, list_ptr = NULL;    
    
    iCSRmat ST; fasp_icsr_trans(S, &ST);
    
    /**************************************************/
    /* Coarsening Phase ONE: find coarse level points */
    /**************************************************/    
    
    // 1. Initialize lambda
    if (row > openmp_holds) {
#pragma omp for parallel private(myid, mybegin,myend,i) 
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
                for (i=mybegin; i<myend; i++) lambda[i]=ST.IA[i+1]-ST.IA[i];
            }
    }
    else {
        for (i=0;i<row;++i) lambda[i]=ST.IA[i+1]-ST.IA[i];
    }
    
    // 2. Before the following algorithm starts, filter out the variables which
    // have no connections at all and assign special F-variables.
    if (row > openmp_holds) {
#pragma omp parallel for reduction(+:num_left) private(myid, mybegin, myend, i) 
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
                for (i=mybegin; i<myend; i++)
                    {
                        if ( (ia[i+1]-ia[i])<=1 ) {
                            vec[i]=ISPT; // set i as an ISOLATED fine node
                            lambda[i]=0;
                        }
                        else {
                            vec[i]=UNPT; // set i as a undecided node
                            num_left++;
                        }
                    }
            }
    }
    else {
        for (i=0;i<row;++i)
            {
                if ( (ia[i+1]-ia[i])<=1 ) {
                    vec[i]=ISPT; // set i as an ISOLATED fine node
                    lambda[i]=0;
                }
                else {
                    vec[i]=UNPT; // set i as a undecided node
                    num_left++;
                }
            }
    }
    
    // 3. Set the variables with nonpositvie measure as F-variables
    for (i=0;i<row;++i)
        {
            measure=lambda[i];
            if (vec[i]!=ISPT) {
                if (measure>0) {
                    enter_list(&LoL_head, &LoL_tail, lambda[i], i, lists, where);
                }
                else {
                    if (measure<0) printf("Warning: negative lambda!\n");
                    vec[i]=FGPT; // set i as fine node
                    for (k=S->IA[i];k<S->IA[i+1];++k) {
                        j=S->JA[k];
                        if (vec[j]!=ISPT) {
                            if (j<i) {
                                new_meas=lambda[j];
                                if (new_meas>0) {
                                    remove_point(&LoL_head, &LoL_tail, new_meas, j, lists, where);
                                }
                                new_meas= ++(lambda[j]);
                                enter_list(&LoL_head, &LoL_tail,  new_meas, j, lists, where);
                            }
                            else {
                                new_meas= ++(lambda[j]);
                            }
                        }
                    }
                    --num_left;
                }
            }
        }
    
    // 4. Main loop
    while (num_left>0) {
    
        // pick $i\in U$ with $\max\lambda_i: C:=C\cup\{i\}, U:=U\\{i\}$
        maxnode=LoL_head->head;
        maxlambda=lambda[maxnode];
    
        vec[maxnode]=CGPT; // set maxnode as coarse node
        lambda[maxnode]=0;
        --num_left;
        remove_point(&LoL_head, &LoL_tail, maxlambda, maxnode, lists, where);
        col++;
    
        // for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$
        for (i=ST.IA[maxnode];i<ST.IA[maxnode+1];++i) {
            j=ST.JA[i];
    
            /** if j is unkown */
            if (vec[j]==UNPT) {
                vec[j]=FGPT;  // set j as fine node
                remove_point(&LoL_head, &LoL_tail, lambda[j], j, lists, where);
                --num_left;
    
                for (l=S->IA[j];l<S->IA[j+1];l++) {
                    k=S->JA[l];
                    if (vec[k]==UNPT) // k is unkown
                        {
                            remove_point(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
                            new_meas= ++(lambda[k]);
                            enter_list(&LoL_head, &LoL_tail,new_meas, k, lists, where);
                        }
                }
    
            } // if
        } // i
    
        for (i=S->IA[maxnode];i<S->IA[maxnode+1];++i) {
            j=S->JA[i];
            if (vec[j]==UNPT) // j is unkown
                {
                    measure=lambda[j];
                    remove_point(&LoL_head, &LoL_tail,measure, j, lists, where);
                    lambda[j]=--measure;
                    if (measure>0) {
                        enter_list(&LoL_head, &LoL_tail,measure, j, lists, where);
                    }
                    else {
                        vec[j]=FGPT; // set j as fine variable
                        --num_left;
    
                        for (l=S->IA[j];l<S->IA[j+1];l++) {
                            k=S->JA[l];
                            if (vec[k]==UNPT) // k is unkown
                                {
                                    remove_point(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
                                    new_meas= ++(lambda[k]);
                                    enter_list(&LoL_head, &LoL_tail,new_meas, k, lists, where);
                                }
                        } // end for l
                    } // end if
                } // end if
        } // end for
    } // while
    
    if (LoL_head) {
        list_ptr=LoL_head;
        LoL_head->prev_elt=NULL;
        LoL_head->next_elt=NULL;
        LoL_head = list_ptr->next_elt;
        fasp_mem_free(list_ptr);
    }
    
    /****************************************************************/
    /* Coarsening Phase TWO: check fine points for coarse neighbors */
    /****************************************************************/
    if (row > openmp_holds) {
#pragma omp for parallel private(myid, mybegin,myend,i) 
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
                for (i=mybegin; i<myend; ++i) graph_array[i] = -1;
            }
    }
    else {
        for (i=0; i < row; ++i) graph_array[i] = -1;
    }
    
    for (i=0; i < row; ++i)
        {
            if (ci_tilde_mark |= i) ci_tilde = -1;
    
            if (vec[i] == FGPT) // F-variable
                {
                    for (ji = S->IA[i]; ji < S->IA[i+1]; ++ji) {
                        j = S->JA[ji];
                        if (vec[j] == CGPT) // C-variable
                            graph_array[j] = i;
                    } // end for ji
    
                    for (ji = S->IA[i]; ji < S->IA[i+1]; ++ji) {
    
                        j = S->JA[ji];
    
                        if (vec[j] == FGPT) {
                            set_empty = 1;
    
                            for (jj = S->IA[j]; jj < S->IA[j+1]; ++jj) {
                                index = S->JA[jj];
                                if (graph_array[index] == i) {
                                    set_empty = 0;
                                    break;
                                }
                            } // end for jj
    
                            if (set_empty) {
                                if (C_i_nonempty) {
                                    vec[i] = CGPT;
                                    col++;
                                    if (ci_tilde > -1) {
                                        vec[ci_tilde] = FGPT;
                                        col--;
                                        ci_tilde = -1;
                                    }
                                    C_i_nonempty = 0;
                                    break;
                                }
                                else {
                                    ci_tilde = j;
                                    ci_tilde_mark = i;
                                    vec[j] = CGPT;
                                    col++;
                                    C_i_nonempty = 1;
                                    i--;
                                    break;
                                } // end if C_i_nonempty
                            } // end if set_empty
    
                        }
                    }
                }
        }
    
    fasp_icsr_free(&ST);
    
    fasp_mem_free(work);
    
#endif
    return col;
}

/**
 * \fn static void generate_sparsity_P_omp(dCSRmat *P, iCSRmat *S, ivector *vertices, int row, int col, int nthreads, int openmp_holds)
 * \brief find coarse level points
 * \param P pointer to the prolongation matrix
 * \param S pointer to the set of all strong couplings matrix
 * \param vertices pointer to the type of vertices (points)
 * \param row integer number of rows of P
 * \param col integer number of cols of P
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author FENG Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 * \date Jan/11/2012 Modified by FENG Chunsheng
 */
static void generate_sparsity_P_omp(dCSRmat *P, iCSRmat *S, ivector *vertices, INT row, INT col, INT nthreads, INT openmp_holds)
{
#if FASP_USE_OPENMP
    int i,j,k,index=0;
    int *vec=vertices->val;
    
    P->row=row; P->col=col;    
    P->IA=(int*)fasp_mem_calloc(row+1, sizeof(INT));    
    
#if CHMEM_MODE
    total_alloc_mem += (row+1)*sizeof(INT);
#endif
    
    // step 1: Find the structure IA of P first
    if (row > openmp_holds) {
        int mybegin,myend,myid;
#pragma omp for parallel private(myid, mybegin,myend,i,j,k) 
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
                for (i=mybegin; i<myend; ++i)
                    {
                        if (vec[i]==FGPT) // if node i is on fine grid
                            {
                                for (j=S->IA[i];j<S->IA[i+1];j++)
                                    {
                                        k=S->JA[j];
                                        if (vec[k]==CGPT)
                                            {
                                                P->IA[i+1]++;
                                            }
                                    }
                            }
                        else if (vec[i]==ISPT) { // if node i is a special fine node
                            P->IA[i+1]=0;
                        }
                        else { // if node i is on coarse grid 
                            P->IA[i+1]=1;
                        }
                    }
            }
    }
    else {
        for (i=0;i<row;++i)
            {
                if (vec[i]==FGPT) // if node i is on fine grid
                    {
                        for (j=S->IA[i];j<S->IA[i+1];j++)
                            {
                                k=S->JA[j];
                                if (vec[k]==CGPT)
                                    {
                                        P->IA[i+1]++;
                                    }
                            }
                    }
                else if (vec[i]==ISPT) { // if node i is a special fine node
                    P->IA[i+1]=0;
                }
                else { // if node i is on coarse grid 
                    P->IA[i+1]=1;
                }
            }
    }
    
    for (i=0;i<P->row;++i) P->IA[i+1]+=P->IA[i];
    
    P->nnz=P->IA[P->row]-P->IA[0];
    
    // step 2: Find the structure JA of P
    P->JA=(int*)fasp_mem_calloc(P->nnz,sizeof(INT));    
    P->val=(double*)fasp_mem_calloc(P->nnz,sizeof(double));
    
#if CHMEM_MODE
    total_alloc_mem += (P->nnz)*sizeof(INT);
    total_alloc_mem += (P->nnz)*sizeof(double);
#endif
    
    for (i=0;i<row;++i)
        {
            if (vec[i]==FGPT) // fine node
                {
                    for (j=S->IA[i];j<S->IA[i+1];j++) {
                        k=S->JA[j];
                        if (vec[k]==CGPT) {
                            P->JA[index]=k;
                            index++;
                        } // end if 
                    } // end for j
                } // end if
            else if (vec[i]==ISPT) 
                {
                    // if node i is a special fine node, do nothing
                }
            else // if node i is on coarse grid 
                {
                    P->JA[index]=i;
                    index++;
                }
        }
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
