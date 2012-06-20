/*! \file coarsening_rs.c
 *  \brief Coarsening with a modified Ruge-Stuben strategy.
 */

#include <omp.h>
#include "fasp.h"
#include "fasp_functs.h"

// These linked list operations are adapted from hypre 2.0
static LinkList create_node (INT Item);
static void dispose_node (LinkList node_ptr);
void enter_list (LinkList *LoL_head_ptr, LinkList *LoL_tail_ptr, INT measure, INT index, INT *lists, INT *where);
void remove_node (LinkList *LoL_head_ptr, LinkList *LoL_tail_ptr, INT measure, INT index, INT *lists, INT *where);

// Private routines for RS coarsening
static INT form_coarse_level (dCSRmat *A, iCSRmat *S, ivector *vertices, INT row, INT interpolation_type);
static INT form_coarse_level_ag (dCSRmat *A, iCSRmat *S, ivector *vertices, INT row,INT interp_type);
static void generate_S (dCSRmat *A, iCSRmat *S, AMG_param *param);
void generate_S_rs (dCSRmat *A, iCSRmat *S, REAL epsilon_str, INT coarsening_type);
static void generate_sparsity_P(dCSRmat *P, iCSRmat *S, ivector *vertices, INT row, INT col);
static void generate_sparsity_P_standard (dCSRmat *P, iCSRmat *S, ivector *vertices, INT row, INT col);
static void generate_S_rs_ag (dCSRmat *A, iCSRmat *S, iCSRmat *Sh, ivector *vertices, ivector *CGPT_index, ivector *CGPT_rindex);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_amg_coarsening_rs (dCSRmat *A, ivector *vertices, dCSRmat *P, iCSRmat *S, 
 *                                 AMG_param *param)
 *
 * \brief RS coarsening
 *
 * \param A          Pointer to the coefficient matrix, the index starts from zero
 * \param vertices   Pointer to the indicator ivector of the CF splitting of the vertices
 *                        0: fine gird points
 *                        1: coarse grid points
 *                        2: isolated grid points
 * \param P          Pointer to the resulted INTerpolation matrix (nonzero pattern only)
 * \param S          Pointer to strength matrix
 * \param param      Pointer to AMG parameters
 *
 * \return           SUCCESS or Error message
 *
 * \note Ref Multigrid by U. Trottenberg, C. W. Oosterlee and A. Schuller 
 *           Appendix P475 A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *           Academic Press Inc., San Diego, CA, 2001. 
 * 
 * \author Xuehai Huang, Chensong Zhang, Xiaozhe Hu, Ludmil Zikatanov
 * \date   09/06/2010
 *
 * Modified by Xiaozhe Hu: add strength matrix as an input/output. -- 05/23/201
 *
 */

INT fasp_amg_coarsening_rs (dCSRmat *A, 
                            ivector *vertices, 
                            dCSRmat *P, 
							iCSRmat *S,
                            AMG_param *param)
{    
    const INT   coarsening_type = param->coarsening_type;
	INT   interpolation_type = param->interpolation_type;
    const INT   row = A->row;
    const REAL epsilon_str = param->strong_threshold;
    INT status = SUCCESS, col = 0;    
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_coarsening_rs ...... [Start]\n");
#endif
    
#if DEBUG_MODE
    printf("### DEBUG: Step 1. form dependent sets ......\n");
#endif
   
	// make sure standard interpolaiton is used for aggressive coarsening
	if (coarsening_type == 4) interpolation_type = INTERP_STD;

    // Step 1: generate S and S transpose
    switch (coarsening_type) {
    case 1: // modified Ruge-Stuben
        generate_S(A, S, param); break;
    case 3: // compatible relaxation still uses S from modified Ruge-Stuben
        generate_S(A, S, param); break;
	case 4: // aggressive coarsening
		generate_S(A, S, param); break;
    default: // classical Ruge-Stuben
        generate_S_rs(A, S, epsilon_str, coarsening_type); break;
    }

    if (S->nnz == 0) return RUN_FAIL;
    
#if DEBUG_MODE
    printf("### DEBUG: Step 2. choose C points ......\n");
#endif
    
    // Step 2: standard coarsening algorithm 
    switch (coarsening_type) {
    case 3: // compatible relaxation (Need to be modified --Chensong)
        col = fasp_amg_coarsening_cr(0,A->row-1, A, vertices, param);    
        break;
	case 4: // aggressive coarsening
		col=form_coarse_level_ag(A, S, vertices,row,interpolation_type);
		break;
    default:
        col = form_coarse_level(A, S, vertices, row, interpolation_type);
        break;
    }
    
#if DEBUG_MODE
    printf("### DEBUG: Step 3. find support of C points ......\n");
#endif
    
    // Step 3: generate sparsity pattern of P
	switch (interpolation_type){
		case INTERP_STD: // standard interpolaiton
			generate_sparsity_P_standard (P, S, vertices, row, col);
			break;
		default: // direct interpolation & energy minimization interpolaiton
			generate_sparsity_P(P, S, vertices, row, col);
			break;
	}
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_amg_coarsening_rs ...... [Finish]\n");
#endif
    
    
#if CHMEM_MODE
    fasp_mem_usage();
#endif
    
    return status;
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

#define LIST_HEAD -1 /**< head of the linked list */
#define LIST_TAIL -2 /**< tail of the linked list */

/**
 * \fn static void dispose_node( LinkList node_ptr )
 *
 * \brief Free memory space used by node_ptr
 *
 * \param node_ptr   Pointer to the node in linked list
 * 
 * \author Xuehai Huang
 * \date   09/06/2009
 */
static void dispose_node (LinkList node_ptr)
{
    if (node_ptr) fasp_mem_free( node_ptr );
}

/**
 * \fn static void remove_node (LinkList *LoL_head_ptr, LinkList *LoL_tail_ptr,
 *                              INT measure, INT index, INT *lists, INT *where)
 *
 * \brief Removes a point from the lists
 * 
 * \author Xuehai Huang
 * \date   09/06/2009
 */
void remove_node (LinkList *LoL_head_ptr, 
                  LinkList *LoL_tail_ptr,
                  INT measure, 
                  INT index, 
                  INT *lists, 
                  INT *where)
{
    LinkList LoL_head = *LoL_head_ptr;
    LinkList LoL_tail = *LoL_tail_ptr;
    LinkList list_ptr = LoL_head;
    
    do {
        if (measure == list_ptr->data) {
            /* point to be removed is only point on list,
               which must be destroyed */
            if (list_ptr->head == index && list_ptr->tail == index) {
                /* removing only list, so num_left better be 0! */
                if (list_ptr == LoL_head && list_ptr == LoL_tail) {
                    LoL_head = NULL;
                    LoL_tail = NULL;
                    dispose_node(list_ptr);
    
                    *LoL_head_ptr = LoL_head;
                    *LoL_tail_ptr = LoL_tail;
                    return;
                }
                else if (LoL_head == list_ptr) { /*removing 1st (max_measure) list */
                    list_ptr -> next_node -> prev_node = NULL;
                    LoL_head = list_ptr->next_node;
                    dispose_node(list_ptr);
    
                    *LoL_head_ptr = LoL_head;
                    *LoL_tail_ptr = LoL_tail;
                    return;
                }
                else if (LoL_tail == list_ptr) { /* removing last list */
                    list_ptr -> prev_node -> next_node = NULL;
                    LoL_tail = list_ptr->prev_node;
                    dispose_node(list_ptr);
    
                    *LoL_head_ptr = LoL_head;
                    *LoL_tail_ptr = LoL_tail;
                    return;
                }
                else {
                    list_ptr -> next_node -> prev_node = list_ptr -> prev_node;
                    list_ptr -> prev_node -> next_node = list_ptr -> next_node;
                    dispose_node(list_ptr);
    
                    *LoL_head_ptr = LoL_head;
                    *LoL_tail_ptr = LoL_tail;
                    return;
                }
            }
            else if (list_ptr->head == index) { /* index is head of list */
                list_ptr->head = lists[index];
                where[lists[index]] = LIST_HEAD;
                return;
            }
            else if (list_ptr->tail == index) { /* index is tail of list */
                list_ptr->tail = where[index];
                lists[where[index]] = LIST_TAIL;
                return;
            }
            else { /* index is in middle of list */
                lists[where[index]] = lists[index];
                where[lists[index]] = where[index];
                return;
            }
        }
        list_ptr = list_ptr -> next_node;
    } while (list_ptr != NULL);
    
    printf("No such list!\n");
    return;
}

/**
 * \fn static LinkList create_node (INT Item)
 *
 * \brief Create an node using Item for its data field
 * 
 * \author Xuehai Huang
 * \date   09/06/2009
 */
static LinkList create_node (INT Item)
{
    LinkList new_node_ptr;
    
    /* Allocate memory space for the new node. 
     * return with error if no space available
     */    
    new_node_ptr = (LinkList) fasp_mem_calloc(1,sizeof(ListElement));
    
    new_node_ptr -> data = Item;
    new_node_ptr -> next_node = NULL;
    new_node_ptr -> prev_node = NULL;
    new_node_ptr -> head = LIST_TAIL;
    new_node_ptr -> tail = LIST_HEAD;
    
    return (new_node_ptr);
}

/**
 * \fn static void enter_list (LinkList *LoL_head_ptr, LinkList *LoL_tail_ptr,
 *                             INT measure, INT index, INT *lists, INT *where)
 *
 * \brief Places point in new list
 * 
 * \author Xuehai Huang
 * \date   09/06/2009
 */
void enter_list (LinkList *LoL_head_ptr, 
                 LinkList *LoL_tail_ptr,
                 INT measure, 
                 INT index, 
                 INT *lists, 
                 INT *where)
{
    LinkList   LoL_head = *LoL_head_ptr;
    LinkList   LoL_tail = *LoL_tail_ptr;    
    LinkList   list_ptr;
    LinkList   new_ptr;
    
    INT        old_tail;
    
    list_ptr =  LoL_head;
    
    if (LoL_head == NULL) { /* no lists exist yet */
        new_ptr = create_node(measure);
        new_ptr->head = index;
        new_ptr->tail = index;
        lists[index] = LIST_TAIL;
        where[index] = LIST_HEAD; 
        LoL_head = new_ptr;
        LoL_tail = new_ptr;
    
        *LoL_head_ptr = LoL_head;
        *LoL_tail_ptr = LoL_tail;
        return;
    }
    else {
        do {
            if (measure > list_ptr->data) {
                new_ptr = create_node(measure);
    
                new_ptr->head = index;
                new_ptr->tail = index;
    
                lists[index] = LIST_TAIL;
                where[index] = LIST_HEAD;
    
                if ( list_ptr->prev_node != NULL) { 
                    new_ptr->prev_node            = list_ptr->prev_node;
                    list_ptr->prev_node->next_node = new_ptr;   
                    list_ptr->prev_node           = new_ptr;
                    new_ptr->next_node            = list_ptr;
                }
                else {
                    new_ptr->next_node  = list_ptr;
                    list_ptr->prev_node = new_ptr;
                    new_ptr->prev_node  = NULL;
                    LoL_head = new_ptr;
                }
    
                *LoL_head_ptr = LoL_head;
                *LoL_tail_ptr = LoL_tail; 
                return;
            }
            else if (measure == list_ptr->data) {
                old_tail = list_ptr->tail;
                lists[old_tail] = index;
                where[index] = old_tail;
                lists[index] = LIST_TAIL;
                list_ptr->tail = index;
                return;
            }
            
            list_ptr = list_ptr->next_node;
        } while (list_ptr != NULL);
    
        new_ptr = create_node(measure);   
        new_ptr->head = index;
        new_ptr->tail = index;
        lists[index] = LIST_TAIL;
        where[index] = LIST_HEAD;
        LoL_tail->next_node = new_ptr;
        new_ptr->prev_node = LoL_tail;
        new_ptr->next_node = NULL;
        LoL_tail = new_ptr;
    
        *LoL_head_ptr = LoL_head;
        *LoL_tail_ptr = LoL_tail;
        return;
    }
}

/**
 * \fn static void generate_S (dCSRmat *A, iCSRmat *S, AMG_param *param)
 *
 * \brief Generate the set of all strong couplings S
 *
 * \param A      Pointer to the coefficient matrix
 * \param S      Pointer to the set of all strong couplings matrix
 * \param param  Pointer to AMG parameters
 * 
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/25/2012
 */

static void generate_S ( dCSRmat *A, 
		                 iCSRmat *S, 
						 AMG_param *param )
{
    REAL max_row_sum=param->max_row_sum;
    REAL epsilon_str=param->strong_threshold;
    const INT row=A->row, col=A->col;
    const INT row_plus_one = row+1;
    const INT nnz=A->IA[row]-A->IA[0];
    
    INT index, i, j, begin_row, end_row;
    INT *ia=A->IA, *ja=A->JA;
    REAL *aj=A->val;
    
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || row <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = use_openmp;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    // get the diagnal entry of A
    dvector diag; fasp_dcsr_getdiag(0, A, &diag);
    
    /* Step 1: generate S */
    REAL row_scale, row_sum;
    
    // copy the structure of A to S
    S->row=row; S->col=col; S->nnz=nnz; S->val=NULL;
    
    S->IA=(INT*)fasp_mem_calloc(row_plus_one, sizeof(INT));
    
#if CHMEM_MODE
    total_alloc_mem += (row+1)*sizeof(REAL);
#endif
    
    S->JA=(INT*)fasp_mem_calloc(nnz, sizeof(INT));
    
#if CHMEM_MODE
    total_alloc_mem += (nnz)*sizeof(INT);
#endif
    
    fasp_iarray_cp(row_plus_one, ia, S->IA);
    fasp_iarray_cp(nnz, ja, S->JA);
    
    if (use_openmp) {
        INT mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin,myend,i,row_scale,row_sum,begin_row,end_row,j)
        for (myid = 0; myid < nthreads; myid++ )
            {
                FASP_GET_START_END(myid, nthreads, row, mybegin, myend);
                for (i=mybegin; i<myend; i++)
                    {
                        /* compute scaling factor and row sum */
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
    
                        /* if $|\sum_{j=1}^n a_{ij}|> \theta_2 |a_{ii}|$ */
                        if ((row_sum>max_row_sum)&&(max_row_sum<1)) { 
                            /* make all dependencies weak */
                            for (j=begin_row;j<end_row;j++) S->JA[j]=-1;
                        }
                        /* otherwise */
                        else {
                            for (j=begin_row;j<end_row;j++) {
                                /* if $a_{ij}>=\epsilon_{str}*\min a_{ij}$, the connection $a_{ij}$ 
                                   is set to be weak connection */
                                if (A->val[j]>=epsilon_str*row_scale) S->JA[j]=-1;
                            }
                        }
                    } // end for i
            }
    }
    else {
        for (i=0;i<row;++i) {
            /* compute scaling factor and row sum */
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
    
            /* if $|\sum_{j=1}^n a_{ij}|> \theta_2 |a_{ii}|$ */
            if ((row_sum>max_row_sum)&&(max_row_sum<1)) { 
                /* make all dependencies weak */
                for (j=begin_row;j<end_row;j++) S->JA[j]=-1;
            }
            /* otherwise */
            else {
                for (j=begin_row;j<end_row;j++) {
                    /* if $a_{ij}>=\epsilon_{str}*\min a_{ij}$, the connection $a_{ij}$ is set to be weak connection */
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
        S->JA=(INT*)fasp_mem_realloc(S->JA,index*sizeof(INT));
    }
    else {
        S->nnz = 0;
        S->JA = NULL;
    }
    
    fasp_dvec_free(&diag);
}

/**
 * \fn static void generate_S_rs (dCSRmat *A, iCSRmat *S, REAL epsilon_str, INT coarsening_type)
 *
 * \brief Generate the set of all strong negative or strong absolute couplings
 *
 * \param A                Pointer to the coefficient matrix
 * \param S                Pointer to the set of all strong couplings matrix
 * \param epsilon_str      Strong coupled ratio
 * \param coarsening_type  Coarsening type (2: strong negative couplings, 3: strong absolute couplings)
 * 
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 */
void generate_S_rs (dCSRmat *A, 
                    iCSRmat *S, 
                    REAL epsilon_str, 
                    INT coarsening_type)
{
    INT i,j;
    REAL *amax=(REAL *)fasp_mem_calloc(A->row,sizeof(REAL));
    
    // get the maximum absolute negative value data in each row of A
    if (coarsening_type==2) { // original RS coarsening, just consider negative strong coupled

        for (i=0;i<A->row;++i) {
            amax[i]=0;
            for (j=A->IA[i];j<A->IA[i+1];j++) {
                if ((-A->val[j])>amax[i] && A->JA[j]!=i) {
                    amax[i]=-A->val[j];
                }
            }
        }
    }
    
    if (coarsening_type==3) { // original RS coarsening, consider absolute strong coupled

        for (i=0;i<A->row;++i) {
            amax[i]=0;
            for (j=A->IA[i];j<A->IA[i+1];j++) {
                if (ABS(A->val[j])>amax[i] && A->JA[j]!=i) {
                    amax[i]=ABS(A->val[j]);
                }
            }
        }
    }
    
    // step 1: Find first the structure IA of S
    S->row=A->row;
    S->col=A->col;
    S->val=NULL;
    S->JA=NULL;
    S->IA=(INT*)fasp_mem_calloc(S->row+1, sizeof(INT));
    
    if (coarsening_type==2) {
        for (i=0;i<S->row;++i) {
            for (j=A->IA[i];j<A->IA[i+1];j++) {
                if ((-A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i) {
                    S->IA[i+1]++;
                }
            }    
        }
    }
    
    if (coarsening_type==3) {
        for (i=0;i<S->row;++i) {
            for (j=A->IA[i];j<A->IA[i+1];j++) {
                if (ABS(A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i) {
                    S->IA[i+1]++;
                }
            }
        }
    }
    
    for (i=0;i<S->row;++i) S->IA[i+1]+=S->IA[i];
    
    // step 2: Find the structure JA of S
    INT index=0;
    S->JA=(INT*)fasp_mem_calloc(S->IA[S->row],sizeof(INT));
    
    if (coarsening_type==2) {
        for (i=0;i<S->row;++i) {
            for (j=A->IA[i];j<A->IA[i+1];j++) {
                if ((-A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i) {
                    S->JA[index]= A->JA[j];
                    index++;
                }
            }
        }
    }
    
    if (coarsening_type==3) {
        for (i=0;i<S->row;++i) {
            for (j=A->IA[i];j<A->IA[i+1];j++) {
                if (ABS(A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i) {
                    S->JA[index]= A->JA[j];
                    index++;
                }
            }
        }
    }
    
    fasp_mem_free(amax);
}

/**
 * \fn static void generate_S_rs_ag (dCSRmat *A, iCSRmat *S, iCSRmat *Sh, ivector *vertices, ivector *CGPT_index, ivector *CGPT_rindex)
 *
 * \brief Generate the set of all strong negative or strong absolute couplings, aggressive coarsening A1 or A2 
 *
 * \param A                Pointer to the coefficient matrix
 * \param S                Pointer to the set of all strong couplings matrix
 * \param Sh               Pointer to the set of all strong couplings matrix between coarse grid points 
 * \param vertices         Pointer to the type of variables--C/F splitting
 * \param CGPT_index       Pointer to the index of CGPT from the list of CGPT to the list of all points
 * \param CGPT_rindex      Pointer to the index of CGPT from the list of all points to the list of CGPT
 * 
 * \author Kai Yang, Xiaozhe Hu
 * \date   09/06/2010
 */
static void generate_S_rs_ag (dCSRmat *A, 
                              iCSRmat *S, 
                              iCSRmat *Sh, 
                              ivector *vertices,
                              ivector * CGPT_index,
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
	
    for (i=0; i<num_c; i++) times_visited[i]=0;
	
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
                        if(cp_rindex[ck]>=num_c) printf("generate_S_rs_ag: index exceed bound!\n");
                        
                        cck=cp_rindex[ck];
                        if (times_visited[cck]!=ci+1) {
                            //newly visited
                            times_visited[cck]=ci+1;//marked as strongly connected from ci
                            count++;
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
                        
                        if(times_visited[cck]!=ci+1)
                        {
                            //newly visited
                            times_visited[cck]=ci+1;
                            Sh->JA[count]=cck;
                            count++;
                        }
                        
                    }//end if
                }//end for k
                
            }//end if
        }//end for j 
         
	    if(count!=Sh->IA[ci+1]) printf("generate_S_rs_ag: inconsistency in number of nonzeros values\n ");
	    
    }//end for ci
    
	fasp_mem_free(times_visited);
}


/**
 * \fn static INT form_coarse_level (dCSRmat *A, iCSRmat *S, ivector *vertices, INT row)
 *
 * \brief Find coarse level points
 *
 * \param A         Pointer to the coefficient matrix
 * \param S         Pointer to the set of all strong couplings matrix
 * \param vertices  Pointer to the type of variables
 * \param row       Number of rows of P
 *
 * \return Number of cols of P
 * 
 * \author Xuehai Huang, Chensong Zhang
 * \date   09/06/2010
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/24/2012
 */

static INT form_coarse_level ( dCSRmat *A, 
		                       iCSRmat *S, 
							   ivector *vertices, 
							   INT row,
							   INT interpolation_type)
{
    INT col = 0; 
    unsigned INT maxlambda, maxnode, num_left=0;
    REAL measure, new_meas;
    
    INT *ia=A->IA, *vec=vertices->val;
    INT ci_tilde = -1, ci_tilde_mark = -1;
    INT set_empty = 1, C_i_nonempty = 0;
    INT ji,jj,i,j,k,l,index;
    INT myid;
    INT mybegin;
    INT myend;
    INT stride_i;
    
    INT *work = (INT*)fasp_mem_calloc(4*row,sizeof(INT));
    INT *lists = work, *where = lists+row, *lambda = where+row, *graph_array = lambda+row;
    
    LinkList LoL_head = NULL, LoL_tail = NULL, list_ptr = NULL;    
    
    iCSRmat ST; fasp_icsr_trans(S, &ST);

	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || row <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
    /**************************************************/
    /* Coarsening Phase ONE: find coarse level points */
    /**************************************************/    
    
    // 1. Initialize lambda
    if (use_openmp) {
#pragma omp parallel for private(myid, mybegin,myend,i) 
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
    if (use_openmp) {
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
                                    remove_node(&LoL_head, &LoL_tail, new_meas, j, lists, where);
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
        remove_node(&LoL_head, &LoL_tail, maxlambda, maxnode, lists, where);
        col++;
    
        // for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$
        for (i=ST.IA[maxnode];i<ST.IA[maxnode+1];++i) {
            j=ST.JA[i];
    
            /* if j is unkown */
            if (vec[j]==UNPT) {
                vec[j]=FGPT;  // set j as fine node
                remove_node(&LoL_head, &LoL_tail, lambda[j], j, lists, where);
                --num_left;
    
                for (l=S->IA[j];l<S->IA[j+1];l++) {
                    k=S->JA[l];
                    if (vec[k]==UNPT) // k is unkown
                        {
                            remove_node(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
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
                    remove_node(&LoL_head, &LoL_tail,measure, j, lists, where);
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
                                    remove_node(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
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
        LoL_head->prev_node=NULL;
        LoL_head->next_node=NULL;
        LoL_head = list_ptr->next_node;
	
        fasp_mem_free(list_ptr);
    }
    
	if (interpolation_type != INTERP_STD) {
    /****************************************************************/
    /* Coarsening Phase TWO: check fine points for coarse neighbors */
    /****************************************************************/
    if (use_openmp) {
#pragma omp parallel for private(myid, mybegin,myend,i) 
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
        }
    
    fasp_icsr_free(&ST);
    
    fasp_mem_free(work);
    
    return col;
}

/**
 * \fn static void generate_sparsity_P (dCSRmat *P, iCSRmat *S, ivector *vertices, 
 *                                      INT row, INT col)
 *
 * \brief Find coarse level points
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
 * Modified by Chunsheng Feng, Xiaoqiang Yue
 * \date   05/23/2012
 */

static void generate_sparsity_P ( dCSRmat *P, 
		                          iCSRmat *S, 
								  ivector *vertices, 
								  INT row, 
								  INT col )
{
    INT i,j,k,index=0;
    INT *vec=vertices->val;
    
    P->row=row; P->col=col;    
    P->IA=(INT*)fasp_mem_calloc(row+1, sizeof(INT));    
    
	INT nthreads, use_openmp;

	if(!FASP_USE_OPENMP || row <= OPENMP_HOLDS){
		use_openmp = FALSE;
	}
	else{
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
    
#if CHMEM_MODE
    total_alloc_mem += (row+1)*sizeof(INT);
#endif
    
    // step 1: Find the structure IA of P first
    if (use_openmp) {
        INT mybegin,myend,myid;
#pragma omp parallel for private(myid, mybegin,myend,i,j,k) 
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
    P->JA=(INT*)fasp_mem_calloc(P->nnz,sizeof(INT));    
    P->val=(REAL*)fasp_mem_calloc(P->nnz,sizeof(REAL));
    
#if CHMEM_MODE
    total_alloc_mem += (P->nnz)*sizeof(INT);
    total_alloc_mem += (P->nnz)*sizeof(REAL);
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
}

/**
 * \fn static void generate_sparsity_P_standard (dCSRmat *P, iCSRmat *S, ivector *vertices, INT row, INT col)
 *
 * \brief Find sparsity pattern of P for standard interpolation
 *
 * \param P         Pointer to the prolongation matrix
 * \param S         Pointer to the set of all strong couplings matrix
 * \param vertices  Pointer to the type of variables
 * \param row       Number of rows of P
 * \param col       Number of cols of P
 *
 * \author Kai Yang, Xiaozhe Hu
 * \date   05/21/2012
 */
static void generate_sparsity_P_standard (dCSRmat *P,
                                          iCSRmat *S,
                                          ivector *vertices,
                                          INT row,
                                          INT col)
{
    INT i,j,k,l,h,index=0;
    INT *vec=vertices->val;
    INT *times_visited;
    
    P->row=row; P->col=col;
    P->IA=(INT*)fasp_mem_calloc(row+1, sizeof(INT));
    
    
    times_visited=(INT*)fasp_mem_calloc(row,sizeof(INT)); // to record the number of times a coarse point is visited.
    
    for (i=0; i<row; i++) times_visited[i] = -1;  // initilize
    
#if CHMEM_MODE
    total_alloc_mem += (row+1)*sizeof(INT);
#endif
    
    // step 1: Find the structure IA of P first
    for (i=0;i<row;++i)
    {
        
        if (vec[i]==FGPT)  // if node i is a F point
        {
            
            for (j=S->IA[i];j<S->IA[i+1];j++) {  
                
                k=S->JA[j];
                
                if (vec[k]==CGPT) // if neighbor k of i is a C point
                {   
                    if (times_visited[k] == i) {
                        // already visited, do nothing
                    }
                    else {
                        // have not visited, do something
                        times_visited[k] = i;
                        P->IA[i+1]++;
                    }
                    
                }
                else if( (vec[k]==FGPT)&&(k!=i) ) // if neighbor k of i is a F point and k is not i
                {
                    for(l=S->IA[k];l<S->IA[k+1];l++) // look at the neighbors of k
                    {
                        h=S->JA[l];
                        if(vec[h]==CGPT) // only care about the C point neighbors of k
                        {
                            
                            if (times_visited[h] == i){
                                // already visited, do nothing
                            }
                            else {
                                // have not visited, do something
                                times_visited[h] = i;
                                P->IA[i+1]++;
                            }
                        }
                        
                    }  // end for(l=S->IA[k];l<S->IA[k+1];l++)
                    
                }  // end if (vec[k]==CGPT)
                
            } // end for (j=S->IA[i];j<S->IA[i+1];j++)
        } 
        else if (vec[i]==ISPT) { // if node i is a special F point
            P->IA[i+1]=0;
        }
        else { // if node i is on coarse grid
            P->IA[i+1]=1;
        }  // end if (vec[i]==FGPT)
    } // end for (i=0;i<row;++i)
    
    for (i=0;i<P->row;++i) P->IA[i+1]+=P->IA[i];
    
    P->nnz=P->IA[P->row]-P->IA[0];
    
    // step 2: Find the structure JA of P
    P->JA=(INT*)fasp_mem_calloc(P->nnz,sizeof(INT));
    P->val=(REAL*)fasp_mem_calloc(P->nnz,sizeof(REAL));
    
#if CHMEM_MODE
    total_alloc_mem += (P->nnz)*sizeof(INT);
    total_alloc_mem += (P->nnz)*sizeof(REAL);
#endif
    
    for (i=0; i<row; i++) times_visited[i] = -1;  // reinitilize
    
    for (i=0;i<row;++i)
    {
        
        index = 0;
        
        if (vec[i]==FGPT)  // if node i is a F point
        {
            
            for (j=S->IA[i];j<S->IA[i+1];j++) {  
                
                k=S->JA[j];
                
                if (vec[k]==CGPT) // if neighbor k of i is a C point
                {   
                    if (times_visited[k] == i) {
                        // already visited, do nothing
                    }
                    else {
                        // have not visited, do something
                        times_visited[k] = i;
                        P->JA[P->IA[i]+index] = k;
                        index++;
                    }
                    
                }
                else if( (vec[k]==FGPT)&&(k!=i) ) // if neighbor k of i is a F point and k is not i
                {
                    for(l=S->IA[k];l<S->IA[k+1];l++) // look at the neighbors of k
                    {
                        h=S->JA[l];
                        if(vec[h]==CGPT) // only care about the C point neighbors of k
                        {
                            
                            if (times_visited[h] == i){
                                // already visited, do nothing
                            }
                            else {
                                // have not visited, do something
                                times_visited[h] = i;
                                P->JA[P->IA[i]+index] = h;
                                index++;
                            }
                        }
                        
                    }  // end for(l=S->IA[k];l<S->IA[k+1];l++)
                    
                }  // end if (vec[k]==CGPT)
                
            } // end for (j=S->IA[i];j<S->IA[i+1];j++)
        } 
        else if (vec[i]==ISPT)
        {
            // if node i is a special F point, do nothing
        }
        else // if node i is a C point
        {
            P->JA[P->IA[i]]=i;
        }
    }
    
    // clea up
    fasp_mem_free(times_visited);
    
}

/**
 * \fn static INT form_coarse_level_ag (dCSRmat *A, iCSRmat *S, ivector *vertices, INT row,INT interp_type)
 *
 * \brief Find coarse level points by aggressive coarsening
 *
 * \param A             Pointer to the coefficient matrix
 * \param S             Pointer to the set of all strong couplings matrix
 * \param vertices      Pointer to the type of variables
 * \param row           Number of rows of P
 * \param interp_type   type of interpolation
 *
 * \return Number of cols of P
 * 
 * \author Kai Yang, Xiaozhe Hu
 * \date   09/06/2010
 */
static INT form_coarse_level_ag (dCSRmat *A, 
                                 iCSRmat *S,
                                 ivector *vertices,
                                 INT row,
                                 INT interp_type)
{
    
    
    INT col = 0; // initialize col(P): returning output 
	unsigned INT maxlambda, maxnode, num_left=0;	
	REAL measure, new_meas;
    
	INT *ia=A->IA, *vec=vertices->val;
	INT ci_tilde = -1, ci_tilde_mark = -1;
	INT set_empty = 1, C_i_nonempty = 0;
	INT ji,jj,i,j,k,l,m,flag,index,ci,cj,ck,cl,num_c,count,fj;
	
	INT *work = (INT*)fasp_mem_calloc(4*row,sizeof(INT));	
	INT *lists = work, *where = lists+row, *lambda = where+row, *graph_array = lambda+row;
    INT *cp_index,*cp_rindex;
    
	
	LinkList LoL_head = NULL, LoL_tail = NULL, list_ptr = NULL;	
	
	iCSRmat ST,Sh,ShT; //Sh is for the strong coupling matrix between temporary CGPTs, ShT is the transpose of Sh, Snew is for combining the information from S and Sh
    
	ivector CGPT_index, CGPT_rindex;
	fasp_icsr_trans(S, &ST);
	
    vertices->row=A->row;
    
	/**************************************************/
	/* Coarsening Phase ONE: find temporary coarse level points */
	/**************************************************/	
    num_c=form_coarse_level(A, S, vertices, row, interp_type);//return the number of temporary CGPT 
    
    /**************************************************/
	/* Coarsening Phase TWO: find real coarse level points  */
	/**************************************************/	
	//find Sh, the strong coupling between coarse grid points w.r.t. (p,2)
	generate_S_rs_ag (A, S, &Sh, vertices, &CGPT_index, &CGPT_rindex);
    	
	fasp_icsr_trans(&Sh, &ShT);	
    
    CGPT_index.row=num_c;
    CGPT_rindex.row=row;
	cp_index=CGPT_index.val;
    cp_rindex=CGPT_rindex.val;
	
    
	// 1. Initialize lambda
	for (ci=0;ci<num_c;++ci) lambda[ci]=ShT.IA[ci+1]-ShT.IA[ci];
	
	// 2. Set the variables with nonpositvie measure as F-variables
	num_left=0; // number of CGPT in the list LoL
	for (ci=0;ci<num_c;++ci)
	{
        i=cp_index[ci];
        measure=lambda[ci];
        if (vec[i]!=ISPT) {
            if (measure>0) {
                enter_list(&LoL_head, &LoL_tail, lambda[ci], ci, lists, where);
                num_left++;
			}
            else {
                if (measure<0) printf("### WARNING(form_coarse_level_ag): negative lambda!\n");
                
                vec[i]=FGPT; // set i as fine node
                
                //update the lambda value in the CGPT neighbor of i
                for (ck=Sh.IA[ci];ck<Sh.IA[ci+1];++ck) {
                    cj=Sh.JA[ck];
                    j=cp_index[cj];			
                    if (vec[j]!=ISPT) {
						
                        if (cj<ci) {
                            new_meas=lambda[cj];
                            if (new_meas>0) {
                                remove_node(&LoL_head, &LoL_tail, new_meas, cj, lists, where);
                                num_left--;
                            }
                            new_meas= ++(lambda[cj]);
                            enter_list(&LoL_head, &LoL_tail,  new_meas, cj, lists, where);
                            num_left++;
                        }
                        else{
                            new_meas= ++(lambda[cj]);
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
		if (maxlambda==0) { printf("head of the list has measure of 0\n");}
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
					if (vec[k]==CGPT) // k is temporary CGPT
					{
                        
						remove_node(&LoL_head, &LoL_tail, lambda[ck], ck, lists, where);
						new_meas= ++(lambda[ck]);
						enter_list(&LoL_head, &LoL_tail,new_meas, ck, lists, where);
					}
				}
				
			} // if
		} // ci
        
		for (ci=Sh.IA[maxnode];ci<Sh.IA[maxnode+1];++ci) {
			cj=Sh.JA[ci];
			j=cp_index[cj];
            
			if (vec[j]==CGPT) // j is temporary CGPT
			{
				measure=lambda[cj];
                
				remove_node(&LoL_head, &LoL_tail,measure, cj, lists, where);
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
						if (vec[k]==CGPT) // k is temporary CGPT
						{
                            
							remove_node(&LoL_head, &LoL_tail, lambda[ck], ck, lists, where);
							new_meas= ++(lambda[ck]);
							enter_list(&LoL_head, &LoL_tail,new_meas, ck, lists, where);
						}
					} // end for l
				} // end if
			} // end if
		} // end for
	} // while
    
    //organize the variable type
    // make temporary CGPT--1 and fake CGPT--4 become FGPT
    // make real CGPT--3 to be CGPT
    for(i=0;i<row;i++)
    {
	    if(vec[i]==CGPT||vec[i]==4)//change all the temporary or fake CGPT into FGPT
            vec[i]=FGPT;
    }
    
	for(i=0;i<row;i++)
    {
	    if(vec[i]==3)// if i is real CGPT
            vec[i]=CGPT;
    }
    
	/* Coarsening Phase THREE: Find all the FGPTs which 
       have not CGPT neighbors within distance of 2. 
       Change them into CGPT to make standard interpolation
       work  */    
    for (i=0; i<row; i++) {
        if (vec[i]==FGPT) {
            flag=0; //flag for whether there is any CGPT neighbor within distance of 2
            for (j=S->IA[i]; j<S->IA[i+1]; j++) {
                k=S->JA[j];
                if (flag==1) {
                    break;
                }
                else if(vec[k]==CGPT)
                {
                    flag=1;
                    break;
                }
                else if(vec[k]==FGPT)
                {
                    for (l=S->IA[k]; l<S->IA[k+1]; l++) {
                        m=S->JA[l];
                        if (vec[m]==CGPT) {
                            flag=1;
                            break;
                        }
                    }//end for l
                }
            }//end for j
            
            if (flag==0) {
                vec[i]=CGPT;
                col++;
            }
            
            
        }//end if 
    }//end for i
    
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


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
