/*! \file AuxSort.c
 *
 *  \brief Subroutines for ordering, merging, removing duplicated integers
 *
 *  \note This file contains Level-0 (Aux) functions. It requires
 *        AuxMemory.c
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/

static void dSwapping(REAL *w, const INT i, const INT j);
static void iSwapping(INT *w, const INT i, const INT j);
static void CMK_ordering (const dCSRmat *, INT, INT, INT, INT, INT *, INT *);
static void multicoloring (AMG_data *, REAL, INT *, INT *);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn INT fasp_BinarySearch (INT *list, const INT value, const INT nlist)
 *
 * \brief Binary Search
 *
 * \param list     Pointer to a set of values
 * \param value    The target
 * \param nlist    Length of the array list
 *
 * \return  The location of value in array list if succeeded; otherwise, return -1.
 *
 * \author Chunsheng Feng
 * \date   03/01/2011
 */
INT fasp_BinarySearch (INT       *list,
                       const INT  value,
                       const INT  nlist)
{
    INT not_found = 1;
    INT low, high, m;
    
    low = 0;
    high = nlist - 1;
    while (not_found && low <= high)
    {
        m = (low + high) / 2;
        if (value < list[m])
        {
            high = m - 1;
        }
        else if (value > list[m])
        {
            low = m + 1;
        }
        else
        {
            not_found = 0;
            return m;
        }
    }
    
    return -1;
}

/*!
 * \fn INT fasp_aux_unique (INT numbers[], const INT size)
 *
 * \brief Remove duplicates in an sorted (ascending order) array
 *
 * \param numbers   Pointer to the array needed to be sorted (in/out)
 * \param size      Length of the target array
 *
 * \return          New size after removing duplicates
 *
 * \author Chensong Zhang
 * \date   11/21/2010
 *
 * \note Operation is in place. Does not use any extra or temporary storage.
 */
INT fasp_aux_unique (INT        numbers[],
                     const INT  size)
{
    INT i, newsize;
    
    if (size==0) return(0);
    
    for (newsize=0, i=1; i<size; ++i) {
        if (numbers[newsize] < numbers[i]) {
            newsize++;
            numbers[newsize] = numbers[i];
        }
    }
    
    return(newsize+1);
}

/*!
 * \fn void fasp_aux_merge (INT numbers[], INT work[], INT left, INT mid, INT right)
 *
 * \brief Merge two sorted arrays
 *
 * \param numbers   Pointer to the array needed to be sorted
 * \param work      Pointer to the work array with same size as numbers
 * \param left      Starting index of array 1
 * \param mid       Starting index of array 2
 * \param right     Ending index of array 1 and 2
 *
 * \author Chensong Zhang
 * \date   11/21/2010
 *
 * \note Both arrays are stored in numbers! Arrays should be pre-sorted!
 */
void fasp_aux_merge (INT numbers[],
                     INT work[],
                     INT left,
                     INT mid,
                     INT right)
{
    INT i, left_end, num_elements, tmp_pos;
    
    left_end = mid - 1;
    tmp_pos = left;
    num_elements = right - left + 1;
    
    while ((left <= left_end) && (mid <= right)) {
        
        if (numbers[left] <= numbers[mid]) // first branch <=
        {
            work[tmp_pos] = numbers[left];
            tmp_pos = tmp_pos + 1;
            left = left +1;
        }
        else // second branch >
        {
            work[tmp_pos] = numbers[mid];
            tmp_pos = tmp_pos + 1;
            mid = mid + 1;
        }
    }
    
    while (left <= left_end) {
        work[tmp_pos] = numbers[left];
        left = left + 1;
        tmp_pos = tmp_pos + 1;
    }
    
    while (mid <= right) {
        work[tmp_pos] = numbers[mid];
        mid = mid + 1;
        tmp_pos = tmp_pos + 1;
    }
    
    for (i = 0; i < num_elements; ++i) {
        numbers[right] = work[right];
        right = right - 1;
    }
    
}

/*!
 * \fn void fasp_aux_msort (INT numbers[], INT work[], INT left, INT right)
 *
 * \brief Sort the INT array in ascending order with the merge sort algorithm
 *
 * \param numbers   Pointer to the array needed to be sorted
 * \param work      Pointer to the work array with same size as numbers
 * \param left      Starting index
 * \param right     Ending index
 *
 * \author Chensong Zhang
 * \date   11/21/2010
 *
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively
 */
void fasp_aux_msort (INT numbers[],
                     INT work[],
                     INT left,
                     INT right)
{
    INT mid;
    
    if (right > left) {
        mid = (right + left) / 2;
        fasp_aux_msort(numbers, work, left, mid);
        fasp_aux_msort(numbers, work, mid+1, right);
        fasp_aux_merge(numbers, work, left, mid+1, right);
    }
    
}

/*!
 * \fn void fasp_aux_iQuickSort (INT *a, INT left, INT right)
 *
 * \brief Sort the array (INT type) in ascending order  with the quick sorting algorithm
 *
 * \param a      Pointer to the array needed to be sorted
 * \param left   Starting index
 * \param right  Ending index
 *
 * \author Zhiyang Zhou
 * \date   11/28/2009
 *
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively
 *        where n is the length of 'a'.
 */
void fasp_aux_iQuickSort (INT *a,
                          INT  left,
                          INT  right)
{
    INT i, last;
    
    if (left >= right) return;
    
    iSwapping(a, left, (left+right)/2);
    
    last = left;
    for (i = left+1; i <= right; ++i) {
        if (a[i] < a[left]) {
            iSwapping(a, ++last, i);
        }
    }
    
    iSwapping(a, left, last);
    
    fasp_aux_iQuickSort(a, left, last-1);
    fasp_aux_iQuickSort(a, last+1, right);
}

/*!
 * \fn void fasp_aux_dQuickSort (REAL *a, INT left, INT right)
 *
 * \brief Sort the array (REAL type) in ascending order  with the quick sorting algorithm
 *
 * \param a      Pointer to the array needed to be sorted
 * \param left   Starting index
 * \param right  Ending index
 *
 * \author Zhiyang Zhou
 * \date   2009/11/28
 *
 * \note 'left' and 'right' are usually set to be 0 and n-1, respectively
 *        where n is the length of 'a'.
 */
void fasp_aux_dQuickSort (REAL *a,
                          INT  left,
                          INT  right)
{
    INT i, last;
    
    if (left >= right) return;
    
    dSwapping(a, left, (left+right)/2);
    
    last = left;
    for (i = left+1; i <= right; ++i) {
        if (a[i] < a[left]) {
            dSwapping(a, ++last, i);
        }
    }
    
    dSwapping(a, left, last);
    
    fasp_aux_dQuickSort(a, left, last-1);
    fasp_aux_dQuickSort(a, last+1, right);
}

/*!
 * \fn void fasp_aux_iQuickSortIndex (INT *a, INT left, INT right, INT *index)
 *
 * \brief Reorder the index of (INT type) so that 'a' is in ascending order
 *
 * \param a       Pointer to the array
 * \param left    Starting index
 * \param right   Ending index
 * \param index   Index of 'a' (out)
 *
 * \author Zhiyang Zhou
 * \date   2009/12/02
 *
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the
 *       length of 'a'. 'index' should be initialized in the nature order and it has the
 *       same length as 'a'.
 */
void fasp_aux_iQuickSortIndex (INT *a,
                               INT  left,
                               INT  right,
                               INT *index)
{
    INT i, last;
    
    if (left >= right) return;
    
    iSwapping(index, left, (left+right)/2);
    
    last = left;
    for (i = left+1; i <= right; ++i) {
        if (a[index[i]] < a[index[left]]) {
            iSwapping(index, ++last, i);
        }
    }
    
    iSwapping(index, left, last);
    
    fasp_aux_iQuickSortIndex(a, left, last-1, index);
    fasp_aux_iQuickSortIndex(a, last+1, right, index);
}

/*!
 * \fn void fasp_aux_dQuickSortIndex(REAL *a, INT left, INT right, INT *index)
 *
 * \brief Reorder the index of (REAL type) so that 'a' is ascending in such order
 *
 * \param a       Pointer to the array
 * \param left    Starting index
 * \param right   Ending index
 * \param index   Index of 'a' (out)
 *
 * \author Zhiyang Zhou
 * \date   2009/12/02
 *
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the
 *       length of 'a'. 'index' should be initialized in the nature order and it has the
 *       same length as 'a'.
 */
void fasp_aux_dQuickSortIndex (REAL *a,
                               INT   left,
                               INT   right,
                               INT  *index)
{
    INT i, last;
    
    if (left >= right) return;
    
    iSwapping(index, left, (left+right)/2);
    
    last = left;
    for (i = left+1; i <= right; ++i) {
        if (a[index[i]] < a[index[left]]) {
            iSwapping(index, ++last, i);
        }
    }
    
    iSwapping(index, left, last);
    
    fasp_aux_dQuickSortIndex(a, left, last-1, index);
    fasp_aux_dQuickSortIndex(a, last+1, right, index);
}

/**
 * \fn void fasp_dcsr_CMK_order (const dCSRmat *A, INT *order, INT *oindex)
 *
 * \brief Ordering vertices of matrix graph corresponding to A.
 *
 * \param A       Pointer to matrix
 * \param oindex  Pointer to index of vertices in order
 * \param order   Pointer to vertices with increasing degree
 *
 * \author Zheng Li, Chensong Zhang
 * \date   05/28/2014
 */
void fasp_dcsr_CMK_order (const dCSRmat *A,
                          INT           *order,
                          INT           *oindex)
{
    const INT *ia = A->IA;
    const INT row= A->row;
    
    INT i, loc, s, vt, mindg, innz;
    
    s = 0;
    vt = 0;
    mindg = row+1;
    
    // select node with minimal degree
    for (i=0; i<row; ++i) {
        innz = ia[i+1] - ia[i];
        if (innz > 1) {
            oindex[i] = -innz;
            if (innz < mindg) {
                mindg = innz;
                vt = i;
            }
        }
        else { // order those diagonal rows first
            oindex[i] = s;
            order[s] = i;
            s ++;
        }
    }
    
    loc = s;
    
    // start to order
    CMK_ordering (A, loc, s, vt, mindg, oindex, order);
}

/**
 * \fn void fasp_dcsr_RCMK_order (const dCSRmat *A, INT *order, INT *oindex, INT *rorder)
 *
 * \brief  Resverse CMK ordering
 *
 * \param A       Pointer to matrix
 * \param order   Pointer to vertices with increasing degree
 * \param oindex  Pointer to index of vertices in order
 * \param rorder  Pointer to reverse order
 *
 * \author Zheng Li, Chensong Zhang
 * \date   10/10/2014
 */
void fasp_dcsr_RCMK_order (const dCSRmat *A,
                           INT           *order,
                           INT           *oindex,
                           INT           *rorder)
{
    INT i;
    INT row = A->row;
    
    // Form CMK order
    fasp_dcsr_CMK_order(A, order, oindex);
    
    // Reverse CMK order
    for (i=0; i<row; ++i) rorder[i] = order[row-1-i];
}

/**
 * \fn void fasp_topological_sorting_ilu(ILU_data *iludata)
 *
 * \brief Reordering vertices according to level schedule strategy
 *
 * \param iludata  Pointer to iludata
 *
 * \author Zheng Li, Chensong Zhang
 * \date   12/04/2016
 */
void fasp_topological_sorting_ilu (ILU_data *iludata)
{
    int i, j, k, l;
    int nlevL, nlevU;
    
    int n = iludata->row;
    int *ijlu = iludata->ijlu;
    
    int *level = (int*)fasp_mem_calloc(n, sizeof(int));
    int *jlevL = (int*)fasp_mem_calloc(n, sizeof(int));
    int *ilevL = (int*)fasp_mem_calloc(n+1, sizeof(int));
    
    nlevL = 0;
    ilevL[0] = 0;
    
    // form level for each row of lower trianguler matrix.
    for (i=0; i<n; i++) {
        l = 0;
        for(j=ijlu[i]; j<ijlu[i+1]; j++) if (ijlu[j]<=i) l = MAX(l, level[ijlu[j]]);
        level[i] = l+1;
        ilevL[l+1] ++;
        nlevL = MAX(nlevL, l+1);
    }
    
    for (i=1; i<=nlevL; i++) ilevL[i] += ilevL[i-1];
    
    for (i=0; i<n; i++) {
        k = ilevL[level[i]-1];
        jlevL[k] = i;
        ilevL[level[i]-1]++;
    }
    
    for (i=nlevL-1; i>0; i--) ilevL[i] = ilevL[i-1];
    
    // form level for each row of upper trianguler matrix.
    nlevU = 0;
    ilevL[0] = 0;
    
    int *jlevU = (int*)fasp_mem_calloc(n, sizeof(int));
    int *ilevU = (int*)fasp_mem_calloc(n+1, sizeof(int));
    
    for (i=0; i<n; i++) level[i] = 0;
    
    ilevU[0] = 0;
    
    for (i=n-1; i>=0; i--) {
        l = 0;
        for (j=ijlu[i]; j<ijlu[i+1]; j++) if (ijlu[j]>=i) l = MAX(l, level[ijlu[j]]);
        level[i] = l+1;
        ilevU[l+1] ++;
        nlevU = MAX(nlevU, l+1);
    }
    
    for (i=1; i<=nlevU; i++) ilevU[i] += ilevU[i-1];
    
    for (i=n-1; i>=0; i--) {
        k = ilevU[level[i]-1];
        jlevU[k] = i;
        ilevU[level[i]-1]++;
    }
    
    for (i=nlevU-1; i>0; i--) ilevU[i] = ilevU[i-1];
    
    ilevU[0] = 0;
    
    iludata->nlevL = nlevL+1; iludata->ilevL = ilevL;iludata->jlevL = jlevL;
    iludata->nlevU = nlevU+1; iludata->ilevU = ilevU;iludata->jlevU = jlevU;
    
    fasp_mem_free(level);
}

/**
 * \fn void fasp_multicolors_independent_set (AMG_data *mgl, INT gslvl)
 *
 * \brief Coloring vertices of adjacency graph of A
 *
 * \param mgl      Pointer to input matrix
 * \param gslvl    Used to specify levels of AMG using multicolor smoothing
 *
 * \author Zheng Li, Chunsheng Feng
 * \date   12/04/2016
 */
void fasp_multicolors_independent_set (AMG_data *mgl,
                                       INT       gslvl)
{
    
    INT Colors, rowmax, level, prtlvl = 0;
    
    REAL theta = 0.00;
    
    INT maxlvl = MIN(gslvl, mgl->num_levels-1);
    
#ifdef _OPENMP
#pragma omp parallel for private(level,rowmax,Colors) schedule(static, 1)
#endif
    for ( level=0; level<maxlvl; level++ ) {
        
        multicoloring(&mgl[level], theta, &rowmax, &Colors);
        
        // print
        if ( prtlvl > PRINT_MIN )
            printf("mgl[%3d].A.row = %12d rowmax = %5d rowavg = %7.2lf colors = %5d theta = %le\n",
                   level, mgl[level].A.row, rowmax, (double)mgl[level].A.nnz/mgl[level].A.row,
                   mgl[level].colors, theta);
    }
}

/*---------------------------------*/
/*--      Private Functions       --*/
/*---------------------------------*/

/**
 * \fn static void iSwapping (INT *w, const INT i, const INT j)
 *
 * \brief swap the i-th and j-th element in the array 'w' (INT type)
 *
 * \param w    Pointer to the array
 * \param i    One entry in w
 * \param j    Another entry in w
 *
 * \author Zhiyang Zhou
 * \date   2009/11/28
 */
static void iSwapping (INT       *w,
                       const INT  i,
                       const INT  j)
{
    INT temp = w[i];
    w[i] = w[j];
    w[j] = temp;
}

/**
 * \fn static void dSwapping (REAL *w, const INT i, const INT j)
 *
 * \brief swap the i-th and j-th element in the array 'w' (REAL type)
 *
 * \param w    Pointer to the array
 * \param i    One entry in w
 * \param j    Another entry in w
 *
 * \author Zhiyang Zhou
 * \date   2009/11/28
 */
static void dSwapping (REAL      *w,
                       const INT  i,
                       const INT  j)
{
    REAL temp = w[i];
    w[i] = w[j];
    w[j] = temp;
}

/**
 * \fn static void CMK_ordering (const dCSRmat *A, INT loc, INT s, INT jj,
 *                               INT mindg, INT *oindex, INT *order)
 *
 * \brief CMK ordering by increasing degree of vertices.
 *
 * \param A       Pointer to matrix
 * \param loc     Main order loop variable
 * \param s       Number of ordered vertices
 * \param jj      Vertices with minimal degree
 * \param mindg   Minimal degree
 * \param oindex  Pointer to index of vertices in order
 * \param order   Pointer to vertices with increasing degree
 *
 * \author Zheng Li, Chensong Zhang
 * \date   05/28/2014
 */
static void CMK_ordering (const dCSRmat *A,
                          INT            loc,
                          INT            s,
                          INT            jj,
                          INT            mindg,
                          INT           *oindex,
                          INT           *order)
{
    const INT *ia = A->IA;
    const INT *ja = A->JA;
    const INT row= A->row;
    
    INT i, j, sp1, k, flag;
    
    if (s < row) {
        order[s] = jj;
        oindex[jj] = s;
    }
    
    while (loc <= s && s < row) {
        i = order[loc];
        sp1 = s+1;
        // neighbor nodes are priority.
        for (j=ia[i]+1; j<ia[i+1]; ++j) {
            k = ja[j];
            if (oindex[k] < 0){
                s++;
                order[s] = k;
            }
        }
        // ordering neighbor nodes by increasing degree
        if (s > sp1) {
            while (flag) {
                flag = 0;
                for (i=sp1+1; i<=s; ++i) {
                    if (oindex[order[i]] > oindex[order[i-1]]) {
                        j = order[i];
                        order[i] = order[i-1];
                        order[i-1] = j;
                        flag = 1;
                    }
                }
            }
        }
        
        for (i=sp1; i<=s; ++i) oindex[order[i]] = i;
        
        loc ++;
    }
    
    // deal with remainder
    if (s < row) {
        jj = 0;
        i  = 0;
        while (jj == 0) {
            i ++;
            if (i >= row) {
                mindg++;
                i = 0;
            }
            if (oindex[i] < 0 && (ia[i+1]-ia[i] == mindg)) {
                jj = i;
            }
        }
        
        s ++;
        
        CMK_ordering (A, loc, s, jj, mindg, oindex, order);
    }
}

/**
 * \fn void generate_S_theta (dCSRmat *A, iCSRmat *S, REAL theta)
 *
 * \brief Generate strong sparsity pattern of A
 *
 * \param A      Pointer to input matrix
 * \param S      Pointer to strong sparsity pattern matrix
 * \param theta  Threshold
 *
 * \author Zheng Li, Chunsheng Feng
 * \date   12/04/2016
 */
static void generate_S_theta (dCSRmat *A,
                              iCSRmat *S,
                              REAL     theta)
{
    const INT row=A->row, col=A->col;
    const INT row_plus_one = row+1;
    const INT nnz=A->IA[row]-A->IA[0];
    
    INT index, i, j, begin_row, end_row;
    INT *ia=A->IA, *ja=A->JA;
    REAL *aj=A->val;
    
    // get the diagnal entry of A
    //dvector diag; fasp_dcsr_getdiag(0, A, &diag);
    
    /* generate S */
    REAL row_abs_sum;
    
    // copy the structure of A to S
    S->row=row; S->col=col; S->nnz=nnz; S->val=NULL;
    
    S->IA=(INT*)fasp_mem_calloc(row_plus_one, sizeof(INT));
    
    S->JA=(INT*)fasp_mem_calloc(nnz, sizeof(INT));
    
    fasp_iarray_cp(row_plus_one, ia, S->IA);
    fasp_iarray_cp(nnz, ja, S->JA);
    
    for (i=0;i<row;++i) {
        /* compute scaling factor and row sum */
        row_abs_sum=0;
        
        begin_row=ia[i]; end_row=ia[i+1];
        
        for (j=begin_row;j<end_row;j++) row_abs_sum+=ABS(aj[j]);
        
        row_abs_sum = row_abs_sum*theta;
        
        /* deal with the diagonal element of S */
        //  for (j=begin_row;j<end_row;j++) {
        //     if (ja[j]==i) {S->JA[j]=-1; break;}
        //  }
        
        /* deal with  the element of S */
        for (j=begin_row;j<end_row;j++){
            /* if $\sum_{j=1}^n |a_{ij}|*theta>= |a_{ij}|$ */
            if ( (row_abs_sum >= ABS(aj[j])) && (ja[j] !=i) ) S->JA[j]=-1;
        }
    } // end for i
    
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
}

/**
 * \fn void multicoloring (AMG_data *mgl, REAL theta, INT *rowmax, INT *groups)
 *
 * \brief Coloring vertices of adjacency graph of A
 *
 * \param mgl      Pointer to input matrix
 * \param theta    Threshold
 * \param rowmax   Pointer to number of each color
 * \param groups   Pointer to index array
 *
 * \author Zheng Li, Chunsheng Feng
 * \date   12/04/2016
 */
static void multicoloring (AMG_data *mgl,
                           REAL      theta,
                           INT      *rowmax,
                           INT      *groups)
{
    INT k, i, j, pre, group, iend;
    INT icount;
    INT front, rear;
    INT *IA,*JA;

    const INT n = mgl->A.row;
    dCSRmat   A = mgl->A;
    iCSRmat   S;
    
    if (theta > 0 && theta < 1.0) {
        generate_S_theta(&A, &S, theta);
        IA = S.IA;
        JA = S.JA;
    } else if (theta == 1.0 ) {
        
        mgl->ic = (INT*)malloc(sizeof(INT)*2);
        mgl->icmap = (INT *)malloc(sizeof(INT)*(n+1));
        mgl->ic[0] = 0;
        mgl->ic[1] = n;
        for(k=0; k<n; k++)  mgl->icmap[k]= k;
        
        mgl->colors = 1;
        *groups = 1;
        *rowmax = 1;
        
        printf("### WARNING: Theta = %lf \n", theta);
        
        return;
        
    } else{
        IA = A.IA;
        JA = A.JA;
    }
    
    INT *cq = (INT *)malloc(sizeof(INT)*(n+1));
    INT *newr = (INT *)malloc(sizeof(INT)*(n+1));
    
#ifdef _OPENMP
#pragma omp parallel for private(k)
#endif
    for(k=0;k<n;k++) {
        cq[k]= k;
    }
    group = 0;
    for(k=0;k<n;k++) {
        if ((A.IA[k+1] - A.IA[k]) > group ) group = A.IA[k+1] - A.IA[k];
    }
    *rowmax = group;
    
    mgl->ic = (INT *)malloc(sizeof(INT)*(group+2));
    mgl->icmap = (INT *)malloc(sizeof(INT)*(n+1));
    
    front = n-1;
    rear = n-1;
    
    memset(newr, -1, sizeof(INT)*(n+1));
    memset(mgl->icmap, 0, sizeof(INT)*n);
    
    group=0;
    icount = 0;
    mgl->ic[0] = 0;
    pre=0;
    
    do {
        //front = (front+1)%n;
        front ++;
        if (front == n ) front =0; // front = front < n ? front : 0 ;
        i = cq[front];
        
        if(i <= pre) {
            mgl->ic[group] = icount;
            mgl->icmap[icount] = i;
            group++;
            icount++;
#if 0
            if ((IA[i+1]-IA[i]) > igold)
                iend = MIN(IA[i+1], (IA[i] + igold));
            else
#endif
                iend = IA[i+1];
            
            for(j= IA[i]; j< iend; j++)  newr[JA[j]] = group;
        }
        else if (newr[i] == group) {
            //rear = (rear +1)%n;
            rear ++;
            if (rear == n) rear = 0;
            cq[rear] = i;
        }
        else {
            mgl->icmap[icount] = i;
            icount++;
#if  0
            if ((IA[i+1] - IA[i]) > igold)  iend =MIN(IA[i+1], (IA[i] + igold));
            else
#endif
                iend = IA[i+1];
            for(j = IA[i]; j< iend; j++)  newr[JA[j]] =  group;
        }
        pre=i;
        
    } while(rear != front);
    
    mgl->ic[group] = icount;
    mgl->colors = group;
    *groups = group;
    
    free(cq);
    free(newr);
    
    if (theta >0 ){
        fasp_mem_free(S.IA);
        fasp_mem_free(S.JA);
    }
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
