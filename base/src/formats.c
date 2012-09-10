/*! \file formats.c
 *  \brief Matrix format conversion routines
 *
 */

#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_format_dcoo_dcsr (dCOOmat *A, dCSRmat *B)
 *
 * \brief Transform a REAL matrix from its IJ format to its CSR format.
 *
 * \param A   Pointer to IJ matrix
 * \param B   Pointer to CSR matrix
 *
 * \return    SUCCESS if succeed
 * 
 * \author Xuehai Huang
 * \date   08/10/2009
 */
SHORT fasp_format_dcoo_dcsr (dCOOmat *A, 
                             dCSRmat *B)
{
    const INT m=A->row, n=A->col, nnz=A->nnz;
    
    fasp_dcsr_alloc(m,n,nnz,B);
    
    INT * ia = B->IA, i;
    INT   iind, jind;
    
    INT * ind = (INT *) fasp_mem_calloc(m+1,sizeof(INT));
    
    for (i=0; i<=m; ++i) ind[i]=0; // initialize    
    
    for (i=0; i<nnz; ++i) ind[A->I[i]+1]++; // count nnz in each row
    
    ia[0] = 0; // first index starting from zero
    for (i=1; i<=m; ++i) {
        ia[i] = ia[i-1]+ind[i]; // set row_idx
        ind[i] = ia[i];
    }
    
    // loop over nnz and set col_idx and val
    for (i=0; i<nnz; ++i) {
        iind = A->I[i]; jind = ind[iind]; 
        B->JA [jind] = A->J[i];
        B->val[jind] = A->val[i];
        ind[iind] = ++jind;
    }
    
    return SUCCESS;
}


/**
 * \fn SHORT fasp_format_dcsr_dcoo (dCSRmat *A, dCOOmat *B)
 *
 * \brief Transform a REAL matrix from its CSR format to its IJ format.
 *
 * \param A   Pointer to CSR matrix
 * \param B   Pointer to IJ matrix
 *
 * \return    SUCCESS if succeed
 * 
 * \author Xuehai Huang
 * \date   08/10/2009
 */
SHORT fasp_format_dcsr_dcoo (dCSRmat *A, 
                             dCOOmat *B)
{
    const INT m=A->row, nnz=A->nnz;
    INT i, j;
    SHORT status = SUCCESS;
    
    B->I=(INT *)fasp_mem_calloc(nnz,sizeof(INT));    
    B->J=(INT *)fasp_mem_calloc(nnz,sizeof(INT));    
    B->val=(REAL *)fasp_mem_calloc(nnz,sizeof(REAL));
    
    for (i=0;i<m;++i) {
        for (j=A->IA[i];j<A->IA[i+1];++j) B->I[N2C(j)]=C2N(i);
    }
    
    memcpy(B->J,  A->JA, nnz*sizeof(INT));
    memcpy(B->val,A->val,nnz*sizeof(REAL));
    
    return status;
}

/**
 * \fn SHORT fasp_format_dstr_dcsr(dSTRmat *A, dCSRmat *B_ptr)
 *
 * \brief Transfer a 'dSTRmat' type matrix into a 'dCSRmat' type matrix.
 *
 * \param A      Pointer to a 'dSTRmat' type matrix
 * \param B_ptr  Pointer to the 'dCSRmat' type matrix
 *
 * \return    SUCCESS if succeed
 * 
 * \author Zhiyang Zhou
 * \date   2010/04/29
 */
SHORT fasp_format_dstr_dcsr (dSTRmat *A, 
                             dCSRmat *B_ptr)
{
    // some members of A
    const INT nc    = A->nc;   
    const INT ngrid = A->ngrid;
    const INT nband = A->nband;
    const INT *offsets = A->offsets;
    
    REAL  *diag = A->diag;
    REAL **offdiag = A->offdiag;
    
    // some members of B
    const INT glo_row = nc*ngrid;
    INT glo_nnz;
    INT *ia = NULL;
    INT *ja = NULL;
    REAL *a = NULL;
    
    dCSRmat B;
    
    // local variables
    INT width;
    INT nc2 = nc*nc;
    INT BAND,ROW,COL;
    INT ncb,nci;
    INT row_start,col_start;
    INT block; // how many blocks in the current ROW
    INT i,j;
    INT pos;
    INT start;
    INT val_L_start,val_R_start;
    INT row;
    INT tmp_col;
    REAL tmp_val;
    
    // allocate for 'ia' array
    ia = (INT *)fasp_mem_calloc(glo_row+1,sizeof(INT));
    
    // Generate the 'ia' array
    ia[0] = 0;
    for (ROW = 0; ROW < ngrid; ++ROW) {
        block = 1; // diagonal block
        for (BAND = 0; BAND < nband; ++BAND) {
            width = offsets[BAND];
            COL   = ROW + width;
            if (width < 0) {
                if (COL >= 0) ++block;
            }
            else {
                if (COL < ngrid) ++block;
            }
        } // end for BAND
        
        ncb = nc*block;
        row_start = ROW*nc;
        
        for (i = 0; i < nc; i ++) {
            row = row_start + i; 
            ia[row+1] = ia[row] + ncb;
        }
    } // end for ROW
    
    // allocate for 'ja' and 'a' arrays
    glo_nnz = ia[glo_row];
    ja = (INT *)fasp_mem_calloc(glo_nnz,sizeof(INT));
    a = (REAL *)fasp_mem_calloc(glo_nnz,sizeof(REAL));
    
    // Generate the 'ja' and 'a' arrays at the same time 
    for (ROW = 0; ROW < ngrid; ++ROW) {
        row_start = ROW*nc;    
        val_L_start = ROW*nc2;
        
        // deal with the diagonal band
        for (i = 0; i < nc; i ++) {
            nci   = nc*i;
            row   = row_start + i;
            start = ia[row];
            for (j = 0; j < nc; j ++) {
                pos     = start + j;
                ja[pos] = row_start + j;
                a[pos]  = diag[val_L_start+nci+j];
            }
        }
        block = 1;
        
        // deal with the off-diagonal bands
        for (BAND = 0; BAND < nband; ++BAND) {
            width     = offsets[BAND];
            COL       = ROW + width;
            ncb       = nc*block;
            col_start = COL*nc;
            
            if (width < 0) {
                if (COL >= 0) {
                    val_R_start = COL*nc2;
                    for (i = 0; i < nc; i ++) {
                        nci = nc*i;
                        row = row_start + i;
                        start = ia[row];
                        for (j = 0 ; j < nc; j ++) {
                            pos     = start + ncb + j;
                            ja[pos] = col_start + j;
                            a[pos]  = offdiag[BAND][val_R_start+nci+j];
                        }
                    }
                    ++block;
                }
            }
            else {
                if (COL < ngrid) {
                    for (i = 0; i < nc; i ++) {
                        nci = nc*i;
                        row = row_start + i;
                        start = ia[row];
                        for (j = 0; j < nc; j ++) {
                            pos = start + ncb + j;
                            ja[pos] = col_start + j;
                            a[pos]  = offdiag[BAND][val_L_start+nci+j];
                        }
                    }
                    ++block;
                }
            }
        }
    }
    
    // Reordering in such manner that every diagonal element 
    // is firstly stored in the corresponding row
    if (nc > 1) {
        for (ROW = 0; ROW < ngrid; ++ROW) {
            row_start = ROW*nc;
            for (j = 1; j < nc; j ++) {
                row   = row_start + j;
                start = ia[row];
                pos   = start + j;
                
                // swap in 'ja'
                tmp_col   = ja[start];
                ja[start] = ja[pos];
                ja[pos]   = tmp_col;
                
                // swap in 'a'
                tmp_val  = a[start];
                a[start] = a[pos];
                a[pos]   = tmp_val;            
            }
        }
    } 
    
    /* fill all the members of B */   
    B.row = glo_row;
    B.col = glo_row;
    B.nnz = glo_nnz;
    B.IA = ia;
    B.JA = ja;
    B.val = a;
    
    *B_ptr = B;
    
    return SUCCESS;
}

/**
 * \fn dCSRmat fasp_format_bdcsr_dcsr (block_dCSRmat *Ab)
 *
 * \brief Form the whole dCSRmat A using blocks given in Ab
 *
 * \param Ab   Pointer to the blocks
 *
 * \return     dCSRmat matrix if succeed, NULL if fail
 * 
 * \author Shiquan Zhang
 * \date   08/10/2010
 */
dCSRmat fasp_format_bdcsr_dcsr (block_dCSRmat *Ab)
{
    INT m=0,n=0,nnz=0;
    const INT mb=Ab->brow, nb=Ab->bcol, nbl=mb*nb;
    dCSRmat **blockptr=Ab->blocks, *blockptrij, A;
    
    INT i,j,ij,ir,i1,length,ilength,start,irmrow,irmrow1;
    INT *row, *col;
    
    row = (INT *)fasp_mem_calloc(mb+1,sizeof(INT));
    col = (INT *)fasp_mem_calloc(nb+1,sizeof(INT));
    
    // count the size of A
    row[0]=0; col[0]=0;
    for (i=0;i<mb;++i) { m+=blockptr[i*nb]->row; row[i+1]=m; }
    for (i=0;i<nb;++i) { n+=blockptr[i]->col;    col[i+1]=n; }
    for (i=0;i<nbl;++i) { nnz+=blockptr[i]->nnz; }
    
    // memory space allocation
    A = fasp_dcsr_create(m,n,nnz);
    
    // set dCSRmat for A
    A.IA[0]=0;
    for (i=0;i<mb;++i) {
        
        for (ir=row[i];ir<row[i+1];ir++) {
            
            for (length=j=0;j<nb;++j) {
                ij=i*nb+j; blockptrij=blockptr[ij];
                if (blockptrij->nnz>0) {
                    start=A.IA[ir]+length;
                    irmrow=ir-row[i]; irmrow1=irmrow+1;
                    ilength=blockptrij->IA[irmrow1]-blockptrij->IA[irmrow];
                    if (ilength>0) {
                        memcpy(&(A.val[start]),&(blockptrij->val[blockptrij->IA[irmrow]]),ilength*sizeof(REAL));
                        memcpy(&(A.JA[start]), &(blockptrij->JA[blockptrij->IA[irmrow]]), ilength*sizeof(INT));
                        for (i1=0;i1<ilength;i1++) A.JA[start+i1]+=col[j];
                        length+=ilength;     
                    }
                }
            } // end for j
            
            A.IA[ir+1]=A.IA[ir]+length;        
        } // end for ir
        
    } // end for i
    
    fasp_mem_free(row);
    fasp_mem_free(col);
    
    return(A);    
}

/**
 * \fn dCSRLmat * fasp_format_dcsrl_dcsr (dCSRmat *A)
 *
 * \brief Convert a dCSRmat into a dCSRLmat
 *
 * \param A   Pointer to the dCSRLmat type matrix
 *
 * \author Zhiyang Zhou
 * \date   2011/01/07
 */
dCSRLmat * fasp_format_dcsrl_dcsr (dCSRmat *A)
{
    REAL   *DATA         = A -> val;
    INT    *IA           = A -> IA;
    INT    *JA           = A -> JA;
    INT     num_rows     = A -> row;
    INT     num_cols     = A -> col;
    INT     num_nonzeros = A -> nnz;
    
    dCSRLmat *B        = NULL;
    INT       dif;
    INT      *nzdifnum = NULL;
    INT      *rowstart = NULL;
    INT      *rowindex = (INT *)fasp_mem_calloc(num_rows, sizeof(INT));
    INT      *ja       = (INT *)fasp_mem_calloc(num_nonzeros, sizeof(INT));
    REAL     *data     = (REAL *)fasp_mem_calloc(num_nonzeros, sizeof(REAL));
    
    /* auxiliary arrays */
    INT *nzrow    = (INT *)fasp_mem_calloc(num_rows, sizeof(INT));
    INT *counter  = NULL;
    INT *invnzdif = NULL;
    
    INT i,j,k,cnt,maxnzrow;
    
    //-----------------------------------------
    //  Generate 'nzrow' and 'maxnzrow'
    //-----------------------------------------
    
    maxnzrow = 0;
    for (i = 0; i < num_rows; i ++) {
        nzrow[i] = IA[i+1] - IA[i];
        if (nzrow[i] > maxnzrow) {
            maxnzrow = nzrow[i];
        }
    }
    /* generate 'counter' */
    counter = (INT *)fasp_mem_calloc(maxnzrow + 1, sizeof(INT));
    for (i = 0; i < num_rows; i ++) {
        counter[nzrow[i]] ++;
    }
    
    //--------------------------------------------
    //  Determine 'dif'
    //-------------------------------------------- 
    
    for (dif = 0, i = 0; i < maxnzrow + 1; i ++) {
        if (counter[i] > 0) dif ++;
    }
    
    //--------------------------------------------
    //  Generate the 'nzdifnum' and 'rowstart'
    //-------------------------------------------- 
    
    nzdifnum = (INT *)fasp_mem_calloc(dif, sizeof(INT));
    invnzdif = (INT *)fasp_mem_calloc(maxnzrow + 1, sizeof(INT));
    rowstart = (INT *)fasp_mem_calloc(dif + 1, sizeof(INT));
    rowstart[0] = 0;
    for (cnt = 0, i = 0; i < maxnzrow + 1; i ++) {
        if (counter[i] > 0) {
            nzdifnum[cnt] = i;
            invnzdif[i] = cnt;
            rowstart[cnt+1] = rowstart[cnt] + counter[i];
            cnt ++;
        }
    }
    
    //--------------------------------------------
    //  Generate the 'rowindex'
    //-------------------------------------------- 
    
    for (i = 0; i < num_rows; i ++) {
        j = invnzdif[nzrow[i]];
        rowindex[rowstart[j]] = i;
        rowstart[j] ++;
    }
    /* recover 'rowstart' */
    for (i = dif; i > 0; i --) {
        rowstart[i] = rowstart[i-1];
    }
    rowstart[0] = 0;
    
    //--------------------------------------------
    //  Generate the 'data' and 'ja'
    //-------------------------------------------- 
    
    for (cnt = 0, i = 0; i < num_rows; i ++) {
        k = rowindex[i];
        for (j = IA[k]; j < IA[k+1]; j ++) {
            data[cnt] = DATA[j];
            ja[cnt] = JA[j];
            cnt ++;
        }
    }
    
    //------------------------------------------------------------
    //  Create and fill a dCSRLmat B
    //------------------------------------------------------------ 
    
    B = fasp_dcsrl_create(num_rows, num_cols, num_nonzeros);
    B -> dif      = dif;
    B -> nz_diff  = nzdifnum;
    B -> index    = rowindex;
    B -> start    = rowstart;
    B -> ja       = ja;
    B -> val      = data;   
    
    //----------------------------
    //  Free the auxiliary arrays
    //----------------------------  
    
    free(nzrow);
    free(counter);
    free(invnzdif);
    
    return B;
}

/*!
 * \fn dCSRmat fasp_format_dbsr_dcsr (dBSRmat *B)
 *
 * \brief Transfer a 'dBSRmat' type matrix into a dCSRmat.
 *
 * \param B Pointer to the 'dBSRmat' type matrix
 *
 * \author  Zhiyang Zhou
 * \date    10/23/2010 
 *
 * Modified by Chunsheng Feng, Xiaoqiang Yue 
 * \date    05/24/2012    
 *
 * \note Works for general nb (Xiaozhe)
 */

dCSRmat fasp_format_dbsr_dcsr(dBSRmat *B)
{
    dCSRmat A;
    
    /* members of B */
    INT     ROW = B->ROW;
    INT     COL = B->COL;
    INT     NNZ = B->NNZ;    
    INT     nb  = B->nb;
    INT    *IA  = B->IA;
    INT    *JA  = B->JA;
    INT     storage_manner = B->storage_manner;
    REAL *val = B->val;
    
    INT jump = nb*nb;
    INT rowA = ROW*nb;
    INT colA = COL*nb;
    INT nzA  = NNZ*jump;
    
    INT     *ia = NULL;
    INT     *ja = NULL;
    REAL  *a  = NULL;
    
    INT i,j,k;
    INT mr,mc;
    INT rowstart0,rowstart,colstart0,colstart;
    INT colblock,nzperrow; 
    REAL *vp = NULL;
    REAL *ap = NULL;
    INT    *jap = NULL;

#ifdef _OPENMP 
    INT stride_i,mybegin,myend,myid;
#endif

	INT nthreads = 1, use_openmp = FALSE;

#ifdef _OPENMP  
	if ( ROW > OPENMP_HOLDS ) {
		use_openmp = TRUE;
        nthreads = FASP_GET_NUM_THREADS();
	}
#endif

    //--------------------------------------------------------
    // Create a CSR Matrix 
    //--------------------------------------------------------
    A = fasp_dcsr_create(rowA, colA, nzA);
    ia = A.IA;
    ja = A.JA;
    a = A.val;
    
    //--------------------------------------------------------------------------
    // Compute the number of nonzeros per row, and after this loop,
    // ia[i],i=1:rowA, will be the number of nonzeros of the (i-1)-th row.
    //--------------------------------------------------------------------------
    
    if (use_openmp) {
#ifdef _OPENMP 
        stride_i = ROW/nthreads;
#pragma omp parallel private(myid, mybegin, myend, i, rowstart, colblock, nzperrow, j) 
        {
            myid = omp_get_thread_num();
            mybegin = myid*stride_i;
            if(myid < nthreads-1) myend = mybegin+stride_i;
            else myend = ROW;
            for (i=mybegin; i<myend; ++i)
            {
                rowstart = i*nb + 1;
                colblock = IA[i+1] - IA[i];
                nzperrow = colblock*nb;
                for (j = 0; j < nb; ++j)
                {
                    ia[rowstart+j] = nzperrow;
                }
            }
        }
#endif
    }
    else {
        for (i = 0; i < ROW; ++i)
        {
            rowstart = i*nb + 1;
            colblock = IA[i+1] - IA[i];
            nzperrow = colblock*nb;
            for (j = 0; j < nb; ++j)
            {
                ia[rowstart+j] = nzperrow;
            }
        }
    }
    
    //-----------------------------------------------------
    // Generate the real 'ia' for CSR of A
    //-----------------------------------------------------
    
    ia[0] = 0;
    for (i = 1; i <= rowA; ++i)
    {
        ia[i] += ia[i-1];
    }
    
    //-----------------------------------------------------
    // Generate 'ja' and 'a' for CSR of A
    //-----------------------------------------------------
    
    switch (storage_manner)
    {
        case 0: // each non-zero block elements are stored in row-major order
        {
            if (use_openmp) {
#ifdef _OPENMP 
#pragma omp parallel private(myid, mybegin, myend, i, k, j, rowstart, colstart, vp, mr, ap, jap, mc) ////num_threads(nthreads)
                {
                    myid = omp_get_thread_num();
                    mybegin = myid*stride_i;
                    if(myid < nthreads-1) myend = mybegin+stride_i;
                    else myend = ROW;
                    for (i=mybegin; i<myend; ++i)
                    {
                        for (k = IA[i]; k < IA[i+1]; ++k)
                        {
                            j = JA[k];
                            rowstart = i*nb;
                            colstart = j*nb;
                            vp = &val[k*jump];
                            for (mr = 0; mr < nb; mr ++)
                            {
                                ap  = &a[ia[rowstart]];
                                jap = &ja[ia[rowstart]];
                                for (mc = 0; mc < nb; mc ++)
                                {
                                    *ap = *vp;
                                    *jap = colstart + mc;
                                    vp ++; ap ++; jap ++;
                                }
                                ia[rowstart] += nb;
                                rowstart ++;
                            }
                        }
                    }
                }
#endif
            }
            else {
                for (i = 0; i < ROW; ++i)
                {
                    for (k = IA[i]; k < IA[i+1]; ++k)
                    {
                        j = JA[k];
                        rowstart = i*nb;
                        colstart = j*nb;
                        vp = &val[k*jump];
                        for (mr = 0; mr < nb; mr ++)
                        {
                            ap  = &a[ia[rowstart]];
                            jap = &ja[ia[rowstart]];
                            for (mc = 0; mc < nb; mc ++)
                            {
                                *ap = *vp;
                                *jap = colstart + mc;
                                vp ++; ap ++; jap ++;
                            }
                            ia[rowstart] += nb;
                            rowstart ++;
                        }
                    }
                }
            }
        }
            break;
            
        case 1: // each non-zero block elements are stored in column-major order
        {
            for (i = 0; i < ROW; ++i)
            {
                for (k = IA[i]; k < IA[i+1]; ++k)
                {
                    j = JA[k];
                    rowstart0 = i*nb;
                    colstart0 = j*nb;
                    vp = &val[k*jump];
                    for (mc = 0; mc < nb; mc ++)
                    {
                        rowstart = rowstart0;
                        colstart = colstart0 + mc;
                        for (mr = 0; mr < nb; mr ++)
                        {
                            a[ia[rowstart]] = *vp; 
                            ja[ia[rowstart]] = colstart; 
                            vp ++; ia[rowstart]++; rowstart++;
                        }
                    }
                }
            }
        }
            break;
    }
    
    //-----------------------------------------------------
    // Map back the real 'ia' for CSR of A
    //-----------------------------------------------------
    
    for (i = rowA; i > 0; i --) {
        ia[i] = ia[i-1];
    }
    ia[0] = 0; 
    
    
    return (A);   
}

/*!
 * \fn dBSRmat fasp_format_dcsr_dbsr ( dCSRmat *B, INT nb )
 *
 * \brief Transfer a dCSRmat type matrix into a dBSRmat.
 *
 * \param B   Pointer to the dCSRmat type matrix
 * \param nb  size of each block 
 *
 * \return    Pointer to the 'dBSRmat' type matrix
 *
 * \author Changhe Qiao
 * \date   03/12/2012
 *
 * Modified by Xiaozhe Hu on 03/13/2012
 */

dBSRmat fasp_format_dcsr_dbsr (dCSRmat *B, 
                               const INT nb)
{
    INT *Is, *Js;
    INT i,j,num, k;
    INT nRow, nCol;
    
    dCSRmat tmpMat;
    dBSRmat A;
    
    // Safe-guard check --Chensong 05/27/2012
    if ((B->row)%nb!=0) {
        printf("### ERROR: B.row=%d is not a multiplication of nb=%d!!!\n", B->row, nb);
        exit(0);
    }
    
    if ((B->col)%nb!=0) {
        printf("### ERROR: B.col=%d is not a multiplication of nb=%d!!!\n", B->col, nb);
        exit(0);
    }
    
    nRow=B->row/nb; //we must ensure this is a integer
    nCol=B->col/nb;
    
    Is=(INT *)fasp_mem_calloc(nRow, sizeof(INT));
    Js=(INT *)fasp_mem_calloc(nCol, sizeof(INT));
    for(i=0;i<nRow;i++) {
        Is[i]=i*nb;
    }
    for(i=0;i<nCol;i++) {
        Js[i]=i*nb;
    }
    
    fasp_dcsr_getblk(B,Is,Js,nRow,nCol,&tmpMat);
    
    //here we have tmpmat as the submatrix
    A.ROW=nRow;
    A.COL=nCol;
    A.NNZ=tmpMat.nnz;
    A.nb=nb;
    A.storage_manner=0;//row majored
    A.IA=(INT*)fasp_mem_calloc(A.ROW+1, sizeof(INT));
    A.JA=(INT*)fasp_mem_calloc(A.NNZ, sizeof(INT));
    A.val=(REAL*)fasp_mem_calloc(A.NNZ*nb*nb, sizeof(REAL));
    for(i=0;i<tmpMat.row+1;i++) {
        A.IA[i]=tmpMat.IA[i];
    }
    for(i=0;i<tmpMat.nnz;i++) {
        A.JA[i]=tmpMat.JA[i];
    }
    for(i=0;i<tmpMat.nnz;i++) {
        A.val[i*nb*nb]=tmpMat.val[i];
    }
    fasp_mem_free(tmpMat.IA);
    fasp_mem_free(tmpMat.JA);
    fasp_mem_free(tmpMat.val);
    
    for(i=0;i<nb;i++) {
        for(j=0;j<nb;j++) {
            num=i*nb+j;
            if (i==0 && j==0) {
            }
            else {
                for(k=0;k<nRow;k++) {
                    Is[k]=k*nb+i;
                }
                for(k=0;k<nCol;k++) {
                    Js[k]=k*nb+j;
                }
                fasp_dcsr_getblk(B,Is,Js,nRow,nCol,&tmpMat);
                for(k=0;k<tmpMat.nnz;k++) {
                    A.val[k*nb*nb+num]=tmpMat.val[k];
                }    
                fasp_mem_free(tmpMat.IA);
                fasp_mem_free(tmpMat.JA);
                fasp_mem_free(tmpMat.val);
            }
        }
    }
    
    return A;
}

/*!
 * \fn dBSRmat fasp_format_dstr_dbsr ( dSTRmat *B )
 *
 * \brief Transfer a 'dSTRmat' type matrix to a 'dBSRmat' type matrix.
 *
 * \param B   Pointer to the 'dSTRmat' type matrix
 *
 * \author Zhiyang Zhou
 * \date   2010/10/26 
 */
dBSRmat fasp_format_dstr_dbsr (dSTRmat *B)
{
    // members of 'B'
    INT      nc      = B->nc;
    INT      ngrid   = B->ngrid;
    REAL    *diag    = B->diag;
    INT      nband   = B->nband;
    INT     *offsets = B->offsets;
    REAL   **offdiag = B->offdiag;
    
    // members of 'A'
    dBSRmat  A;
    INT      NNZ;
    INT     *IA  = NULL;
    INT     *JA  = NULL;
    REAL    *val = NULL;
    
    // local variables
    INT i,j,k,m;
    INT nc2 = nc*nc;
    INT ngridplus1 = ngrid + 1;
    
    // compute NNZ
    NNZ = ngrid;
    for (i = 0; i < nband; ++i) {
        NNZ += (ngrid - abs(offsets[i]));
    } 
    
    // Create and Initialize a dBSRmat 'A'
    A = fasp_dbsr_create(ngrid, ngrid, NNZ, nc, 0);
    IA = A.IA;
    JA = A.JA;
    val = A.val;
    
    // Generate 'IA'
    for (i = 1; i < ngridplus1; ++i) IA[i] = 1; // take the diagonal blocks into account
    for (i = 0; i < nband; ++i) {
        k = offsets[i];
        if (k < 0) {
            for (j = -k+1; j < ngridplus1; ++j) {
                IA[j] ++;
            }
        }
        else {
            m = ngridplus1 - k;
            for (j = 1; j < m; ++j)
            {
                IA[j] ++;
            }
        }
    }
    IA[0] = 0;
    for (i = 1; i < ngridplus1; ++i) {
        IA[i] += IA[i-1];
    }
    
    // Generate 'JA' and 'val' at the same time
    for (i = 0 ; i < ngrid; ++i) {
        memcpy(val + IA[i]*nc2, diag + i*nc2, nc2*sizeof(REAL));
        JA[IA[i]] = i;
        IA[i] ++;
    }
    
    for (i = 0; i < nband; ++i) {
        k = offsets[i];
        if (k < 0) {
            for (j = -k; j < ngrid; ++j) {
                m = j + k;
                memcpy(val+IA[j]*nc2, offdiag[i]+m*nc2, nc2*sizeof(REAL));
                JA[IA[j]] = m;
                IA[j] ++;
            }
        }
        else {
            m = ngrid - k;
            for (j = 0; j < m; ++j) {
                memcpy(val + IA[j]*nc2, offdiag[i] + j*nc2, nc2*sizeof(REAL));
                JA[IA[j]] = k + j;
                IA[j] ++;
            }
        }
    }
    
    // Map back the real 'IA' for BSR of A
    for (i = ngrid; i > 0; i --) {
        IA[i] = IA[i-1];
    }
    IA[0] = 0; 
    
    return (A);
}

/*!
 * \fn dCOOmat * fasp_format_dbsr_dcoo ( dBSRmat *B )
 *
 * \brief Transfer a 'dBSRmat' type matrix to a 'dCOOmat' type matrix.
 *
 * \param B   Pointer to the 'dBSRmat' type matrix
 *
 * \author Zhiyang Zhou
 * \date   2010/10/26 
 */
dCOOmat * fasp_format_dbsr_dcoo (dBSRmat *B)
{
    /* members of B */
    INT     ROW = B->ROW;
    INT     COL = B->COL;
    INT     NNZ = B->NNZ;    
    INT     nb  = B->nb;
    INT    *IA  = B->IA;
    INT    *JA  = B->JA;
    REAL   *val = B->val;
    
    dCOOmat *A = NULL;
    INT      nb2 = nb*nb;
    INT      num_nonzeros = NNZ*nb2;
    INT     *rowA = NULL;
    INT     *colA = NULL;
    REAL    *valA = NULL;
    
    INT      i,j,k,inb;
    INT      row_start, col_start;
    INT      cnt,mr,mc;
    REAL    *pt = NULL;
    
    // Create and Initialize a dBSRmat 'A'
    A = (dCOOmat *)fasp_mem_calloc(1, sizeof(dCOOmat));
    A->row = ROW*nb;
    A->col = COL*nb;
    A->nnz = num_nonzeros;    
    rowA   = (INT *)fasp_mem_calloc(num_nonzeros, sizeof(INT)); 
    colA   = (INT *)fasp_mem_calloc(num_nonzeros, sizeof(INT));
    valA   = (REAL *)fasp_mem_calloc(num_nonzeros, sizeof(REAL));
    A->I   = rowA;
    A->J   = colA;
    A->val = valA;
    
    cnt = 0;
    for (i = 0; i < ROW; ++i) {
        inb = i*nb;
        for (k = IA[i]; k < IA[i+1]; ++k) {
            j  = JA[k];
            pt = &val[k*nb2];
            row_start = inb; 
            col_start = j*nb;
            for (mr = 0; mr < nb; mr ++) {
                for (mc = 0; mc < nb; mc ++) {
                    rowA[cnt] = row_start;
                    colA[cnt] = col_start + mc;
                    valA[cnt] = (*pt);
                    pt ++;
                    cnt ++;
                }
                row_start ++;
            }
        }
    }
    
    return (A);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
