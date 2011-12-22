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
 * \fn int fasp_format_dcoo_dcsr(dCOOmat *A, dCSRmat *B)
 * \brief Transform a double matrix from its IJ format to its CSR format.
 *
 * \param *A   pointer to IJ matrix
 * \param *B   pointer to CSR matrix
 * \return     SUCCESS if succeed
 * 
 * \author Xuehai Huang
 * \date 08/10/2009
 */
int fasp_format_dcoo_dcsr (dCOOmat *A, 
                           dCSRmat *B)
{
	const int m=A->row, n=A->col, nnz=A->nnz;
	fasp_dcsr_alloc(m,n,nnz,B);
    
	int* ia  = B->IA, i;
    int  iind, jind;
    int* ind = (int *) fasp_mem_calloc(m+1,sizeof(int));
	
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
 * \fn int fasp_format_dcsr_dcoo(dCSRmat *A, dCOOmat *B)
 * \brief Transform a double matrix from its CSR format to its IJ format.
 *
 * \param *A   pointer to CSR matrix
 * \param *B   pointer to IJ matrix
 * \return     SUCCESS if succeed 
 * 
 * \author Xuehai Huang
 * \date 08/10/2009
 */
int fasp_format_dcsr_dcoo (dCSRmat *A, 
                           dCOOmat *B)
{
	const int m=A->row, nnz=A->nnz;
	int i, j;
	int status = SUCCESS;
	
	B->I=(int *)fasp_mem_calloc(nnz,sizeof(int));	
	B->J=(int *)fasp_mem_calloc(nnz,sizeof(int));	
	B->val=(double *)fasp_mem_calloc(nnz,sizeof(double));
	
	for (i=0;i<m;++i) {
		for (j=A->IA[i];j<A->IA[i+1];++j) B->I[N2C(j)]=C2N(i);
	}
	
	memcpy(B->J,  A->JA, nnz*sizeof(int));
	memcpy(B->val,A->val,nnz*sizeof(double));
	
	return status;
}

/**
 * \fn int fasp_format_dstr_dcsr(dSTRmat *A, dCSRmat *B_ptr)
 * \brief Transfer a 'dSTRmat' type matrix into a 'dCSRmat' type matrix.
 *
 * \param *A      pointer to a 'dSTRmat' type matrix
 * \param *B_ptr  pointer to the 'dCSRmat' type matrix
 * 
 * \author Zhou Zhiyang
 * \date 2010/04/29
 */
int fasp_format_dstr_dcsr (dSTRmat *A, 
                           dCSRmat *B_ptr)
{
	//! some members of A
	const int nc    = A->nc;   
	const int ngrid = A->ngrid;
	const int nband = A->nband;
	const int *offsets = A->offsets;
	
	double  *diag = A->diag;
	double **offdiag = A->offdiag;
	
	//! some members of B
	const int glo_row = nc*ngrid;
	int glo_nnz;
	int *ia = NULL;
	int *ja = NULL;
	double *a = NULL;
	
	dCSRmat B;
	
	//! local variables
	int width;
	int nc2 = nc*nc;
	int BAND,ROW,COL;
	int ncb,nci;
	int row_start,col_start;
	int block; //! how many blocks in the current ROW
	int i,j;
	int pos;
	int start;
	int val_L_start,val_R_start;
	int row;
	int tmp_col;
	double tmp_val;
	
	//! allocate for 'ia' array
	ia = (int *)fasp_mem_calloc(glo_row+1,sizeof(int));
	
	//! Generate the 'ia' array
	ia[0] = 0;
	for (ROW = 0; ROW < ngrid; ++ROW)
	{
		block = 1; // diagonal block
		for (BAND = 0; BAND < nband; ++BAND)
		{
			width = offsets[BAND];
			COL   = ROW + width;
			if (width < 0)
			{
				if (COL >= 0) ++block;
			}
			else
			{
				if (COL < ngrid) ++block;
			}
		} //! end for BAND
		
		ncb = nc*block;
		row_start = ROW*nc;
		
		for (i = 0; i < nc; i ++)
		{
			row = row_start + i; 
			ia[row+1] = ia[row] + ncb;
		}
	} //! end for ROW
	
	//! allocate for 'ja' and 'a' arrays
	glo_nnz = ia[glo_row];
	ja = (int *)fasp_mem_calloc(glo_nnz,sizeof(int));
	a = (double *)fasp_mem_calloc(glo_nnz,sizeof(double));
	
	//! Generate the 'ja' and 'a' arrays at the same time 
	for (ROW = 0; ROW < ngrid; ++ROW)
	{
		row_start = ROW*nc;    
		val_L_start = ROW*nc2;
		
		//! deal with the diagonal band
		for (i = 0; i < nc; i ++)
		{
			nci   = nc*i;
			row   = row_start + i;
			start = ia[row];
			for (j = 0; j < nc; j ++)
			{
				pos     = start + j;
				ja[pos] = row_start + j;
				a[pos]  = diag[val_L_start+nci+j];
			}
		}
		block = 1;
		
		//! deal with the off-diagonal bands
		for (BAND = 0; BAND < nband; ++BAND)
		{
			width     = offsets[BAND];
			COL       = ROW + width;
			ncb       = nc*block;
			col_start = COL*nc;
			
			if (width < 0)
			{
				if (COL >= 0)
				{
					val_R_start = COL*nc2;
					for (i = 0; i < nc; i ++)
					{
						nci = nc*i;
						row = row_start + i;
						start = ia[row];
						for (j = 0 ; j < nc; j ++)
						{
							pos     = start + ncb + j;
							ja[pos] = col_start + j;
							a[pos]  = offdiag[BAND][val_R_start+nci+j];
						}
					}
					++block;
				}
			}
			else
			{
				if (COL < ngrid)
				{
					for (i = 0; i < nc; i ++)
					{
						nci = nc*i;
						row = row_start + i;
						start = ia[row];
						for (j = 0; j < nc; j ++)
						{
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
	
	//! Reordering in such manner that every diagonal element 
	//! is firstly stored in the corresponding row
	if (nc > 1)
	{
		for (ROW = 0; ROW < ngrid; ++ROW)
		{
			row_start = ROW*nc;
			for (j = 1; j < nc; j ++)
			{
				row   = row_start + j;
				start = ia[row];
				pos   = start + j;
				
				//! swap in 'ja'
				tmp_col   = ja[start];
				ja[start] = ja[pos];
				ja[pos]   = tmp_col;
				
				//! swap in 'a'
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
	
	return (0);
}

/**
 * \fn dCSRmat fasp_format_bdcsr_dcsr(block_dCSRmat *Ab)
 * \brief Form the whole dCSRmat A using blocks given in Ab
 *
 * \param *Ab   pointer to the blocks
 * \return      the whole dCSRmat matrix if succeed, NULL if fail
 * 
 * \author Shiquan Zhang
 * \date 08/10/2010
 */
dCSRmat fasp_format_bdcsr_dcsr (block_dCSRmat *Ab)
{
    int m=0,n=0,nnz=0;
    const int mb=Ab->brow, nb=Ab->bcol, nbl=mb*nb;
    dCSRmat **blockptr=Ab->blocks, *blockptrij, A;
	
    int i,j,ij,ir,i1,length,ilength,start,irmrow,irmrow1;
    int *row, *col;
	
    row = (int *)fasp_mem_calloc(mb+1,sizeof(int));
    col = (int *)fasp_mem_calloc(nb+1,sizeof(int));
	
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
						memcpy(&(A.val[start]),&(blockptrij->val[blockptrij->IA[irmrow]]),ilength*sizeof(double));
						memcpy(&(A.JA[start]), &(blockptrij->JA[blockptrij->IA[irmrow]]), ilength*sizeof(int));
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
 * \brief Convert a dCSRmat into a dCSRLmat
 * \param *A pointer to the dCSRLmat type matrix
 * \author Zhou Zhiyang
 * \date 2011/01/07
 */
dCSRLmat * fasp_format_dcsrl_dcsr (dCSRmat *A)
{
	double *DATA         = A -> val;
	int    *IA           = A -> IA;
	int    *JA           = A -> JA;
	int     num_rows     = A -> row;
	int     num_cols     = A -> col;
	int     num_nonzeros = A -> nnz;
	
	dCSRLmat *B        = NULL;
	int       dif;
	int      *nzdifnum = NULL;
	int      *rowstart = NULL;
	int      *rowindex = (int *)fasp_mem_calloc(num_rows, sizeof(int));
	int      *ja       = (int *)fasp_mem_calloc(num_nonzeros, sizeof(int));
	double   *data     = (double *)fasp_mem_calloc(num_nonzeros, sizeof(double));
	
	/* auxiliary arrays */
	int *nzrow    = (int *)fasp_mem_calloc(num_rows, sizeof(int));
	int *counter  = NULL;
	int *invnzdif = NULL;
	
	int i,j,k,cnt,maxnzrow;
	
	
	//-----------------------------------------
	//  Generate 'nzrow' and 'maxnzrow'
	//-----------------------------------------
	
	maxnzrow = 0;
	for (i = 0; i < num_rows; i ++)
	{
		nzrow[i] = IA[i+1] - IA[i];
		if (nzrow[i] > maxnzrow) 
		{
			maxnzrow = nzrow[i];
		}
	}
	/* generate 'counter' */
	counter = (int *)fasp_mem_calloc(maxnzrow + 1, sizeof(int));
	for (i = 0; i < num_rows; i ++)
	{
		counter[nzrow[i]] ++;
	}
	
	//--------------------------------------------
	//  Determine 'dif'
	//-------------------------------------------- 
	
	for (dif = 0, i = 0; i < maxnzrow + 1; i ++)
	{
		if (counter[i] > 0) dif ++;
	}
	
	//--------------------------------------------
	//  Generate the 'nzdifnum' and 'rowstart'
	//-------------------------------------------- 
	
	nzdifnum = (int *)fasp_mem_calloc(dif, sizeof(int));
	invnzdif = (int *)fasp_mem_calloc(maxnzrow + 1, sizeof(int));
	rowstart = (int *)fasp_mem_calloc(dif + 1, sizeof(int));
	rowstart[0] = 0;
	for (cnt = 0, i = 0; i < maxnzrow + 1; i ++)
	{
		if (counter[i] > 0) 
		{
			nzdifnum[cnt] = i;
			invnzdif[i] = cnt;
			rowstart[cnt+1] = rowstart[cnt] + counter[i];
			cnt ++;
		}
	}
	
	//--------------------------------------------
	//  Generate the 'rowindex'
	//-------------------------------------------- 
	
	for (i = 0; i < num_rows; i ++)
	{
		j = invnzdif[nzrow[i]];
		rowindex[rowstart[j]] = i;
		rowstart[j] ++;
	}
	/* recover 'rowstart' */
	for (i = dif; i > 0; i --)
	{
		rowstart[i] = rowstart[i-1];
	}
	rowstart[0] = 0;
	
	//--------------------------------------------
	//  Generate the 'data' and 'ja'
	//-------------------------------------------- 
	
	for (cnt = 0, i = 0; i < num_rows; i ++)
	{
		k = rowindex[i];
		for (j = IA[k]; j < IA[k+1]; j ++)
		{
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
	B -> nzdifnum = nzdifnum;
	B -> rowindex = rowindex;
	B -> rowstart = rowstart;
	B -> ja       = ja;
	B -> data     = data;   
	
	//----------------------------
	//  Free the auxiliary arrays
	//----------------------------  
	
	free(nzrow);
	free(counter);
	free(invnzdif);
	
	return B;
}

/*!
 * \fn dCSRmat fasp_format_dbsr_dcsr ( dBSRmat *B )
 * \brief Transfer a 'dBSRmat' type matrix into a dCSRmat.
 * \param dBSRmat *B pointer to the 'dBSRmat' type matrix
 * \author Zhou Zhiyang
 * \date 2010/10/23 
 *
 * \note: works for general nb (Xiaozhe)
 */
dCSRmat fasp_format_dbsr_dcsr (dBSRmat *B)
{
	/* members of B */
	int     ROW = B->ROW;
	int     COL = B->COL;
	int     NNZ = B->NNZ;    
	int     nb  = B->nb;
	int    *IA  = B->IA;
	int    *JA  = B->JA;
	int     storage_manner = B->storage_manner;
	double *val = B->val;
	
	int jump = nb*nb;
	int rowA = ROW*nb;
	int colA = COL*nb;
	int nzA  = NNZ*jump;
	
	dCSRmat A;
	int     *ia = NULL;
	int     *ja = NULL;
	double  *a  = NULL;
	
	int i,j,k;
	int mr,mc;
	int rowstart0,rowstart,colstart0,colstart;
	int colblock,nzperrow; 
	double *vp = NULL;
	double *ap = NULL;
	int    *jap = NULL;
    
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
	
	for (i = rowA; i > 0; i --)
	{
		ia[i] = ia[i-1];
	}
	ia[0] = 0; 
	
	return (A);   
}

/*!
 * \fn dBSRmat fasp_format_dstr_dbsr ( dSTRmat *B )
 * \brief Transfer a 'dSTRmat' type matrix to a 'dBSRmat' type matrix.
 * \param dSTRmat *B pointer to the 'dSTRmat' type matrix
 * \author Zhou Zhiyang
 * \date 2010/10/26 
 */
dBSRmat fasp_format_dstr_dbsr (dSTRmat *B)
{
	// members of 'B'
	int      nc      = B->nc;
	int      ngrid   = B->ngrid;
	double  *diag    = B->diag;
	int      nband   = B->nband;
	int     *offsets = B->offsets;
	double **offdiag = B->offdiag;
	
	// members of 'A'
	dBSRmat A;
	int      NNZ;
	int     *IA  = NULL;
	int     *JA  = NULL;
	double  *val = NULL;
	
	// local variables
	int i,j,k,m;
	int nc2 = nc*nc;
	int ngridplus1 = ngrid + 1;
	
	// compute NNZ
	NNZ = ngrid;
	for (i = 0; i < nband; ++i)
	{
		NNZ += (ngrid - abs(offsets[i]));
	} 
	
	// Create and Initialize a dBSRmat 'A'
	A = fasp_dbsr_create(ngrid, ngrid, NNZ, nc, 0);
	IA = A.IA;
	JA = A.JA;
	val = A.val;
	
	// Generate 'IA'
	for (i = 1; i < ngridplus1; ++i) IA[i] = 1; // take the diagonal blocks into account
	for (i = 0; i < nband; ++i)
	{
		k = offsets[i];
		if (k < 0)
		{
			for (j = -k+1; j < ngridplus1; ++j)
			{
				IA[j] ++;
			}
		}
		else
		{
			m = ngridplus1 - k;
			for (j = 1; j < m; ++j)
			{
				IA[j] ++;
			}
		}
	}
	IA[0] = 0;
	for (i = 1; i < ngridplus1; ++i)
	{
		IA[i] += IA[i-1];
	}
	
	// Generate 'JA' and 'val' at the same time
	for (i = 0 ; i < ngrid; ++i)
	{
		memcpy(val + IA[i]*nc2, diag + i*nc2, nc2*sizeof(double));
		JA[IA[i]] = i;
		IA[i] ++;
	}
	
	for (i = 0; i < nband; ++i)
	{
		k = offsets[i];
		if (k < 0)
		{
			for (j = -k; j < ngrid; ++j)
			{
				m = j + k;
				memcpy(val+IA[j]*nc2, offdiag[i]+m*nc2, nc2*sizeof(double));
				JA[IA[j]] = m;
				IA[j] ++;
			}
		}
		else
		{
			m = ngrid - k;
			for (j = 0; j < m; ++j)
			{
				memcpy(val + IA[j]*nc2, offdiag[i] + j*nc2, nc2*sizeof(double));
				JA[IA[j]] = k + j;
				IA[j] ++;
			}
		}
	}
	
	// Map back the real 'IA' for BSR of A
	for (i = ngrid; i > 0; i --)
	{
		IA[i] = IA[i-1];
	}
	IA[0] = 0; 
	
	return (A);
}


/*!
 * \fn dCOOmat * fasp_format_dbsr_dcoo ( dBSRmat *B )
 * \brief Transfer a 'dBSRmat' type matrix to a 'dCOOmat' type matrix.
 * \param dSTRmat *B pointer to the 'dBSRmat' type matrix
 * \author Zhou Zhiyang
 * \date 2010/10/26 
 */
dCOOmat * fasp_format_dbsr_dcoo (dBSRmat *B)
{
	/* members of B */
	int     ROW = B->ROW;
	int     COL = B->COL;
	int     NNZ = B->NNZ;    
	int     nb  = B->nb;
	int    *IA  = B->IA;
	int    *JA  = B->JA;
	double *val = B->val;
	
	dCOOmat *A = NULL;
	int nb2 = nb*nb;
	int num_nonzeros = NNZ*nb2;
	int *rowA = NULL;
	int *colA = NULL;
	double *valA = NULL;
	
	int i,j,k,inb;
	int row_start, col_start;
	int cnt,mr,mc;
	double *pt = NULL;
	
	// Create and Initialize a dBSRmat 'A'
	A = (dCOOmat *)fasp_mem_calloc(1, sizeof(dCOOmat));
	A->row = ROW*nb;
	A->col = COL*nb;
	A->nnz = num_nonzeros;    
	rowA   = (int *)fasp_mem_calloc(num_nonzeros, sizeof(int)); 
	colA   = (int *)fasp_mem_calloc(num_nonzeros, sizeof(int));
	valA   = (double *)fasp_mem_calloc(num_nonzeros, sizeof(double));
	A->I   = rowA;
	A->J   = colA;
	A->val = valA;
	
	cnt = 0;
	for (i = 0; i < ROW; ++i)
	{
		inb = i*nb;
		for (k = IA[i]; k < IA[i+1]; ++k)
		{
			j  = JA[k];
			pt = &val[k*nb2];
			row_start = inb; 
			col_start = j*nb;
			for (mr = 0; mr < nb; mr ++)
			{
				for (mc = 0; mc < nb; mc ++)
				{
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

/*-----------------------------------omp--------------------------------------*/

#if FASP_USE_OPENMP

/*!
 * \fn dCSRmat dBSR2dCSRMatrix_omp( dBSRmat *B, int nthreads, int openmp_holds )
 * \brief Transfer a 'dBSRmat' type matrix into a dCSRmat.
 * \param dBSRmat *B pointer to the 'dBSRmat' type matrix
 * \param nthreads number of threads
 * \param openmp_holds threshold of parallelization
 *
 * \author Feng Chunsheng, Yue Xiaoqiang
 * \date 03/01/2011
 */
dCSRmat fasp_format_dbsr_dcsr_omp (dBSRmat *B, 
                                   int nthreads, 
                                   int openmp_holds )
{
	dCSRmat A;
    
	/* members of B */
	int     ROW = B->ROW;
	int     COL = B->COL;
	int     NNZ = B->NNZ;    
	int     nb  = B->nb;
	int    *IA  = B->IA;
	int    *JA  = B->JA;
	int     storage_manner = B->storage_manner;
	double *val = B->val;
	
	int jump = nb*nb;
	int rowA = ROW*nb;
	int colA = COL*nb;
	int nzA  = NNZ*jump;
	
	int     *ia = NULL;
	int     *ja = NULL;
	double  *a  = NULL;
	
	int i,j,k;
	int mr,mc;
	int rowstart0,rowstart,colstart0,colstart;
	int colblock,nzperrow; 
	double *vp = NULL;
	double *ap = NULL;
	int    *jap = NULL;
	int myid, mybegin, myend;
    
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
	
	if (ROW > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i, rowstart, colblock, nzperrow, j) ////num_threads(nthreads)
		{
			myid = omp_get_thread_num();
			FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
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
			if (ROW > openmp_holds) {
#pragma omp parallel private(myid, mybegin, myend, i, k, j, rowstart, colstart, vp, mr, ap, jap, mc) ////num_threads(nthreads)
				{
					myid = omp_get_thread_num();
					FASP_GET_START_END(myid, nthreads, ROW, mybegin, myend);
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
	
	for (i = rowA; i > 0; i --)
	{
		ia[i] = ia[i-1];
	}
	ia[0] = 0; 
    
	return (A);   
}
#endif // OMP

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
