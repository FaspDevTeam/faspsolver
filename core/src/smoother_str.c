/*! \file smoother_str.c
 *  \brief Smoothers for sparse matrix in STR format
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

static void blkcontr2_str(INT start_data, INT start_vecx, INT start_vecy, INT nc, REAL *data, REAL *x, REAL *y);
static void saAxpby(REAL alpha, REAL beta, INT size, REAL *A, REAL *x, REAL *y);

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dstr_jacobi (dSTRmat *A, dvector *b, dvector *u)
 *
 * \brief Jacobi method as the smoother
 *
 * \param *A   pointer to stiffness matrix
 * \param *b   pointer to right hand side vector
 * \param *u   initial guess and new approximation to the solution obtained after one iteration 
 */
void fasp_smoother_dstr_jacobi (dSTRmat *A, 
                                dvector *b, 
                                dvector *u)
{
	INT     nc    = A->nc;    // size of each block (number of components)
	INT     ngrid = A->ngrid; // number of grids
	REAL   *diag  = A->diag;  // Diagonal entries		
	REAL   *diaginv = NULL;   // Diagonal inverse, same size and storage scheme as A->diag
	
	INT nc2   = nc*nc;
	INT size  = nc2*ngrid;
	INT block = 0;
	INT start = 0;
    
	if (nc > 1)
	{
		// allocate memory
		diaginv = (REAL *)fasp_mem_calloc(size,sizeof(REAL));
		
		// diaginv = diag;
		fasp_array_cp(size,diag,diaginv);
		
		// generate diaginv
		for (block = 0; block < ngrid; block ++)
		{          
			fasp_blas_smat_inv(diaginv+start, nc);
			start += nc2;
		}
	}
	
	fasp_smoother_dstr_jacobi1(A, b, u, diaginv);
	
	fasp_mem_free(diaginv);
}


/**
 * \fn void fasp_smoother_dstr_jacobi1(dSTRmat *A, dvector *b, dvector *u, REAL *diaginv)
 *
 * \brief Jacobi method as the smoother with diag_inv given
 *
 * \param *A        pointer to stiffness matrix
 * \param *b        pointer to right hand side vector
 * \param *u        pointer to approximation solution 
 * \param *diaginv  all the inverse matrices for all the diagonal block of A when (A->nc)>1,
 *                                                                  and NULL when (A->nc)=1  
 */
void fasp_smoother_dstr_jacobi1 (dSTRmat *A, 
                                 dvector *b, 
                                 dvector *u, 
                                 REAL *diaginv)	
{	
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT size  = nc*ngrid; 
	INT nc2   = nc*nc;
	INT start  = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecb = 0;
	INT start_vecu = 0;
	
	// auxiliary array
	REAL *b_tmp = NULL;
	
	// this should be done once and for all!!
	b_tmp = (REAL *)fasp_mem_calloc(size,sizeof(REAL));
	
	// b_tmp = b_val
	fasp_array_cp(size,b_val,b_tmp);	
	
	// It's not necessary to assign the smoothing order since the results doesn't depend on it
	if (nc == 1)
	{
		for (point = 0; point < ngrid; point ++)
		{
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						b_tmp[point] -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						b_tmp[point] -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
		} // end for point
		
		for (point = 0; point < ngrid; point ++)
		{
			// zero-diagonal should be tested previously
			u_val[point] = b_tmp[point] / diag[point]; 
		}
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (block = 0; block < ngrid; block ++)
		{
			start_DATA = nc2*block;
			start_vecb = nc*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, start_vecb,
                                      nc, offdiag[band], u_val, b_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, start_vecb,
                                      nc, offdiag[band], u_val, b_tmp );
					}
				}
			} // end for band     
		} // end for block
		
		for (block = 0; block < ngrid; block ++)
		{
			start = nc*block;
			fasp_blas_smat_mxv(diaginv+nc2*block, b_tmp+start, u_val+start, nc);
		}
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(b_tmp);
}

/**
 * \fn void fasp_smoother_dstr_gs (dSTRmat *A, dvector *b, dvector *u, INT order, INT *mark)
 *
 * \brief Gauss-Seidel method as the smoother
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending manner
 *           DESCEND      21: in descending manner
 *        when mark != NULL
 *           USERDEFINED  0 : in the user-defined manner
 *           CPFIRST      1 : C-points first and then F-points
 *           FPFIRST     -1 : F-points first and then C-points
 * \param *mark pointer to the user-defined ordering(when order=0) 
 *        or CF_marker array(when order!=0)
 */
void fasp_smoother_dstr_gs (dSTRmat *A, 
                            dvector *b, 
                            dvector *u,
                            INT order, 
                            INT *mark)
{
	INT     nc    = A->nc;    // size of each block (number of components)
	INT     ngrid = A->ngrid; // number of grids
	REAL *diag  = A->diag;  // Diagonal entries		
	REAL *diaginv = NULL;   // Diagonal inverse(when nc>1),same size and storage scheme as A->diag
	
	INT nc2   = nc*nc;
	INT size  = nc2*ngrid;
	INT block = 0;
	INT start = 0;
    
	if (nc > 1)
	{
		// allocate memory
		diaginv = (REAL *)fasp_mem_calloc(size,sizeof(REAL));
		
		// diaginv = diag;
		fasp_array_cp(size,diag,diaginv);
		
		// generate diaginv
		for (block = 0; block < ngrid; block ++)
		{    
			fasp_blas_smat_inv(diaginv+start, nc);
			start += nc2;
		}
	}
	
	fasp_smoother_dstr_gs1(A, b, u, order, mark, diaginv);
	
	fasp_mem_free(diaginv);
}

/**
 * \fn void fasp_smoother_dstr_gs1(dSTRmat *A, dvector *b, dvector *u, INT order, INT *mark, REAL *diaginv)
 *
 * \brief Gauss-Seidel method as the smoother with diag_inv given
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending manner
 *           DESCEND      21: in descending manner
 *        when mark != NULL
 *           USERDEFINED  0 : in the user-defined manner
 *           CPFIRST      1 : C-points first and then F-points
 *           FPFIRST     -1 : F-points first and then C-points
 * \param *mark pointer to the user-defined ordering or CF_marker array
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1  
 */
void fasp_smoother_dstr_gs1 (dSTRmat *A, 
                             dvector *b, 
                             dvector *u, 
                             INT order, 
                             INT *mark, 
                             REAL *diaginv)
{	
	
	if (!mark)
	{
		if (order == ASCEND)       // smooth ascendingly
		{
			fasp_smoother_dstr_gs_ascend(A, b, u, diaginv);
		}
		else if (order == DESCEND) // smooth descendingly
		{
			fasp_smoother_dstr_gs_descend(A, b, u, diaginv);
		}
	}
	else
	{
		if (order == USERDEFINED)  // smooth according to the order 'mark' defined by user
		{
			fasp_smoother_dstr_gs_order(A, b, u, diaginv, mark);
		}
		else // smooth according to 'mark', where 'mark' is a CF_marker array
		{
			fasp_smoother_dstr_gs_cf(A, b, u, diaginv, mark, order);
		}
	}
}

/**
 * \fn void fasp_smoother_dstr_gs_ascend (dSTRmat *A, dvector *b, dvector *u, REAL *diaginv)
 *
 * \brief Gauss-Seidel method as the smoother in the ascending manner
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1 
 */
void fasp_smoother_dstr_gs_ascend (dSTRmat *A, 
                                   dvector *b,
                                   dvector *u, 
                                   REAL *diaginv)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		for (point = 0; point < ngrid; point ++)
		{
			rhs = b_val[point];
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						rhs -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						rhs -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
			
			// zero-diagonal should be tested previously
			u_val[point] = rhs / diag[point];	          
			
		} // end for point
		
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (block = 0; block < ngrid; block ++)
		{
			ncb = nc*block;
			for (point = 0; point < nc; point ++)
			{
				vec_tmp[point] = b_val[ncb+point];
			}
			start_DATA = nc2*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
			} // end for band
			
			// subblock smoothing
			fasp_blas_smat_mxv(diaginv+start_DATA, vec_tmp, u_val+nc*block, nc);
			
		} // end for block
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_gs_descend (dSTRmat *A, dvector *b, dvector *u, REAL *diaginv)
 *
 * \brief Gauss-Seidel method as the smoother in the descending manner
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1 
 */
void fasp_smoother_dstr_gs_descend (dSTRmat *A, 
                                    dvector *b, 
                                    dvector *u, 
                                    REAL *diaginv)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		for (point = ngrid-1; point >= 0; point --)
		{
			rhs = b_val[point];
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						rhs -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						rhs -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
			
			// zero-diagonal should be tested previously
			u_val[point] = rhs / diag[point];	          
			
		} // end for point
		
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (block = ngrid-1; block >= 0; block --)
		{
			ncb = nc*block;
			for (point = 0; point < nc; point ++)
			{
				vec_tmp[point] = b_val[ncb+point];
			}
			start_DATA = nc2*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
			} // end for band
			
			// subblock smoothing
			fasp_blas_smat_mxv(diaginv+start_DATA, vec_tmp, u_val+nc*block, nc);
			
		} // end for block
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_gs_order(dSTRmat *A, dvector *b, dvector *u, REAL *diaginv, INT *mark)
 *
 * \brief Gauss method as the smoother in the user-defined order
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1
 * \param *mark pointer to the user-defined order array 
 */
void fasp_smoother_dstr_gs_order (dSTRmat *A, 
                                  dvector *b, 
                                  dvector *u,
                                  REAL *diaginv, 
                                  INT *mark)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT index = 0; 	
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		for (index = 0; index < ngrid; index ++)
		{
			point = mark[index];
			rhs = b_val[point];
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						rhs -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						rhs -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
			
			// zero-diagonal should be tested previously
			u_val[point] = rhs / diag[point];	          
			
		} // end for index
		
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (index = 0; index < ngrid; index ++)
		{
			block = mark[index];
			ncb = nc*block;
			for (point = 0; point < nc; point ++)
			{
				vec_tmp[point] = b_val[ncb+point];
			}
			start_DATA = nc2*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
			} // end for band
			
			// subblock smoothing
			fasp_blas_smat_mxv(diaginv+start_DATA, vec_tmp, u_val+nc*block, nc);
			
		} // end for index
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_gs_cf (dSTRmat *A, dvector *b, dvector *u, REAL *diaginv, INT *mark, INT order)
 *
 * \brief Gauss method as the smoother in the C-F manner
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1
 * \param *mark pointer to the user-defined order array 
 * \param order a flag to indicate the order for smoothing
 *           CPFIRST  1 : C-points first and then F-points
 *           FPFIRST -1 : F-points first and then C-points   
 */
void fasp_smoother_dstr_gs_cf (dSTRmat *A, 
                               dvector *b, 
                               dvector *u, 
                               REAL *diaginv, 
                               INT *mark, 
                               INT order)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	INT FIRST  = order;  // which kind of points to be smoothed firstly?
	INT SECOND = -order; // which kind of points to be smoothed secondly?
	
	REAL rhs = 0.0;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		// deal with the points marked FIRST
		for (point = 0; point < ngrid; point ++)
		{
			if (mark[point] == FIRST)
			{
				rhs = b_val[point];
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = point + width;
					if (width < 0)
					{
						if (column >= 0)
						{
							rhs -= offdiag[band][column]*u_val[column];
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							rhs -= offdiag[band][point]*u_val[column];
						}
					}
				} // end for band
				
				// zero-diagonal should be tested previously
				u_val[point] = rhs / diag[point];	          
			} // end if (mark[point] == FIRST)
		} // end for point
		
		// deal with the points marked SECOND
		for (point = 0; point < ngrid; point ++)
		{
			if (mark[point] == SECOND)
			{
				rhs = b_val[point];
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = point + width;
					if (width < 0)
					{
						if (column >= 0)
						{
							rhs -= offdiag[band][column]*u_val[column];
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							rhs -= offdiag[band][point]*u_val[column];
						}
					}
				} // end for band
				
				// zero-diagonal should be tested previously
				u_val[point] = rhs / diag[point];	          
			} // end if (mark[point] == SECOND)
		} // end for point
		
	} // end if (nc == 1)  
	
	else if (nc > 1)
	{
		// deal with the blocks marked FIRST
		for (block = 0; block < ngrid; block ++)
		{
			if (mark[block] == FIRST)
			{
				ncb = nc*block;
				for (point = 0; point < nc; point ++)
				{
					vec_tmp[point] = b_val[ncb+point];
				}
				start_DATA = nc2*block;
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = block + width; 
					if (width < 0)
					{
						if (column >= 0)
						{
							start_data = nc2*column;
							start_vecu = nc*column;
							blkcontr2_str( start_data, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							start_vecu = nc*column;
							blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
				} // end for band
				
				// subblock smoothing
				fasp_blas_smat_mxv(diaginv+start_DATA, vec_tmp, u_val+nc*block, nc);
			} // end if (mark[block] == FIRST)
			
		} // end for block
		
		// deal with the blocks marked SECOND
		for (block = 0; block < ngrid; block ++)
		{
			if (mark[block] == SECOND)
			{
				ncb = nc*block;
				for (point = 0; point < nc; point ++)
				{
					vec_tmp[point] = b_val[ncb+point];
				}
				start_DATA = nc2*block;
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = block + width; 
					if (width < 0)
					{
						if (column >= 0)
						{
							start_data = nc2*column;
							start_vecu = nc*column;
							blkcontr2_str( start_data, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							start_vecu = nc*column;
							blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
				} // end for band
				
				// subblock smoothing
				fasp_blas_smat_mxv(diaginv+start_DATA, vec_tmp, u_val+nc*block, nc);
			} // end if (mark[block] == SECOND)
			
		} // end for block
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_sor (dSTRmat *A, dvector *b, dvector *u, INT order, INT *mark, REAL weight)
 *
 * \brief SOR method as the smoother
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending manner
 *           DESCEND      21: in descending manner
 *        when mark != NULL
 *           USERDEFINED  0 : in the user-defined manner
 *           CPFIRST      1 : C-points first and then F-points
 *           FPFIRST     -1 : F-points first and then C-points
 * \param *mark pointer to the user-defined ordering(when order=0) 
 *        or CF_marker array(when order!=0)
 * \param weight over-relaxation parameter
 */
void fasp_smoother_dstr_sor (dSTRmat *A, 
                             dvector *b, 
                             dvector *u, 
                             INT order, 
                             INT *mark, 
                             REAL weight)
{
	INT     nc    = A->nc;    // size of each block (number of components)
	INT     ngrid = A->ngrid; // number of grids
	REAL *diag  = A->diag;  // Diagonal entries		
	REAL *diaginv = NULL;   // Diagonal inverse(when nc>1),same size and storage scheme as A->diag
	
	INT nc2   = nc*nc;
	INT size  = nc2*ngrid;
	INT block = 0;
	INT start = 0;
	
	if (nc > 1)
	{
		// allocate memory
		diaginv = (REAL *)fasp_mem_calloc(size,sizeof(REAL));
		
		// diaginv = diag;
		fasp_array_cp(size,diag,diaginv);
		
		// generate diaginv
		for (block = 0; block < ngrid; block ++)
		{    
			fasp_blas_smat_inv(diaginv+start, nc);
			start += nc2;
		}
	}
	
	fasp_smoother_dstr_sor1(A, b, u, order, mark, diaginv, weight);
	
	fasp_mem_free(diaginv);
}

/**
 * \fn void fasp_smoother_dstr_sor1(dSTRmat *A, dvector *b, dvector *u, INT order, 
 *                     INT *mark, REAL *diaginv, REAL weight)
 *
 * \brief SOR method as the smoother
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param order a flag to indicate the order for smoothing
 *        when mark = NULL       
 *           ASCEND       12: in ascending manner
 *           DESCEND      21: in descending manner
 *        when mark != NULL
 *           USERDEFINED  0 : in the user-defined manner
 *           CPFIRST      1 : C-points first and then F-points
 *           FPFIRST     -1 : F-points first and then C-points
 * \param *mark pointer to the user-defined ordering or CF_marker array
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1
 * \param weight over-relaxation parameter   
 */
void fasp_smoother_dstr_sor1(dSTRmat *A, 
                             dvector *b,
                             dvector *u, 
                             INT order, 
                             INT *mark, 
                             REAL *diaginv, 
                             REAL weight)
{	
	if (!mark)
	{
		if (order == ASCEND)       // smooth ascendingly
		{
			fasp_smoother_dstr_sor_ascend(A, b, u, diaginv, weight);
		}
		else if (order == DESCEND) // smooth descendingly
		{
			fasp_smoother_dstr_sor_descend(A, b, u, diaginv, weight);
		}
	}
	else
	{
		if (order == USERDEFINED)  // smooth according to the order 'mark' defined by user
		{
			fasp_smoother_dstr_sor_order(A, b, u, diaginv, mark, weight);
		}
		else // smooth according to 'mark', where 'mark' is a CF_marker array
		{
			fasp_smoother_dstr_sor_cf(A, b, u, diaginv, mark, order, weight);
		}
	}
}

/**
 * \fn void fasp_smoother_dstr_sor_ascend(dSTRmat *A, dvector *b, dvector *u, REAL *diaginv, REAL weight)
 *
 * \brief SOR method as the smoother in the ascending manner
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1 
 * \param weight over-relaxation parameter  
 */
void fasp_smoother_dstr_sor_ascend(dSTRmat *A,
                                   dvector *b,
                                   dvector *u,
                                   REAL *diaginv, 
                                   REAL weight)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	REAL one_minus_weight = 1.0 - weight;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		for (point = 0; point < ngrid; point ++)
		{
			rhs = b_val[point];
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						rhs -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						rhs -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
			
			// zero-diagonal should be tested previously
			u_val[point] = one_minus_weight*u_val[point] + 
			weight*(rhs / diag[point]);	          
			
		} // end for point
		
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (block = 0; block < ngrid; block ++)
		{
			ncb = nc*block;
			for (point = 0; point < nc; point ++)
			{
				vec_tmp[point] = b_val[ncb+point];
			}
			start_DATA = nc2*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
			} // end for band
			
			// subblock smoothing
			saAxpby(weight, one_minus_weight, nc, 
                    diaginv+start_DATA, vec_tmp, u_val+nc*block);
			
		} // end for block
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_sor_descend (dSTRmat *A, dvector *b, dvector *u, REAL *diaginv, REAL weight)
 *
 * \brief SOR method as the smoother in the descending manner
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1 
 * \param weight over-relaxation parameter  
 */
void fasp_smoother_dstr_sor_descend (dSTRmat *A, 
                                     dvector *b, 
                                     dvector *u, 
                                     REAL *diaginv, 
                                     REAL weight)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	REAL one_minus_weight = 1.0 - weight;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		for (point = ngrid-1; point >= 0; point --)
		{
			rhs = b_val[point];
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						rhs -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						rhs -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
			
			// zero-diagonal should be tested previously
			u_val[point] = one_minus_weight*u_val[point] + 
			weight*(rhs / diag[point]);	          
			
		} // end for point
		
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (block = ngrid-1; block >= 0; block --)
		{
			ncb = nc*block;
			for (point = 0; point < nc; point ++)
			{
				vec_tmp[point] = b_val[ncb+point];
			}
			start_DATA = nc2*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
			} // end for band
			
			// subblock smoothing
			saAxpby(weight, one_minus_weight, nc, 
                    diaginv+start_DATA, vec_tmp, u_val+nc*block);
			
		} // end for block
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_sor_order (dSTRmat *A, dvector *b, dvector *u, 
 *                             REAL *diaginv, INT *mark, REAL weight)
 *
 * \brief SOR method as the smoother in the user-defined order
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1
 * \param *mark pointer to the user-defined order array   
 * \param weight over-relaxation parameter  
 */
void fasp_smoother_dstr_sor_order (dSTRmat *A, 
                                   dvector *b,
                                   dvector *u, 
                                   REAL *diaginv, 
                                   INT *mark, 
                                   REAL weight)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT index  = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	REAL one_minus_weight = 1.0 - weight;
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		for (index = 0; index < ngrid; index ++)
		{
			point = mark[index];
			rhs = b_val[point];
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = point + width;
				if (width < 0)
				{
					if (column >= 0)
					{
						rhs -= offdiag[band][column]*u_val[column];
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						rhs -= offdiag[band][point]*u_val[column];
					}
				}
			} // end for band
			
			// zero-diagonal should be tested previously
			u_val[point] = one_minus_weight*u_val[point] + 
			weight*(rhs / diag[point]);	          
			
		} // end for index
		
	} // end if (nc == 1)  
	else if (nc > 1)
	{
		for (index = 0; index < ngrid; index ++)
		{
			block = mark[index];
			ncb = nc*block;
			for (point = 0; point < nc; point ++)
			{
				vec_tmp[point] = b_val[ncb+point];
			}
			start_DATA = nc2*block;
			for (band = 0; band < nband; band ++)
			{
				width  = offsets[band];
				column = block + width; 
				if (width < 0)
				{
					if (column >= 0)
					{
						start_data = nc2*column;
						start_vecu = nc*column;
						blkcontr2_str( start_data, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
				else // width > 0
				{
					if (column < ngrid)
					{
						start_vecu = nc*column;
						blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                      offdiag[band], u_val, vec_tmp );
					}
				}
			} // end for band
			
			// subblock smoothing
			saAxpby(weight, one_minus_weight, nc, 
                    diaginv+start_DATA, vec_tmp, u_val+nc*block);
			
		} // end for index
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_smoother_dstr_sor_cf (dSTRmat *A, dvector *b, dvector *u, REAL *diaginv, 
 *                        INT *mark, INT order, REAL weight)
 *
 * \brief SOR method as the smoother in the C-F manner
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration
 * \param *diaginv all the inverse matrices for all the diagonal block of A when (A->nc) > 1,
 *        and NULL when (A->nc)=1
 * \param *mark pointer to the user-defined order array 
 * \param order a flag to indicate the order for smoothing
 *           CPFIRST  1 : C-points first and then F-points
 *           FPFIRST -1 : F-points first and then C-points  
 * \param weight over-relaxation parameter  
 */
void fasp_smoother_dstr_sor_cf (dSTRmat *A, 
                                dvector *b, 
                                dvector *u, 
                                REAL *diaginv,
                                INT *mark, 
                                INT order, 
                                REAL weight)
{
	// information of A
	INT ngrid = A->ngrid;  // number of grids
	INT nc = A->nc;        // size of each block (number of components)
	INT nband = A->nband ; // number of off-diag band
	INT *offsets = A->offsets; // offsets of the off-diagals
	REAL  *diag = A->diag;       // Diagonal entries
	REAL **offdiag = A->offdiag; // Off-diagonal entries	
	
	// values of dvector b and u
	REAL *b_val = b->val;
	REAL *u_val = u->val;
	
	// local variables
	INT block = 0;
	INT point = 0;
	INT band  = 0;
	INT width = 0;
	INT nc2   = nc*nc;
	INT ncb   = 0;
	INT column = 0;
	INT start_data = 0;
	INT start_DATA = 0;
	INT start_vecu = 0;
	REAL rhs = 0.0;
	REAL one_minus_weight = 1.0 - weight;
	INT FIRST  = order;  // which kind of points to be smoothed firstly?
	INT SECOND = -order; // which kind of points to be smoothed secondly?
	
	// auxiliary array(nc*1 vector)
	REAL *vec_tmp = NULL;
	
	vec_tmp = (REAL *)fasp_mem_calloc(nc,sizeof(REAL));
	
	if (nc == 1)
	{
		// deal with the points marked FIRST
		for (point = 0; point < ngrid; point ++)
		{
			if (mark[point] == FIRST)
			{
				rhs = b_val[point];
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = point + width;
					if (width < 0)
					{
						if (column >= 0)
						{
							rhs -= offdiag[band][column]*u_val[column];
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							rhs -= offdiag[band][point]*u_val[column];
						}
					}
				} // end for band
				
				// zero-diagonal should be tested previously
				u_val[point] = one_minus_weight*u_val[point] + 
				weight*(rhs / diag[point]);
				
			} // end if (mark[point] == FIRST)
		} // end for point
		
		// deal with the points marked SECOND
		for (point = 0; point < ngrid; point ++)
		{
			if (mark[point] == SECOND)
			{
				rhs = b_val[point];
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = point + width;
					if (width < 0)
					{
						if (column >= 0)
						{
							rhs -= offdiag[band][column]*u_val[column];
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							rhs -= offdiag[band][point]*u_val[column];
						}
					}
				} // end for band
				
				// zero-diagonal should be tested previously
				u_val[point] = rhs / diag[point];	          
			} // end if (mark[point] == SECOND)
		} // end for point
		
	} // end if (nc == 1)  
	
	else if (nc > 1)
	{
		// deal with the blocks marked FIRST
		for (block = 0; block < ngrid; block ++)
		{
			if (mark[block] == FIRST)
			{
				ncb = nc*block;
				for (point = 0; point < nc; point ++)
				{
					vec_tmp[point] = b_val[ncb+point];
				}
				start_DATA = nc2*block;
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = block + width; 
					if (width < 0)
					{
						if (column >= 0)
						{
							start_data = nc2*column;
							start_vecu = nc*column;
							blkcontr2_str( start_data, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							start_vecu = nc*column;
							blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
				} // end for band
				
				// subblock smoothing
				saAxpby(weight, one_minus_weight, nc, 
                        diaginv+start_DATA, vec_tmp, u_val+nc*block);	            
			} // end if (mark[block] == FIRST)
			
		} // end for block
		
		// deal with the blocks marked SECOND
		for (block = 0; block < ngrid; block ++)
		{
			if (mark[block] == SECOND)
			{
				ncb = nc*block;
				for (point = 0; point < nc; point ++)
				{
					vec_tmp[point] = b_val[ncb+point];
				}
				start_DATA = nc2*block;
				for (band = 0; band < nband; band ++)
				{
					width  = offsets[band];
					column = block + width; 
					if (width < 0)
					{
						if (column >= 0)
						{
							start_data = nc2*column;
							start_vecu = nc*column;
							blkcontr2_str( start_data, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
					else // width > 0
					{
						if (column < ngrid)
						{
							start_vecu = nc*column;
							blkcontr2_str( start_DATA, start_vecu, 0, nc,
                                          offdiag[band], u_val, vec_tmp );
						}
					}
				} // end for band
				
				// subblock smoothing
				saAxpby(weight, one_minus_weight, nc, 
                        diaginv+start_DATA, vec_tmp, u_val+nc*block);
			} // end if (mark[block] == SECOND)
			
		} // end for block
		
	} // end else if (nc > 1)
	else
	{
		printf("\n nc is illegal!\n\n");
		return;
	}
	fasp_mem_free(vec_tmp);
}

/**
 * \fn void fasp_generate_diaginv_block(dSTRmat *A, ivector *neigh, dvector *diaginv, divector *pivot)
 *
 * \brief generate inverse of diagonal block for block smoothers
 *
 * \param *A pointer to stiffness matrix
 * \param *neigh pointer to neighborhoods
 * \param *diaginv pointer to the inverse of the diagonals
 *
 * \author Xiaozhe Hu
 */
void fasp_generate_diaginv_block (dSTRmat *A, 
                                  ivector *neigh, 
                                  dvector *diaginv, 
                                  ivector *pivot)
{	
	// information about A
	const INT nc = A->nc;
	const INT ngrid = A->ngrid;
	const INT nband = A->nband;
	
	INT *offsets = A->offsets;
	REAL *diag = A->diag;
	REAL **offdiag = A->offdiag;
	
	// information about neighbors
	INT nneigh; 
	if (!neigh)
	{
		nneigh = 0;
	}
	else
	{
		nneigh= neigh->row/ngrid;
	}
	
	// local variable
	INT i, j, k, l, m, n, nbd, p;
	INT count;
	INT block_size;
	INT mem_inv = 0;
	INT mem_pivot = 0;
	
	// allocation 
	REAL *temp = (REAL *)fasp_mem_calloc(((nneigh+1)*nc)*((nneigh+1)*nc)*ngrid, sizeof(REAL));
	INT *tmp = (INT *)fasp_mem_calloc(((nneigh+1)*nc)*ngrid, sizeof(INT));
	
	// main loop
	for (i=0; i<ngrid; ++i)
	{
		// count number of neighbors of node i
		count = 1;
		for (l=0; l<nneigh; ++l)
		{
			if (neigh->val[i*nneigh+l] >= 0) count++ ;
		}
		
		// prepare the inverse of diagonal block i
		block_size = count*nc;
		
		diaginv[i].row = block_size*block_size;
		diaginv[i].val = temp + mem_inv;
		mem_inv += diaginv[i].row;
		
		pivot[i].row = block_size;
		pivot[i].val = tmp + mem_pivot;
		mem_pivot += pivot[i].row;
		
		// put the diagonal block corresponding to node i
		for (j=0; j<nc; ++j)
		{
			for (k=0; k<nc; ++k)
			{
				diaginv[i].val[j*block_size+k] = diag[i*nc*nc + j*nc + k];
			}
		}
		
		// put the blocks corresponding to the neighbor of node i
		count = 1;
		for (l=0; l<nneigh; ++l)
		{
			p = neigh->val[i*nneigh+l];
			if (p >= 0){
				// put the diagonal block corresponding to this neighbor
				for (j=0; j<nc; ++j)
				{
					for (k=0; k<nc; ++k)
					{
						m = count*nc + j; n = count*nc+k;
						diaginv[i].val[m*block_size+n] = diag[p*nc*nc+j*nc+k];
					}
				}
				
				for (nbd=0; nbd<nband; nbd++)
				{
					if ( offsets[nbd] == (p-i) )// put the block corresponding to (i, p)
					{
						for (j=0; j<nc; ++j)
						{
							for(k=0; k<nc; ++k)
							{
								m = j; n = count*nc + k;
								diaginv[i].val[m*block_size+n] = offdiag[nbd][(p-MAX(p-i, 0))*nc*nc+j*nc+k];
							}
						}
					}  
					
					if ( offsets[nbd] == (i-p) )// put the block corresponding to (p, i)
					{
						for (j=0; j<nc; ++j)
						{
							for(k=0; k<nc; ++k)
							{
								m = count*nc + j; n = k;
								diaginv[i].val[m*block_size+n] = offdiag[nbd][(i-MAX(i-p, 0))*nc*nc+j*nc+k];
							}
						}
					}
				}
				count++;
			} //end if
		}// end for (l=0; l<nneigh; ++l)
		
		//fasp_blas_smat_inv(diaginv[i].val, block_size);
		fasp_smat_lu_decomp(diaginv[i].val, pivot[i].val, block_size);
		
	} // end of main loop
}

/**
 * \fn void fasp_smoother_dstr_schwarz(dSTRmat *A, dvector *b, dvector *u, dvector *diaginv, ivector *pivot, ivector *neigh, INT *order)
 *
 * \brief Schwarz method as the smoother
 *
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side vector
 * \param *u initial guess and new approximation to the solution obtained after one iteration 
 * \param *diaginv pointer to the inverse of diagonal blocks
 * \param *pivot pointer to the pivot of diagonal blocks
 * \param *neigh pointer to the neighbors
 * \param *order pointer to the smoothing order 
 *
 * \author Xiaozhe Hu
 */
void fasp_smoother_dstr_schwarz (dSTRmat *A, 
                                 dvector *b, 
                                 dvector *u, 
                                 dvector *diaginv, 
                                 ivector *pivot, 
                                 ivector *neigh, 
                                 ivector *order)
{
	// information about A
	const INT ngrid = A->ngrid;
	const INT nc = A->nc;
	
	// information about neighbors
	INT nneigh; 
	if (!neigh)
	{
		nneigh = 0;
	}
	else
	{
		nneigh= neigh->row/ngrid;
	}
	
	// local variable
	INT i, j, k, l, p, ti;
	
	// work space
	REAL *temp = (REAL *)fasp_mem_calloc(b->row + (nneigh+1)*nc + (nneigh+1)*nc, sizeof(REAL));
	dvector r, e, ri;
	r.row = b->row; r.val = temp;
	e.row = (nneigh+1)*nc; e.val = temp + b->row;
	ri.row = (nneigh+1)*nc; ri.val = temp + b->row + (nneigh+1)*nc;
	
	// initial residual 
	fasp_dvec_cp(b,&r);fasp_blas_dstr_aAxpy(-1.0,A,u->val,r.val);
	
	// main loop
	if (!order)
	{
		for (i=0; i<ngrid; ++i)
		{
			//-----------------------------------------------------
			// right hand side for A_ii e_i = r_i
			// rhs corresponding to node i
			for (j=0; j<nc; ++j)
			{
				ri.val[j] = r.val[i*nc + j];
			}
			// rhs corrsponding to the neighbors of node i
			k = 1;
			for (l=0; l<nneigh; ++l)
			{
				p=neigh->val[nneigh*i+l];
				if ( p>=0 )
				{
					for (j=0; j<nc; ++j)
					{
						ri.val[k*nc+j] = r.val[p*nc+j];
					}
					
					++k;
				} // end if
			}
			
			ri.row = k*nc;
			//----------------------------------------------------
			//----------------------------------------------------
			// solve local problem
			e.row = k*nc;
			//fasp_blas_smat_mxv(diaginv[ti].val, ri.val, k*nc, e.val);
			fasp_smat_lu_solve(diaginv[i].val, ri.val, pivot[i].val, e.val, k*nc);
			//----------------------------------------------------
			//----------------------------------------------------
			// update solution
			// solution corresponding to node i
			for (j=0; j<nc; ++j)
			{
				u->val[i*nc + j] += e.val[j];
			}
			// solution corresponding to the neighbor of node i
			k = 1;
			for (l=0; l<nneigh; ++l)
			{
				p=neigh->val[nneigh*i+l];
				if ( p>=0 )
				{
					for (j=0; j<nc; ++j)
					{
						u->val[p*nc+j] += e.val[k*nc+j];
					}
					
					++k;
				} // end if
			}
			//----------------------------------------------------
			//----------------------------------------------------
			// update residule
			fasp_dvec_cp(b,&r); fasp_blas_dstr_aAxpy(-1.0,A,u->val,r.val);
		}
	}
	else
	{
		for (i=0; i<ngrid; ++i)
		{
			ti = order->val[i];
			//-----------------------------------------------------
			// right hand side for A_ii e_i = r_i
			// rhs corresponding to node i
			for (j=0; j<nc; ++j)
			{
				ri.val[j] = r.val[ti*nc + j];
			}
			// rhs corrsponding to the neighbors of node i
			k = 1;
			for (l=0; l<nneigh; ++l)
			{
				p=neigh->val[nneigh*ti+l];
				if ( p>=0 )
				{
					for (j=0; j<nc; ++j)
					{
						ri.val[k*nc+j] = r.val[p*nc+j];
					}
					
					++k;
				} // end if
			}
			
			ri.row = k*nc;
			//----------------------------------------------------
			//----------------------------------------------------
			// solve local problem
			e.row = k*nc;
			//fasp_blas_smat_mxv(diaginv[ti].val, ri.val, k*nc, e.val);
			fasp_smat_lu_solve(diaginv[ti].val, ri.val, pivot[ti].val, e.val, k*nc);
			//----------------------------------------------------
			//----------------------------------------------------
			// update solution
			// solution corresponding to node i
			for (j=0; j<nc; ++j)
			{
				u->val[ti*nc + j] += e.val[j];
			}
			// solution corresponding to the neighbor of node i
			k = 1;
			for (l=0; l<nneigh; ++l)
			{
				p=neigh->val[nneigh*ti+l];
				if ( p>=0 )
				{
					for (j=0; j<nc; ++j)
					{
						u->val[p*nc+j] += e.val[k*nc+j];
					}
					
					++k;
				} // end if
			}
			//----------------------------------------------------
			//----------------------------------------------------
			// update residule
			fasp_dvec_cp(b,&r); fasp_blas_dstr_aAxpy(-1.0,A,u->val,r.val);
		} 
	}// end of main loop
}


/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void blkcontr2_str(INT start_data, INT start_vecx, INT start_vecy,  
 *                                  INT nc, REAL *data, REAL *x, REAL *y)
 *
 * \brief subtract the block computation 'P*z' from 'y', where 'P' is a nc*nc  
 *        full matrix stored in 'data' from the address 'start_data', and 'z' 
 *        is a nc*1 vector stored in 'x' from the address 'start_vecx'. 
 *
 * \param start_data   starting position in data
 * \param start_vecx   starting position in x
 * \param start_vecy   starting position in y
 * \param nc           the dimension of the submatrix
 * \param *data        pointer to matrix data
 * \param *x           pointer to the REAL vector with length nc
 * \param *y           pointer to the REAL vector with length nc
 *
 * \date 04/24/2010
 */
static void blkcontr2_str(INT start_data, INT start_vecx, INT start_vecy, INT nc, 
                          REAL *data, REAL *x, REAL *y)
{
	INT i,j,k,m;
	if (start_vecy == 0)
	{
		for (i = 0; i < nc; i ++) 
		{
			k = start_data + i*nc;
			for (j = 0; j < nc; j ++) 
			{
				y[i] -= data[k+j]*x[start_vecx+j];
			}
		}
	}
	else
	{
		for (i = 0; i < nc; i ++) 
		{
			k = start_data + i*nc;
			m = start_vecy + i;
			for (j = 0; j < nc; j ++) 
			{
				y[m] -= data[k+j]*x[start_vecx+j];
			}
		}
	}
} 

/**
 * \fn static void saAxpby(REAL alpha, REAL beta, INT size, REAL *A, REAL *x, REAL *y)  
 *
 * \brief Compute y:=alpha*A*x + beta*y, here x,y are vectors, A is a size*size full matrix. 
 *  
 * \param alpha   a real number
 * \param beta    a real number
 * \param size    length of vector x and y
 * \param *A      pointer to the REAL vector which stands for a size*size full matrix 
 * \param *x      pointer to the REAL vector with length size
 * \param *y      pointer to the REAL vector with length size
 *
 * \date 04/27/2010
 */
static void saAxpby(REAL alpha, REAL beta, INT size, REAL *A, REAL *x, REAL *y)
{
	INT i,j;
	REAL tmp = 0.0;
	if (alpha == 0)
	{
		for (i = 0; i < size; i ++)
		{
			y[i] *= beta;
		}
		return;
	}
	tmp = beta / alpha;
	// y:=(beta/alpha)y
	for (i = 0; i < size; i ++)
	{
		y[i] *= tmp;
	}
	// y:=y+Ax
	for (i = 0; i < size; i ++)
	{
		for (j = 0; j < size; j ++)
		{
			y[i] += A[i*size+j]*x[j];
		}
	}
	// y:=alpha*y
	for (i = 0; i < size; i ++)
	{
		y[i] *= alpha;
	}
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
