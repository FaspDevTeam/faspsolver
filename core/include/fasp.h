/*
 *  fasp.h
 *
 *------------------------------------------------------
 *  Created by Chensong Zhang on 08/12/2010.
 *  Modified by Chensong Zhang on 12/13/2011.
 *  Modified by Chensong Zhang on 12/25/2011.
 *------------------------------------------------------
 *
 */

/*! \file fasp.h
 *  \brief Main header file for the FASP package
 *  
 *  Note: Only define macros and data structures, no function decorations. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "messages.h"

#ifndef __FASP_HEADER__		/*-- allow multiple inclusions --*/
#define __FASP_HEADER__

/*---------------------------*/ 
/*---  Macros definition  ---*/
/*---------------------------*/

/**
 * \brief For internal use only
 */
#define DEBUG_MODE       OFF  /**< output DEBUG information */
#define CHMEM_MODE       OFF  /**< output MEMORY usage information */
#define DIAGONAL_PREF    OFF   /**< order each row such that diagonal appears first */

/**
 * \brief For external software support
 */
#define FASP_USE_ILU     OFF  /**< enable ILU or not */
#define FASP_USE_OPENMP  ON  /**< enable OpenMP support or not */

/**
 * \brief For external memory management 
 */
#define DLMALLOC         OFF  /**< use dlmalloc instead of standard malloc */
#define NEDMALLOC        OFF  /**< use nedmalloc instead of standard malloc */

/**
 * \brief For Fortran compatibilty 
 */
#define INT    int
#define REAL   double

/**
 * \brief Some global constants 
 */
#define BIGREAL     1e+15 /**< A large real number */ 
#define SMALLREAL   1e-15 /**< A small real number */ 
#define MAX_REFINE_LVL 20 /**< Maximal refinement level */
#define MAX_AMG_LVL    20 /**< Maximal AMG coarsening level */

/** 
 * \brief Definition of max, min, abs
 */
#define MAX(a,b) (((a)>(b))?(a):(b)) /**< bigger one in a and b */
#define MIN(a,b) (((a)<(b))?(a):(b)) /**< smaller one in a and b */
#define ABS(a) (((a)>=0.0)?(a):-(a)) /**< absolute value of a */

/** 
 * \brief Definition of max, min, abs
 */
#define GT(a,b) (((a)>(b))?(TRUE):(FALSE))  /**< is a > b? */
#define GE(a,b) (((a)>=(b))?(TRUE):(FALSE)) /**< is a >= b? */
#define LS(a,b) (((a)<(b))?(TRUE):(FALSE))  /**< is a < b? */
#define LE(a,b) (((a)<=(b))?(TRUE):(FALSE)) /**< is a <= b? */

/**
 * \brief Index starting point: C convention or Fortran convention
 */

#define ISTART 0 /**< 0 if in Natural index, 1 if data is in C index */
#define N2C(ind) ((ind)-ISTART) /**< map from Natural index 1,2,... to C index 0,1,... */
#define C2N(ind) ((ind)+ISTART) /**< map from C index 0,1,... to Natural index 1,2,... */

/*---------------------------*/ 
/*---  Global variables   ---*/
/*---------------------------*/ 

extern unsigned INT total_alloc_mem;   /**< total allocated memory */
extern unsigned INT total_alloc_count; /**< total allocation times */

/*---------------------------*/ 
/*---   Data structures   ---*/
/*---------------------------*/ 

/** 
 * \struct ddenmat
 * \brief Dense matrix of REAL type.
 *
 * A dense REAL matrix
 */ 
typedef struct ddenmat{
	
	//! number of rows
	INT row;	  
	//! number of columns
	INT col;	
	//! actual matrix entries
	REAL **val;
	
} ddenmat;

/** 
 * \struct idenmat
 * \brief Dense matrix of INT type.
 *
 * A dense INT matrix
 */ 
typedef struct idenmat{
	
	//! number of rows
	INT row;	  
	//! number of columns
	INT col;	
	//! actual matrix entries 
	INT **val;
	
} idenmat;

/** 
 * \struct dCSRmat
 * \brief Sparse matrix of REAL type in CSR format.
 *
 * CSR Format (IA,JA,A) in REAL
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct dCSRmat{
	
	//! row number of matrix A, m
	INT row;   
	//! column of matrix A, n
	INT col;   
	//! number of nonzero entries
	INT nnz;
	//! integer array of row pointers, the size is m+1
	INT *IA;   
	//! integer array of column indexes, the size is nnz
	INT *JA;    
	//! nonzero entries of A
	REAL *val;
	
} dCSRmat;

/** 
 * \struct iCSRmat
 * \brief Sparse matrix of INT type in CSR format.
 *
 * CSR Format (IA,JA,A) in integer
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct iCSRmat{
	
	//! row number of matrix A, m
	INT row;   
	//! column of matrix A, n
	INT col;   
	//! number of nonzero entries
	INT nnz;
	//! integer array of row pointers, the size is m+1
	INT *IA;   
	//! integer array of column indexes, the size is nnz
	INT *JA;    
	//! nonzero entries of A
	INT *val;
	
} iCSRmat;

/** 
 * \struct dCOOmat
 * \brief Sparse matrix of REAL type in COO (or IJ) format.
 *
 * Coordinate Format (I,J,A)
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.   
 */
typedef struct dCOOmat{
	
	//! row number of matrix A, m
	INT row;   
	//! column of matrix A, n
	INT col;   
	//! number of nonzero entries
	INT nnz;
	//! integer array of row indices, the size is nnz
	INT *I;   
	//! integer array of column indices, the size is nnz
	INT *J;    
	//! nonzero entries of A
	REAL *val;
	
} dCOOmat;

/** 
 * \struct iCOOmat
 * \brief Sparse matrix of INT type in COO (or IJ) format.
 *
 * Coordinate Format (I,J,A)
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.   
 */
typedef struct iCOOmat{
	
	//! row number of matrix A, m
	INT row;   
	//! column of matrix A, n
	INT col;   
	//! number of nonzero entries
	INT nnz;
	//! integer array of row indices, the size is nnz
	INT *I;   
	//! integer array of column indices, the size is nnz
	INT *J;    
	//! nonzero entries of A
	INT *val;
	
} iCOOmat;

/*!
 * \struct dCSRLmat
 * \brief Sparse matrix of REAL type in CSRL format.
 */
typedef struct
{
	//! number of rows	
	INT num_rows;
	//! number of cols
	INT num_cols;
	//! number of nonzero entries
	INT num_nonzeros;
	//! number of different values in nzrow[i], i=0:nrows-1, i.e. number of nonzeros in the i-th row
	INT dif;
	//! nzdifnum[i]: the i-th different value in 'nzrow'
	INT *nzdifnum;
	//! row index of the matrix in length-grouped manner, i.e., the rows with same nnz are together
	INT *rowindex;
	//! if j\in V = {rowstart[i],...,rowstart[i+1]-1}, then there are nzdifnum[i] nnz in the rowindex[j]-th row	
	INT *rowstart;
	//! column indices of all the nonzeros
	INT *ja;
	//! values of all the nonzero entries
	REAL *data;
    
} dCSRLmat;

/** 
 * \struct dvector
 * \brief Vector with n entries of REAL type.
 */
typedef struct dvector{
	
    //! number of rows
	INT row;
    //! actual vector entries
	REAL *val;
	
} dvector;

/** 
 * \struct ivector
 * \brief Vector with n entries of INT type.
 */
typedef struct ivector{
	
    //! number of rows
	INT row;
    //! actual vector entries
	INT *val;
	
} ivector;

/** 
 * \struct dSTRmat
 * \brief Structure matrix of REAL type.
 *
 * A structure REAL matrix
 *
 * Note: Every nc^2 entries of the array diag and off-diag[i] store one block.
 *		For 2D matrix, the recommended offsets is [-1,1,-nx,nx];
 *		For 3D matrix, the recommended offsets is [-1,1,-nx,nx,-nxy,nxy]; 
 */ 

typedef struct dSTRmat{
	
	//! number of grids in x direction
	INT nx;	  
	//! number of grids in y direction
	INT ny;	
	//! number of grids in z direction
	INT nz;
	//! number of grids on x-y plane
	INT nxy;
	//! size of each block (number of components)
	INT nc;
	
	//! number of grids
	INT ngrid;
	//! Diagonal entries (length is ngrid*(nc^2))
	REAL *diag;
	
	//! number of off-diag bands
	INT nband;	
	//! offsets of the off-diagals (length is nband)
	INT *offsets;
	//! Off-diagonal entries (dimension is nband * [(ngrid-|offsets|) * nc^2])
	REAL **offdiag;
	
} dSTRmat;

/*---------------------------*/ 
/*--- Parameter structures --*/
/*---------------------------*/ 

/**
 * \struct ILU_param
 * \brief Parameters for ILU
 */
typedef struct {
	
	//! prINT leve
	INT print_level;
	//! ILU type for decomposition
	INT ILU_type;
	//! level of fill-in for ILUk
	INT ILU_lfil;
	//! drop tolerence for ILUt
	REAL ILU_droptol;
	//! add the sum of dropped elements to diagnal element in proportion relax
	REAL ILU_relax;
	//! permuted if permtol*|a(i,j)| > |a(i,i)|
	REAL ILU_permtol;
	
} ILU_param;	

/** 
 * \struct ILU_data
 * \brief Data for ILU setup
 */
typedef struct {
	
	//! row number of matrix LU, m
	INT row;   
	//! column of matrix LU, n
	INT col;   
	//! number of nonzero entries
	INT nzlu;
	//! integer array of row pointers and column indexes, the size is nzlu 
	INT *ijlu;   
	//! nonzero entries of LU
	REAL *luval;
	//! block size for bsr
	INT nb;
	
	//! work space size
	INT nwork;
	//! work space
	REAL *work;
	
} ILU_data;

/** 
 * \struct AMG_param
 * \brief Parameters for AMG solver.
 *
 * This is needed for the AMG solver/preconditioner.
 *
 */
typedef struct {
	
	//! type of AMG method
	INT AMG_type;
	//! print level for AMG
	INT print_level;
	//! max number of iterations of AMG
	INT max_iter;
	//! stopping tolerance for AMG solver
	REAL tol;
	
	//! max number of levels of AMG	
	INT max_levels;
	//! max coarsest level dof
	INT coarse_dof;	
	//! type of AMG cycle
	INT cycle_type;
	//! smoother type
	INT smoother;
	//! Smoother order type
	INT smooth_order;  // 1: nature order 2: C/F order (both are symmetric)
	//! number of presmoothers
	INT presmooth_iter;
	//! number of postsmoothers
	INT postsmooth_iter;
	//! relaxation parameter for SOR smoother
	REAL relaxation;
	//! switch of scaling of the coarse grid correction
	INT coarse_scaling;
	//! degree of the polynomial used by AMLI cycle
	INT amli_degree;
	//! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;
	
	//! coarsening type
	INT coarsening_type;
	//! interpolation type
	INT interpolation_type;
	
	//! strong connection threshold for coarsening
	REAL strong_threshold;
	//! maximal row sum parameter
	REAL max_row_sum;
	//! truncation threshold
	REAL truncation_threshold;
	
	//! strong coupled threshold for aggregate
	REAL strong_coupled;
	//! max size of each aggregate
	INT max_aggregation;
	//! relaxation parameter for smoothing the tentative prolongation
	REAL tentative_smooth;
	//! switch for filtered matrix used for smoothing the tentative prolongation
	INT smooth_filter;
	
	//! number of levels use ILU smoother
	INT ILU_levels;
	//! ILU type for smoothing
	INT ILU_type;
	//! level of fill-in for ILUs and ILUk
	INT ILU_lfil;
	//! drop tolerence for ILUt
	REAL ILU_droptol;
	//! relaxiation for ILUs
	REAL ILU_relax;
	//! permuted if permtol*|a(i,j)| > |a(i,i)|
	REAL ILU_permtol;
	
} AMG_param;

/** 
 * \struct AMG_data
 * \brief Data for multigrid levels.
 *
 * This is needed for the AMG solver/preconditioner.
 *
 */
typedef struct { 
	
	/* Level information */
	
	//! max number of levels
	INT max_levels;
	//! number of levels in use <= max_levels
	INT num_levels;
	
	/* Problem information */	
	
	//! pointer to the matrix at level level_num
	dCSRmat A;
	//! restriction operator at level level_num
	dCSRmat R;
	//! prolongation operator at level level_num
	dCSRmat P;
	//! pointer to the right-hand side at level level_num
	dvector b;
	//! pointer to the iterative solution at level level_num
	dvector x;
	
	/* Extra information */
	
	//! pointer to the CF marker at level level_num
	ivector cfmark; 	
	//! number of levels use ILU smoother
	INT ILU_levels;
	//! ILU matrix for ILU smoother 	
	ILU_data LU;
	//! dimension of the near kernel for SAMG
	INT near_kernel_dim;
	//! basis of near kernel space for SAMG
	REAL **near_kernel_basis;
	// Smoother order information
	
	//! Temporary work space
	dvector w;
	
} AMG_data;

/** 
 * \struct precond_data
 * \brief Data passed to the preconditioners.
 *
 */
typedef struct {
	
    //! type of AMG method
	INT AMG_type;
	//! print level in AMG preconditioner
	INT print_level;
	//! max number of iterations of AMG preconditioner
	INT max_iter;
	//! max number of AMG levels
	INT max_levels;
	//! tolerance for AMG preconditioner
	REAL tol;
	//! AMG cycle type
	INT cycle_type;
	//! AMG smoother type
	INT smoother;
	//! AMG smoother ordering
	INT smooth_order;
	//! number of presmoothing
	INT presmooth_iter;
	//! number of postsmoothing
	INT postsmooth_iter;
	//! coarsening type
	INT coarsening_type;
	//! relaxation parameter for SOR smoother
	REAL relaxation;
	//! switch of scaling of the coarse grid correction
	INT coarse_scaling;
	//! degree of the polynomial used by AMLI cycle
	INT amli_degree;
	//! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;
	//! smooth factor for smoothing the tentative prolongation
	REAL tentative_smooth;
	
	//! AMG preconditioner data
	AMG_data *mgl_data;
	
	//! ILU preconditioner data (needed for CPR type preconditioner)
	ILU_data *LU;
	
	//! Matrix data
	dCSRmat *A;
	
	// temporary work space
	dvector r; /*< temporary dvector used to store and restore the residual */
	REAL *w; /*<  temporary work space for other usage */
	
} precond_data;

/** 
 * \struct precond_data_str
 * \brief Data passed to the preconditioner for STR format.
 *
 */
typedef struct {
	
    //! type of AMG method
	INT AMG_type;
	//! print level in AMG preconditioner
	INT print_level;
	//! max number of iterations of AMG preconditioner
	INT max_iter;
	//! max number of AMG levels
	INT max_levels;
	//! tolerance for AMG preconditioner
	REAL tol;
	//! AMG cycle type
	INT cycle_type;
	//! AMG smoother type
	INT smoother;
	//! number of presmoothing
	INT presmooth_iter;
	//! number of postsmoothing
	INT postsmooth_iter;
	//! coarsening type
	INT coarsening_type;
	//! relaxation parameter for SOR smoother
	REAL relaxation;
	//! switch of scaling of the coarse grid correction
	INT coarse_scaling;
	
	//! AMG preconditioner data
	AMG_data *mgl_data;
	
	//! ILU preconditioner data (needed for CPR type preconditioner)
	ILU_data *LU;
	
	// Matrix data
    
    //! whether the matrx are scaled or not
	INT scaled;
    //! the orginal CSR matrix
	dCSRmat *A;
    //! stor the whole reservoir block in STR format
	dSTRmat *A_str;
    //! store Saturation block in STR format    
	dSTRmat *SS_str;
	
	// data for GS/block GS smoothers (STR format)
    
    //! the inverse of the diagonals for GS/block GS smoother (whole reservoir matrix)
	dvector *diaginv;
    //! the pivot for the GS/block GS smoother (whole reservoir matrix)
	ivector *pivot;
    //! the inverse of the diagonals for GS/block GS smoother (saturation block)
	dvector *diaginvS;
    //! the pivot for the GS/block GS smoother (saturation block)
	ivector *pivotS;
    //! order for smoothing
	ivector *order;
    //! arrary to store neighbor information
	ivector *neigh;
	
	// temporary work space
    
    //! temporary dvector used to store and restore the residual
	dvector r; 
    //! temporary work space for other usage
	REAL *w; 
	
} precond_data_str;


/** 
 * \struct precond_diagstr
 * \brief Data passed to diagnal preconditioner for STR.
 *
 * This is needed for the diagnal preconditioner.
 */

typedef struct {
	
	//! number of components
	INT nc;
	
	//! diagnal elements
	dvector diag;
	
} precond_diagstr;

/** 
 * \struct precond_diagbsr
 * \brief Data passed to diagnal preconditioner for BSR.
 *
 * This is needed for the diagnal preconditioner.
 */

typedef struct {
	
	//! dimension of each sub-block
	INT nb;
	
	//! diagnal elements
	dvector diag;
	
} precond_diagbsr;

/** 
 * \struct precond
 * \brief Preconditioner data and action.
 *
 * This is the preconditioner structure for preconditioned Krylov methods.
 */ 
typedef struct {
	
	//! data for preconditioner, void pointer
	void *data; 	
	//! action for preconditioner, void function pointer
	void (*fct)(REAL *, REAL *, void *);
#if FASP_USE_OPENMP
	//! action for preconditioner, void function pointer
	void (*fct_omp)(REAL *, REAL *, void *, INT, INT);
#endif
	
} precond;

/**
 * \struct input_param
 * \brief Input parameters 
 *
 * Input parameters, reading from disk file
 */
typedef struct {
	
	// output flags
	INT print_level; /**< print level */
	INT output_type; /**< type of output stream */
	
	// problem parameters
	char workdir[256]; /**< working directory for data files */
	INT  problem_num; /**< problem number to solve */
	
	// parameters for iterative solvers
	INT itsolver_type; /**< type of iterative solvers */
	INT precond_type; /**< type of preconditioner for iterative solvers */
	INT stop_type; /**< type of stopping criteria for iterative solvers */
	REAL itsolver_tol; /**< tolerance for iterative linear solver */
	INT itsolver_maxit; /**< maximal number of iterations for iterative solvers */
	INT restart; /**< restart number used in GMRES */
	
	//pamameters for ILU
	INT ILU_type; /**< ILU type for decomposition*/
	INT ILU_lfil; /**< level of fill-in */
	REAL ILU_droptol; /**< drop tolerence */
	REAL ILU_relax; /**< add the sum of dropped elements to diagnal element in proportion relax */
	REAL ILU_permtol; /**< permutation tol */
	
	// parameters for AMG
	INT AMG_type; /**< Type of AMG */
	INT AMG_levels; /**< maximal number of levels */
	INT AMG_cycle_type; /**< type of cycle*/
	INT AMG_smoother; /**< type of smoother */
	REAL AMG_relaxation; /**< over-relaxation parameter for SOR */
	INT AMG_presmooth_iter; /**< number of presmoothing */
	INT AMG_postsmooth_iter; /**< number of postsmoothing */
	INT AMG_coarse_dof;	/**< minimal coarsest level dof */
	REAL AMG_tol; /**< tolerance for AMG if used as preconditioner */
	INT AMG_maxit; /**< max number of iterations for AMG if used as preconditioner */
	INT AMG_ILU_levels; /**< how many levels use ILU smoother */	
	INT AMG_coarse_scaling; /**< switch of scaling of the coarse grid correction */
	INT AMG_amli_degree; /**< degree of the polynomial used by AMLI cycle */
	
	// parameters for classical AMG
	INT AMG_coarsening_type; /**< coarsening type */
	INT AMG_interpolation_type; /**< interpolation type */
	REAL AMG_strong_threshold; /**< strong threshold for coarsening */
	REAL AMG_truncation_threshold; /**< truncation factor for interpolation */
	REAL AMG_max_row_sum; /**< maximal row sum */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled; /**< strong coupled threshold for aggregate */
	INT AMG_max_aggregation; /**< max size of each aggregate */
	REAL AMG_tentative_smooth; /**< relaxation parameter for smoothing the tentative prolongation */
	INT AMG_smooth_filter; /**< switch for filtered matrix used for smoothing the tentative prolongation */
	
} input_param;

/** 
 * \struct itsolver_param
 * \brief Data passed to iterative solvers.
 *
 */
typedef struct {
	
	INT print_level; /**< print level */	
	INT itsolver_type; /**< solver type */
	INT precond_type; /**< preconditioner type */
	INT stop_type; /**< stopping criteria type */
	INT maxit; /**< solve params, max number of iterations */
	INT restart; /**< number of steps for restarting the solver */
	REAL tol;	/**< solve params, tolerance for solver */
	
} itsolver_param; 

/**
 * \struct grid2d
 * \brief 2d grid data structure
 *
 * \note The grid2d structure is simply a list of triangles, 
 *       edges and vertices. 
 *		 edge i has 2 vertices e[i], 
 *		 triangle i has 3 edges s[i], 3 vertices t[i]
 *		 vertex i has two coordinates p[i]
 */
typedef struct grid2d {
	REAL (*p)[2];	/**< Coordinates of vertices */
	INT (*e)[2]; /**< Vertices of edges */
	INT (*t)[3]; /**< Vertices of triangles */
	INT (*s)[3]; /**< Edges of triangles */
	INT *pdiri; /**< Boundary flags (0 <=> interior point) */
	INT *ediri; /**< Boundary flags (0 <=> interior edge) */
	
	INT *pfather; /**< Father point or edge */
	INT *efather; /**< Father edge or triangle */
	INT *tfather; /**< Father triangle */
	
	INT vertices; /**< Number of grid points */
	INT edges; /**< Number of edges */
	INT triangles; /**< Number of triangles */
} grid2d;
typedef grid2d *pgrid2d;
typedef const grid2d *pcgrid2d;

/**
 * \struct Link
 * \brief Struct for linked lists.
 *
 */
typedef struct
{
	//! previous node in the linklist
	INT prev;
	//! next node in the linklist
	INT next;
	
} Link;

/** 
 * \struct linked_list
 * \brief Standard linked list element, adapted from hypre2.0
 */
typedef struct linked_list
{
	
    //! data
	INT data;
	//! next element
	struct linked_list *next_elt;
	//! previous element
	struct linked_list *prev_elt;
	//! starting of the list
	INT head;
	//! ending of the list
	INT tail;
	
} ListElement;

/** 
 * List of links
 */
typedef ListElement *LinkList;  

/*
 * OpenMP definitions and decrorations
 */
#if FASP_USE_OPENMP

#include "omp.h"

extern INT THDs_AMG_GS;
extern INT THDs_CPR_lGS;
extern INT THDs_CPR_gGS;

extern double total_linear_time; /**< total linear times */
extern double total_setup_time;
extern double total_start_time;
extern int total_iter;
extern int fasp_called_times;
extern int  nx_rb ;  // Red Black Gs Smoother 
extern int  ny_rb ;  // Red Black Gs Smoother
extern int  nz_rb ;  // Red Black Gs Smoother
extern  int *IMAP;   // Red Black Gs Smoother
    //! tmp map for the level 0 grid, geometry to algebre dofs.
extern  int MAXIMAP;  // Red Black Gs Smoother  max dofs of reservoir

#define FASP_GET_START_END(procid,nprocs,n,start,end) \
if((procid)<(n)%(nprocs)) \
{ \
(end)=(n)/(nprocs)+1; \
(start)=(end)*(procid); \
} \
else \
{ \
(end)=(n)/(nprocs); \
(start)=(end)*(procid)+(n)%(nprocs); \
} \
(end)=(end)+(start);

#endif

#endif /* end if for __FASP_HEADER__ */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
