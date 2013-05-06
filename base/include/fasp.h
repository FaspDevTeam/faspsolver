/*! \file fasp.h
 *  \brief Main header file for the FASP package
 *  
 *  This header file contains general constants and data structures used in FASP.
 * 
 *  \note Only define macros and data structures, no function decorations. 
 *
 *------------------------------------------------------
 *  Created by Chensong Zhang on 08/12/2010.
 *  Modified by Chensong Zhang on 12/13/2011.
 *  Modified by Chensong Zhang on 12/25/2011.
 *------------------------------------------------------
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "messages.h"

#ifndef __FASP_HEADER__       /*-- allow multiple inclusions --*/
#define __FASP_HEADER__       /**< indicate fasp.h has been included before */

/*---------------------------*/ 
/*---  Macros definition  ---*/
/*---------------------------*/

/**
 * \brief Flags for developer's use only
 */
// When this flag is ON, the extra printing will be enabled for debug purpose
#define DEBUG_MODE       OFF /**< output DEBUG information */
// When this flag is ON, the extra printing will be enabled for debug purpose
#define CHMEM_MODE       OFF /**< output MEMORY usage information */
// When this flag is ON, the matrix rows will be reordered as diagonal entries first
// Use with cautious!!! 
#define DIAGONAL_PREF    OFF /**< order each row such that diagonal appears first */

/**
 * \brief For external software package support
 */
#define FASP_USE_ILU     ON  /**< enable ILU or not */
#define DLMALLOC         OFF /**< use dlmalloc instead of standard malloc */
#define NEDMALLOC        OFF /**< use nedmalloc instead of standard malloc */

/**
 * \brief For Fortran compatibilty 
 */
#define SHORT            short  /**< short integer type */
#define INT              int    /**< regular integer type: int or long */
#define LONG             long   /**< long integer type */ 
#define REAL             double /**< float type */

/**
 * \brief Some global constants 
 */
#define BIGREAL          1e+20 /**< A large real number */ 
#define SMALLREAL        1e-20 /**< A small real number */ 
#define MAX_REFINE_LVL   20    /**< Maximal refinement level */
#define MAX_AMG_LVL      20    /**< Maximal AMG coarsening level */
#define STAG_RATIO       1e-4  /**< Staganation tolerance = tol*STAGRATIO */
#define MAX_STAG         20    /**< Maximal number of staganation times */
#define MAX_RESTART      20    /**< Maximal number of restarting */
#define OPENMP_HOLDS     2000  /**< Switch to sequence or openmp */

/** 
 * \brief Definition of max, min, abs
 */
#define MAX(a,b) (((a)>(b))?(a):(b)) /**< bigger one in a and b */
#define MIN(a,b) (((a)<(b))?(a):(b)) /**< smaller one in a and b */
#define ABS(a) (((a)>=0.0)?(a):-(a)) /**< absolute value of a */

/** 
 * \brief Definition of >, >=, <, <=, and isnan. 
 */
#define GT(a,b) (((a)>(b))?(TRUE):(FALSE))   /**< is a > b? */
#define GE(a,b) (((a)>=(b))?(TRUE):(FALSE))  /**< is a >= b? */
#define LS(a,b) (((a)<(b))?(TRUE):(FALSE))   /**< is a < b? */
#define LE(a,b) (((a)<=(b))?(TRUE):(FALSE))  /**< is a <= b? */
#define ISNAN(a) (((a)!=(a))?(TRUE):(FALSE)) /**< is a == NAN? */

/**
 * \brief Index starting point: C convention or Fortran convention
 */
#define ISTART 0                /**< 0 if in Natural index, 1 if data is in C index */
#define N2C(ind) ((ind)-ISTART) /**< map from Natural index 1,2,... to C index 0,1,... */
#define C2N(ind) ((ind)+ISTART) /**< map from C index 0,1,... to Natural index 1,2,... */

/*---------------------------*/ 
/*---  Global variables   ---*/
/*---------------------------*/ 

extern unsigned INT total_alloc_mem;   /**< total allocated memory */
extern unsigned INT total_alloc_count; /**< total allocation times */

/*---------------------------*/ 
/*---  Matrix and vector  ---*/
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
	
} ddenmat; /**< Dense matrix of REAL type */

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
	
} idenmat; /**< Dense matrix of INT type */

/** 
 * \struct dCSRmat
 * \brief Sparse matrix of REAL type in CSR format.
 *
 * CSR Format (IA,JA,A) in REAL
 *
 * \note The starting index of A is 0.
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
	
} dCSRmat; /**< Sparse matrix of REAL type in CSR format */

/** 
 * \struct iCSRmat
 * \brief Sparse matrix of INT type in CSR format.
 *
 * CSR Format (IA,JA,A) in integer
 *
 * \note The starting index of A is 0.
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
	
} iCSRmat; /**< Sparse matrix of INT type in CSR format */

/** 
 * \struct dCOOmat
 * \brief Sparse matrix of REAL type in COO (or IJ) format.
 *
 * Coordinate Format (I,J,A)
 *
 * \note The starting index of A is 0.
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
	
} dCOOmat; /**< Sparse matrix of REAL type in COO format */

/** 
 * \struct iCOOmat
 * \brief Sparse matrix of INT type in COO (or IJ) format.
 *
 * Coordinate Format (I,J,A)
 *
 * \note The starting index of A is 0.
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
	
} iCOOmat; /**< Sparse matrix of INT type in COO format */

/*!
 * \struct dCSRLmat
 * \brief Sparse matrix of REAL type in CSRL format.
 */
typedef struct dCSRLmat{
    
	//! number of rows	
	INT row;
	
    //! number of cols
	INT col;
	
    //! number of nonzero entries
	INT nnz;
	
    //! number of different values in i-th row, i=0:nrows-1
	INT dif;
	
    //! nz_diff[i]: the i-th different value in 'nzrow'
	INT *nz_diff;
	
    //! row index of the matrix (length-grouped): rows with same nnz are together
	INT *index;
	
    //! j in {start[i],...,start[i+1]-1} means nz_diff[i] nnz in index[j]-row
	INT *start;
	
    //! column indices of all the nonzeros
	INT *ja;
	
    //! values of all the nonzero entries
	REAL *val;
    
} dCSRLmat; /**< Sparse matrix of REAL type in CSRL format */

/** 
 * \struct dSTRmat
 * \brief Structure matrix of REAL type.
 *
 * A structure REAL matrix
 *
 * \note Every nc^2 entries of the array diag and off-diag[i] store one block:
 *		 For 2D matrix, the recommended offsets is [-1,1,-nx,nx];
 *		 For 3D matrix, the recommended offsets is [-1,1,-nx,nx,-nxy,nxy].
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
	
    //! diagonal entries (length is ngrid*(nc^2))
	REAL *diag;
	
	//! number of off-diag bands
	INT nband;	
	
    //! offsets of the off-diagals (length is nband)
	INT *offsets;
	
    //! off-diagonal entries (dimension is nband * [(ngrid-|offsets|) * nc^2])
	REAL **offdiag;
	
} dSTRmat; /**< Structured matrix of REAL type */

/** 
 * \struct dvector
 * \brief Vector with n entries of REAL type.
 */
typedef struct dvector{
	
    //! number of rows
	INT row;
    
    //! actual vector entries
	REAL *val;
	
} dvector; /**< Vector of REAL type */

/** 
 * \struct ivector
 * \brief Vector with n entries of INT type.
 */
typedef struct ivector{
	
    //! number of rows
	INT row;
    
    //! actual vector entries
	INT *val;
	
} ivector; /**< Vector of INT type */

/*---------------------------*/ 
/*--- Parameter structures --*/
/*---------------------------*/ 

/**
 * \struct ILU_param
 * \brief Parameters for ILU
 */
typedef struct {
	
	//! print leve
	SHORT print_level;
    
	//! ILU type for decomposition
	SHORT ILU_type;
	
    //! level of fill-in for ILUk
	INT ILU_lfil;
	
    //! drop tolerence for ILUt
	REAL ILU_droptol;
	
    //! add the sum of dropped elements to diagnal element in proportion relax
	REAL ILU_relax;
	
    //! permuted if permtol*|a(i,j)| > |a(i,i)|
	REAL ILU_permtol;
	
} ILU_param; /**< Parameters for ILU */	

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
	
    //! block size for BSR type only
	INT nb;
	
	//! work space size
	INT nwork;
	
    //! work space
	REAL *work;
	
} ILU_data; /**< Data for ILU */


/**
 * \struct Schwarz_param
 * \brief Parameters for Schwarz method
 *
 * Added on 05/14/2012
 */
typedef struct {
	
	//! print leve
	SHORT print_level;
	
    //! type for Schwarz method
	SHORT schwarz_type;
	
    //! maximal level for constructing the blocks
	INT schwarz_maxlvl;
    
    //! maxiaml size of blocks
    INT schwarz_mmsize;
	
} Schwarz_param; /**< Parameters for ILU */	


/** 
 * \struct Schwarz_data
 * \brief Data for Schwarz methods.
 *
 * This is needed for the schwarz solver, preconditioner/smoother.
 *
 */
typedef struct { 
    
	/* matrix information */
	
	//! pointer to the matrix
	dCSRmat A;  // note: has to be start from 1!! Change later
	
	/* blocks information */	
	
    //! number of blocks
	INT nblk;
	
    //! row index of blocks
	INT *iblock;
	
    //! column index of blocks
	INT *jblock;
	
    //! local right hand side
	REAL *rhsloc;
	
    //! LU decomposition: the U block
	REAL *au;
	
    //! LU decomposition: the L block
	REAL *al;
	
	//! Schwarz method type
	INT schwarz_type;
	
	//! working space size
	INT memt;
    
    //! mask
	INT *mask;
    
    //! maxa
	INT *maxa;
	
} Schwarz_data;

/** 
 * \struct AMG_param
 * \brief Parameters for AMG solver.
 *
 * \note This is needed for the AMG solver/preconditioner.
 */
typedef struct {
	
	//! type of AMG method
	SHORT AMG_type;
	
    //! print level for AMG
	SHORT print_level;
	
    //! max number of iterations of AMG
	INT maxit;
	
    //! stopping tolerance for AMG solver
	REAL tol;
	
	//! max number of levels of AMG	
	SHORT max_levels;
	
    //! max coarsest level dof
	INT coarse_dof;	
	
    //! type of AMG cycle
	SHORT cycle_type;
	
    //! smoother type
	SHORT smoother;
	
    //! smoother order
	SHORT smooth_order;  // 1: nature order 2: C/F order (both are symmetric)
	
    //! number of presmoothers
	SHORT presmooth_iter;
	
    //! number of postsmoothers
	SHORT postsmooth_iter;
	
    //! relaxation parameter for SOR smoother
	REAL relaxation;
    
    //! degree of the polynomial smoother
    SHORT polynomial_degree;
	
    //! switch of scaling of the coarse grid correction
	SHORT coarse_scaling;
	
    //! degree of the polynomial used by AMLI cycle
	SHORT amli_degree;
	
    //! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;
    
    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;
	
	//! coarsening type
	SHORT coarsening_type;
	
    //! interpolation type
	SHORT interpolation_type;
	
	//! strong connection threshold for coarsening
	REAL strong_threshold;
	
    //! maximal row sum parameter
	REAL max_row_sum;
	
    //! truncation threshold
	REAL truncation_threshold;
    
    //! number of levels use aggressive coarsening
    INT aggressive_level;
    
    //! numebr of paths use to determin stongly coupled C points
    INT aggressive_path;
	
	//! strong coupled threshold for aggregate
	REAL strong_coupled;
	
    //! max size of each aggregate
	INT max_aggregation;
	
    //! relaxation parameter for smoothing the tentative prolongation
	REAL tentative_smooth;
	
    //! switch for filtered matrix used for smoothing the tentative prolongation
	SHORT smooth_filter;
	
	//! number of levels use ILU smoother
	SHORT ILU_levels;
	
    //! ILU type for smoothing
	SHORT ILU_type;
	
    //! level of fill-in for ILUs and ILUk
	INT ILU_lfil;
	
    //! drop tolerence for ILUt
	REAL ILU_droptol;
	
    //! relaxiation for ILUs
	REAL ILU_relax;
	
    //! permuted if permtol*|a(i,j)| > |a(i,i)|
	REAL ILU_permtol;
    
	//! number of levels use schwarz smoother
	INT schwarz_levels;
	
    //! maximal block size
	INT schwarz_mmsize;
	
    //! maximal levels
	INT schwarz_maxlvl;
	
    //! type of schwarz method
	INT schwarz_type;
	
} AMG_param; /**< Parameters for AMG */

/** 
 * \struct AMG_data
 * \brief Data for multigrid levels.
 *
 * \note This is needed for the AMG solver/preconditioner.
 */
typedef struct { 
	
	/* Level information */
	
	//! max number of levels
	SHORT max_levels;
	
    //! number of levels in use <= max_levels
	SHORT num_levels;
	
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
    
    //! number of levels use schwarz smoother
	INT schwarz_levels;
	
    //! data of Schwarz smoother
	Schwarz_data schwarz;
	
	//! Temporary work space
	dvector w;
	
} AMG_data; /**< Data for AMG */

/** 
 * \struct precond_data
 * \brief Data passed to the preconditioners.
 *
 */
typedef struct {
	
    //! type of AMG method
	SHORT AMG_type;
	
    //! print level in AMG preconditioner
	SHORT print_level;
	
    //! max number of iterations of AMG preconditioner
	INT maxit;
	
    //! max number of AMG levels
	SHORT max_levels;
	
    //! tolerance for AMG preconditioner
	REAL tol;
	
    //! AMG cycle type
	SHORT cycle_type;
	
    //! AMG smoother type
	SHORT smoother;
	
    //! AMG smoother ordering
	SHORT smooth_order;
	
    //! number of presmoothing
	SHORT presmooth_iter;
	
    //! number of postsmoothing
	SHORT postsmooth_iter;
	
    //! coarsening type
	REAL relaxation;
    
    //! degree of the polynomial smoother
    SHORT polynomial_degree;
	
    //! switch of scaling of the coarse grid correction
	SHORT coarsening_type;
	
    //! relaxation parameter for SOR smoother
	SHORT coarse_scaling;
	
    //! degree of the polynomial used by AMLI cycle
	SHORT amli_degree;
    
    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;
	
    //! smooth factor for smoothing the tentative prolongation
	REAL tentative_smooth;
	
	//! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;

	//! AMG preconditioner data
	AMG_data *mgl_data;
	
	//! ILU preconditioner data (needed for CPR type preconditioner)
	ILU_data *LU;
	
	//! Matrix data
	dCSRmat *A;
	
	// temporary work space
    
    //! temporary dvector used to store and restore the residual
	dvector r;
    
    //! temporary work space for other usage
	REAL *w;
	
} precond_data; /**< Data for general preconditioner */

/** 
 * \struct precond_data_str
 * \brief Data passed to the preconditioner for STR format.
 *
 */
typedef struct {
	
    //! type of AMG method
	SHORT AMG_type;
	
    //! print level in AMG preconditioner
	SHORT print_level;
	
    //! max number of iterations of AMG preconditioner
	INT maxit;
	
    //! max number of AMG levels
	SHORT max_levels;
	
    //! tolerance for AMG preconditioner
	REAL tol;
	
    //! AMG cycle type
	SHORT cycle_type;
	
    //! AMG smoother type
	SHORT smoother;
	
    //! number of presmoothing
	SHORT presmooth_iter;
	
    //! number of postsmoothing
	SHORT postsmooth_iter;
	
    //! coarsening type
	SHORT coarsening_type;
	
    //! relaxation parameter for SOR smoother
	REAL relaxation;
	
    //! switch of scaling of the coarse grid correction
	SHORT coarse_scaling;
	
	//! AMG preconditioner data
	AMG_data *mgl_data;
	
	//! ILU preconditioner data (needed for CPR type preconditioner)
	ILU_data *LU;
	    
    //! whether the matrx are scaled or not
	SHORT scaled;
    
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
	
} precond_data_str; /**< Data for preconditioner of STR matrices */


/** 
 * \struct precond_diagstr
 * \brief Data passed to diagnal preconditioner for STR.
 *
 * \note This is needed for the diagnal preconditioner.
 */
typedef struct {
	
	//! number of components
	INT nc;
    
	//! diagnal elements
	dvector diag;
	
} precond_diagstr; /**< Data for diagonal preconditioner of STR matrices */

/** 
 * \struct precond_diagbsr
 * \brief Data passed to diagnal preconditioner for BSR.
 *
 * \note This is needed for the diagnal preconditioner.
 */
typedef struct {
	
	//! dimension of each sub-block
	INT nb;
    
	//! diagnal elements
	dvector diag;
	
} precond_diagbsr; /**< Data for diagonal preconditioner of BSR matrices */

/** 
 * \struct precond
 * \brief Preconditioner data and action.
 *
 * \note This is the preconditioner structure for preconditioned iterative methods.
 */ 
typedef struct {
	
	//! data for preconditioner, void pointer
	void *data;
	
	//! action for preconditioner, void function pointer
	void (*fct)(REAL *, REAL *, void *);
    
#ifdef _OPENMP 
	//! action for preconditioner, void function pointer
	void (*fct_omp)(REAL *, REAL *, void *, INT, INT);
#endif
	
} precond; /**< Data for general preconditioner passed to iterative solvers */

/**
 * \struct mxv_matfree
 * \brief Matrix-vector multiplication, replace the actual matrix.
 */
typedef struct {
	
	//! data for MxV, can be a Matrix or something else
	void *data;
    
	//! action for MxV, void function pointer
	void (*fct)(void *, REAL *, REAL *);
	
} mxv_matfree; /**< Data for general matrix passed to iterative solvers */

/**
 * \struct input_param
 * \brief Input parameters 
 *
 * Input parameters, reading from disk file
 */
typedef struct {
	
	// output flags
	SHORT print_level;   /**< print level */
	SHORT output_type;   /**< type of output stream */
	
	// problem parameters
	char workdir[256];   /**< working directory for data files */
	INT  problem_num;    /**< problem number to solve */
	
	// parameters for iterative solvers
	SHORT solver_type;   /**< type of iterative solvers */
	SHORT precond_type;  /**< type of preconditioner for iterative solvers */
	SHORT stop_type;     /**< type of stopping criteria for iterative solvers */
	REAL itsolver_tol;   /**< tolerance for iterative linear solver */
	INT itsolver_maxit;  /**< maximal number of iterations for iterative solvers */
	INT restart;         /**< restart number used in GMRES */
	
	//pamameters for ILU
	SHORT ILU_type;      /**< ILU type for decomposition*/
	INT ILU_lfil;        /**< level of fill-in */
	REAL ILU_droptol;    /**< drop tolerance */
	REAL ILU_relax;      /**< scaling factor: add the sum of dropped entries to diagnal */
	REAL ILU_permtol;    /**< permutation tolerance */
    
    // parameter for Schwarz
	INT Schwarz_mmsize;  /**< maximal block size */
	INT Schwarz_maxlvl;  /**< maximal levels */
	INT Schwarz_type;    /**< type of schwarz method */
	
	// parameters for AMG
	SHORT AMG_type;                /**< Type of AMG */
	SHORT AMG_levels;              /**< maximal number of levels */
	SHORT AMG_cycle_type;          /**< type of cycle*/
	SHORT AMG_smoother;            /**< type of smoother */
	REAL AMG_relaxation;           /**< over-relaxation parameter for SOR */
    SHORT AMG_polynomial_degree;   /**< degree of the polynomial smoother */
	SHORT AMG_presmooth_iter;      /**< number of presmoothing */
	SHORT AMG_postsmooth_iter;     /**< number of postsmoothing */
	INT AMG_coarse_dof;	           /**< minimal coarsest level dof */
	REAL AMG_tol;                  /**< tolerance for AMG if used as preconditioner */
	INT AMG_maxit;                 /**< number of iterations for AMG used as preconditioner */
	SHORT AMG_ILU_levels;          /**< how many levels use ILU smoother */	
	SHORT AMG_coarse_scaling;      /**< switch of scaling of the coarse grid correction */
	SHORT AMG_amli_degree;         /**< degree of the polynomial used by AMLI cycle */
    SHORT AMG_nl_amli_krylov_type; /**< type of krylov method used by nonlinear AMLI cycle */
    INT AMG_schwarz_levels;        /**< number of levels use schwarz smoother */
	
	// parameters for classical AMG
	SHORT AMG_coarsening_type;     /**< coarsening type */
	SHORT AMG_interpolation_type;  /**< interpolation type */
	REAL AMG_strong_threshold;     /**< strong threshold for coarsening */
	REAL AMG_truncation_threshold; /**< truncation factor for interpolation */
	REAL AMG_max_row_sum;          /**< maximal row sum */
    INT AMG_aggressive_level;      /**< number of levels use aggressive coarsening */
    INT AMG_aggressive_path;       /**< number of paths used to determine strongly coupled C-set */
	
	//  parameters for smoothed aggregation AMG
	REAL AMG_strong_coupled;       /**< strong coupled threshold for aggregate */
	INT AMG_max_aggregation;       /**< max size of each aggregate */
	REAL AMG_tentative_smooth;     /**< relaxation factor for smoothing the tentative prolongation */
	SHORT AMG_smooth_filter;       /**< use filterfor smoothing the tentative prolongation or not */
	
} input_param; /**< Input parameters */

/** 
 * \struct itsolver_param
 * \brief Parameters passed to iterative solvers.
 *
 */
typedef struct {
	
	SHORT itsolver_type; /**< solver type: see message.h */
	SHORT precond_type;  /**< preconditioner type: see message.h */
	SHORT stop_type;     /**< stopping criteria type */
	INT   maxit;         /**< max number of iterations */
	REAL  tol;           /**< convergence tolerance */
	INT   restart;       /**< number of steps for restarting: for GMRES etc */
	SHORT print_level;   /**< print level: 0--10 */	
	
} itsolver_param; /**< Parameters for iterative solvers */ 

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

	REAL (*p)[2];  /**< Coordinates of vertices */
	INT (*e)[2];   /**< Vertices of edges */
	INT (*t)[3];   /**< Vertices of triangles */
	INT (*s)[3];   /**< Edges of triangles */
	INT *pdiri;    /**< Boundary flags (0 <=> interior point) */
	INT *ediri;    /**< Boundary flags (0 <=> interior edge) */
	
	INT *pfather;  /**< Father point or edge */
	INT *efather;  /**< Father edge or triangle */
	INT *tfather;  /**< Father triangle */
	
	INT vertices;  /**< Number of grid points */
	INT edges;     /**< Number of edges */
	INT triangles; /**< Number of triangles */

} grid2d; /**< 2D grid type for plotting */

typedef grid2d *pgrid2d; /**< Grid in 2d */

typedef const grid2d *pcgrid2d; /**< Grid in 2d */

/**
 * \struct Link
 * \brief Struct for Links
 */
typedef struct
{

	//! previous node in the linklist
	INT prev;
	
    //! next node in the linklist
	INT next;
	
} Link; /**< General data structure for Links */

/** 
 * \struct linked_list
 * \brief A linked list node
 *
 * \note This definition is adapted from hypre 2.0.
 */
typedef struct linked_list
{
    
    //! data
	INT data;
    
    //! starting of the list
	INT head;
	
    //! ending of the list
	INT tail;

    //! next node
	struct linked_list *next_node;
    
	//! previous node
	struct linked_list *prev_node;
        
} ListElement; /**< Linked element in list */

/** 
 * List of links
 */
typedef ListElement *LinkList; /**< linked list */

/*
 * OpenMP definitions and decrorations
 */
#define FASP_GSRB 1  /*mg level 0  use RedBlack Gausss Seidel Smoothing */

#if FASP_GSRB
extern INT  nx_rb ;  /**< Red Black Gs Smoother Nx */
extern INT  ny_rb ;  /**< Red Black Gs Smoother Ny */
extern INT  nz_rb ;  /**< Red Black Gs Smoother Nz */
extern INT *IMAP;    /**< Red Black Gs Smoother imap */
extern INT  MAXIMAP; /**< Red Black Gs Smoother max dofs of reservoir */
#endif

#ifdef _OPENMP 

#include "omp.h"

extern INT THDs_AMG_GS;  /**< AMG GS smoothing threads  */
extern INT THDs_CPR_lGS; /**< Reservoir GS smoothing threads */
extern INT THDs_CPR_gGS; /**< Global matrix GS smoothing threads  */

extern REAL total_linear_time; /**< total linear times */
extern REAL total_setup_time;  /**< ??? */
extern REAL total_start_time;  /**< ??? */
extern INT  total_iter;        /**< ??? */
extern INT  fasp_called_times; /**< ??? */

#endif /* end if for FASP_USE_OPENMP */

#endif /* end if for __FASP_HEADER__ */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
