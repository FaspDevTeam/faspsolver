/*
 *  fasp_block.h
 *
 *------------------------------------------------------
 *  Created by Chensong Zhang on 05/21/2010.
 *  Modified by Xiaozhe Hu on 05/28/2010: added precond_block_reservoir_data.
 *  Modified by Xiaozhe Hu on 06/15/2010: modify precond_block_reservoir_data.
 *  Modified by Chensong Zhang on 10/11/2010: added BSR data.
 *------------------------------------------------------
 *
 */

/*! \file fasp_block.h
 *  \brief Main header file for the FASP_BLK package
 *  
 *  Note: Only define macros and data structures, no function decorations. 
 */

#include "fasp.h"

#ifndef __FASPBLK_HEADER__		/*-- allow multiple inclusions --*/
#define __FASPBLK_HEADER__

/*---------------------------*/ 
/*---   Data structures   ---*/
/*---------------------------*/ 

/** 
 * \struct dBSRmat
 * \brief Block sparse row storage matrix of REAL type.
 *
 * Note: This data structure is adapted from the Intel MKL library. 
 * Refer to http://software.intel.com/sites/products/documentation/hpc/mkl/lin/index.htm
 *
 */
typedef struct dBSRmat{
	
    // Note: Some of the following entries are captialized to stress they are for blocks!

	//! number of rows of sub-blocks in matrix A, M
	INT ROW;
	//! number of cols of sub-blocks in matrix A, N
	INT COL;
	//! number of nonzero sub-blocks in matrix A, NNZ
	INT NNZ;
	
    //! dimension of each sub-block
	INT nb; // for the moment, allow nb*nb full block 
	//! storage manner for each sub-block           
	INT storage_manner; // 0: row-major order, 1: column-major order
	
	//! A real array that contains the elements of the non-zero blocks of 
	//! a sparse matrix. The elements are stored block-by-block in row major 
	//! order. A non-zero block is the block that contains at least one non-zero 
	//! element. All elements of non-zero blocks are stored, even if some of 
	//! them is equal to zero. Within each nonzero block elements are stored 
	//! in row-major order and the size is (NNZ*nb*nb). 
	REAL *val;
	
	//! integer array of row pointers, the size is ROW+1
	INT *IA;
	
	//! Element i of the integer array columns is the number of the column in the
	//! block matrix that contains the i-th non-zero block. The size is NNZ.
	INT *JA;
	
} dBSRmat; /**< Matrix of REAL type in BSR format */

/** 
 * \struct block_dCSRmat
 * \brief Block REAL CSR matrix structure.
 *
 * Block CSR Format in REAL
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct block_dCSRmat{
	
	//! row number of blocks in A, m
	INT brow;   
	//! column number of blocks A, n
	INT bcol;   
	//! blocks of dCSRmat, point to blocks[brow][bcol]
	dCSRmat **blocks;
	
} block_dCSRmat; /**< Matrix of REAL type in Block BSR format */

/** 
 * \struct block_iCSRmat
 * \brief Block integer CSR matrix structure.
 *
 * Block CSR Format in integer
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct block_iCSRmat{
	
	//! row number of blocks in A, m
	INT brow;   
	//! column number of blocks A, n
	INT bcol;   
	//! blocks of iCSRmat, point to blocks[brow][bcol]
	iCSRmat **blocks;
	
} block_iCSRmat; /**< Matrix of INT type in Block BSR format */

/** 
 * \struct block_dvector
 * \brief Block REAL vector structure.
 *
 * Block Vector Format in REAL
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct block_dvector{
	
	//! row number of blocks in A, m
	INT brow;   
	//! blocks of dvector, point to blocks[brow]
	dvector **blocks;
	
} block_dvector; /**< Vector of REAL type in Block format */

/** 
 * \struct block_ivector
 * \brief Block integer vector structure.
 *
 * Block Vector Format in integer
 *
 * Note: The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct block_ivector{
	
	//! row number of blocks in A, m
	INT brow;   
	//! blocks of dvector, point to blocks[brow]
	ivector **blocks;
	
} block_ivector; /**< Vector of INT type in Block format */

/** 
 * \struct block_Reservoir
 * \brief Block REAL matrix structure for reservoir simulation.
 *
 */
typedef struct block_Reservoir{
	
	//! reservoir-reservoir block
	dSTRmat ResRes;
	//! reservoir-well block
	dCSRmat ResWel;	
	//! well-reservoir block
	dCSRmat WelRes;
	//! well-well block
	dCSRmat WelWel;
	
} block_Reservoir; /**< Special bloack matrix for Reservoir Simulation */

/** 
 * \struct block_BSR
 * \brief Block REAL matrix structure for reservoir simulation.
 *
 */
typedef struct block_BSR{
	
	//! reservoir-reservoir block
	dBSRmat ResRes;
	//! reservoir-well block
	dCSRmat ResWel;	
	//! well-reservoir block
	dCSRmat WelRes;
	//! well-well block
	dCSRmat WelWel;
	
} block_BSR; /**< Block of BSR matrices of REAL type */

/*---------------------------*/ 
/*--- Parameter structures --*/
/*---------------------------*/ 

/** 
 * \struct AMG_data_bsr
 * \brief Data for multigrid levels. (BSR format)
 *
 * This is needed for the AMG solver/preconditioner for BSR format
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
	dBSRmat A;
	//! restriction operator at level level_num
	dBSRmat R;
	//! prolongation operator at level level_num
	dBSRmat P;
	//! pointer to the right-hand side at level level_num
	dvector b;
	//! pointer to the iterative solution at level level_num
	dvector x;
	//! pointer to the matrix at level level_num (csr format)
	dCSRmat Ac;
	
	/* Extra information */
	
	//! pointer to the pressure block (only for reservoir simulation)
	dCSRmat PP;
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
	
} AMG_data_bsr; /**< AMG data for BSR matrices */

/** 
 * \struct precond_data_bsr
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
	INT max_levels;
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
	SHORT coarsening_type;
	//! relaxation parameter for SOR smoother
	REAL relaxation;
	//! switch of scaling of the coarse grid correction
	SHORT coarse_scaling;
	//! degree of the polynomial used by AMLI cycle
	SHORT amli_degree;
	//! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;
	//! smooth factor for smoothing the tentative prolongation
	REAL tentative_smooth;
    //! type of krylov method used by Nonlinear AMLI cycle
    SHORT nl_amli_krylov_type;
	
	//! AMG preconditioner data
	AMG_data_bsr *mgl_data;
	
	//! ILU preconditioner data (needed for CPR type preconditioner)
	ILU_data *LU;
	
	//! Matrix data
	dBSRmat *A;
	
	//! temporary work space
	dvector r; /**< temporary dvector used to store and restore the residual */
	REAL *w; /**<  temporary work space for other usage */
	
} precond_data_bsr; /**< Preconditioner data for BSR matrices */


/**
 * \brief Data passed to the preconditioner for preconditioning reservoir simulation problems
 *
 * This is needed for the Black oil model with wells
 */
typedef struct precond_block_reservoir_data {
	
	//! basic data: original matrix
	block_Reservoir *A; /**< problem data in block_Reservoir format */
	block_dCSRmat *Abcsr; /**< problem data in block_dCSRmat format */
	dCSRmat *Acsr; /**< problem data in CSR format */
	
	//! data of ILU decomposition
	INT ILU_lfil; /**< level of fill-in for structured ILU(k) */
	dSTRmat *LU; /**< LU matrix for ResRes block */
	ILU_data *LUcsr; /**< LU matrix for ResRes block */
	
	//! data of AMG
	AMG_data *mgl_data; /**< AMG data for presure-presure block */
	
	//! parameters for AMG solver
	//! prINT level in AMG preconditioner
	SHORT print_level;
	//! max number of iterations of AMG preconditioner
	INT maxit_AMG;
	//! max number of AMG levels
	SHORT max_levels;
	//! tolerance for AMG preconditioner
	REAL amg_tol;
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
	//! switch of scaling of coarse grid correction
	SHORT coarse_scaling;
	
	//! parameters for krylov method used for blocks
	INT  maxit; /**< max number of iterations */
	INT  restart; /**< number of iterations for restart */
	REAL tol; /**< tolerance */ 
	
	//! inverse of the schur complement (-I - Awr*Arr^{-1}*Arw)^{-1}, Arr may be replaced by LU
	REAL *invS;
	
	//! Diag(PS) * inv(Diag(SS)) 
	dvector *DPSinvDSS;
	
	//! Data for FASP
	SHORT scaled; /**< whether the matirx is scaled */
	ivector *perf_idx; /**< variable index for perf */
	
	dSTRmat *RR;  /**< Diagonal scaled reservoir block */
	dCSRmat *WW; /**< Argumented well block */
	dCSRmat *PP; /**< pressure block after diagonal scaling */
	dSTRmat *SS; /**< saturation block after diaogonal scaling */
	
	precond_diagstr *diag; /**< the diagonal inverse for diagonal scaling */
	dvector *diaginv; /**< the inverse of the diagonals for GS/block GS smoother (whole reservoir matrix) */
	ivector *pivot; /**< the pivot for the GS/block GS smoother (whole reservoir matrix) */
	dvector *diaginvS; /**< the inverse of the diagonals for GS/block GS smoother (saturation block) */
	ivector *pivotS; /**< the pivot for the GS/block GS smoother (saturation block) */
	ivector *order; /**< order for smoothing */
	
	//! temporary work space
	dvector r; /**< temporary dvector used to store and restore the residual */
	REAL *w; /**<  temporary work space for other usage */
	
} precond_block_reservoir_data; /**< Precond data for Reservoir Simulation */

/** 
 * \brief Data passed to the preconditioner for block diagonal preconditioning.
 *
 * This is needed for the diagnoal block preconditioner.
 */
typedef struct {
	
	dCSRmat  *A; /**< problem data, the sparse matrix */
	dvector  *r; /**< problem data, the right-hand side vector */
	
	dCSRmat **Ablock; /**< problem data, the blocks */
	ivector **row_idx; /**< problem data, row indices */
	ivector **col_idx; /**< problem data, col indices */
	
	AMG_param *amgparam; /**< parameters for AMG */
	dCSRmat  **Aarray; /**< data generated in the setup phase */	
	
} precond_block_data; /**< Precond data for block matrices */

/**
 * \brief Data passed to the preconditioner for preconditioning reservoir simulation problems
 *
 * This is needed for the Black oil model with wells
 */
typedef struct precond_FASP_blkoil_data{
	
	//-------------------------------------------------------------------------------
	//! Part 1: Basic data
	//-------------------------------------------------------------------------------
	block_BSR *A; /**< whole jacobian system in block_BSRmat */
	
	//-------------------------------------------------------------------------------
	//! Part 2: Data for CPR-like preconditioner for reservoir block
	//-------------------------------------------------------------------------------
	//! diagonal scaling for reservoir block
	SHORT scaled; /**< scaled = 1 means the the following RR block is diagonal scaled */
	dvector *diaginv_noscale; /**< inverse of diagonal blocks for diagonal scaling */
	dBSRmat *RR;  /**< reservoir block */
	
	//! neighborhood and reordering of reservoir block
	ivector *neigh; /**< neighbor information of the reservoir block */
	ivector *order; /**< ordering of the reservoir block */
	
	//! data for GS/bGS smoother for saturation block
	dBSRmat *SS; /**< saturation block */
	dvector *diaginv_S; /**< inverse of the diagonal blocks for GS/block GS smoother for saturation block */
	ivector *pivot_S; /**< pivot for the GS/block GS smoother for saturation block */
	
	//! data of AMG for pressure block
	//dCSRmat *PP; /**< pressure block */
	AMG_data *mgl_data; /**< AMG data for presure-presure block */
	//! parameters for AMG solver
	//! print level in AMG preconditioner
	SHORT print_level;
	//! max number of iterations of AMG preconditioner
	INT maxit_AMG;
	//! max number of AMG levels
	SHORT max_levels;
	//! tolerance for AMG preconditioner
	REAL amg_tol;
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
	//! switch of scaling of coarse grid correction
	SHORT coarse_scaling;
	//! degree of the polynomial used by AMLI cycle
	SHORT amli_degree;
	//! coefficients of the polynomial used by AMLI cycle
	REAL *amli_coef;
	//! relaxation parameter for smoothing the tentative prolongation
	REAL tentative_smooth;
	
	//! data of GS/bGS smoother for reservoir block
	dvector *diaginv; /**< inverse of the diagonal blocks for GS/bGS smoother for reservoir block */
	ivector *pivot; /**< pivot for the GS/bGS smoother for the reservoir matrix */
	//! data of ILU for reservoir block
	ILU_data *LU;
	
	//! data for the argumented well block
	ivector *perf_idx; /**< index of blocks which have perforation */
	ivector *perf_neigh; /**< index of blocks which are neighbors of perforations (include perforations) */
	dCSRmat *WW; /**< Argumented well block */
	//! data for direct solver for argumented well block
	void *Numeric;
	
	//! inverse of the schur complement (-I - Awr*Arr^{-1}*Arw)^{-1}, Arr may be replaced by LU
	REAL *invS;
	
	//! parameters for krylov method used for blocks
	INT  maxit; /**< max number of iterations */
	INT  restart; /**< number of iterations for restart */
	REAL tol; /**< tolerance */ 
	
	//! temporary work space
	dvector r; /**< temporary dvector used to store and restore the residual */
	REAL *w; /**<  temporary work space for other usage */
	
} precond_FASP_blkoil_data; /**< FASP precond data for Reservoir Simulation */

#endif /* end if for __FASPBLK_HEADER__ */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
