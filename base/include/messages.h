/*! \file messages.h
 *  \brief Definition of all kinds of messages, including error messages, 
 *         solver types, etc.
 *
 *  \note  This is internal use only. 
 *
 *---------------------------------------------------------------------------------
 *  Created  by Chensong Zhang on 03/20/2010.
 *  Modified by Chensong Zhang on 12/06/2011.
 *  Modified by Chensong Zhang on 12/25/2011.
 *  Modified by Chensong Zhang on 04/22/2012.
 *  Modified by Ludmil Zikatanov on 02/15/2013: CG -> SMOOTHER_CG.
 *  Modified by Chensong Zhang on 02/16/2013: GS -> SMOOTHER_GS, etc.
 *  Modified by Chensong Zhang on 04/09/2013: Add safe krylov methods.
 *  Modified by Chensong Zhang on 09/22/2013: Clean up Doxygen.
 *---------------------------------------------------------------------------------
 *
 */

#ifndef __FASP_MESSAGES__		   /*-- allow multiple inclusions --*/
#define __FASP_MESSAGES__

/**
 * \brief Definition of return status and error messages
 */
#define ERROR_OPEN_FILE       -10  /**< fail to open a file */
#define ERROR_WRONG_FILE      -11  /**< input contains wrong format */
#define ERROR_INPUT_PAR       -13  /**< wrong input argument */
#define ERROR_REGRESS         -14  /**< regression test fail */
#define ERROR_NUM_BLOCKS      -18  /**< wrong number of blocks */
#define ERROR_MISC            -19  /**< other error */
//---------------------------------------------------------------------------------
#define ERROR_ALLOC_MEM       -20  /**< fail to allocate memory */
#define ERROR_DATA_STRUCTURE  -21  /**< problem with data structures */
#define ERROR_DATA_ZERODIAG   -22  /**< matrix has zero diagonal entries */
#define ERROR_DUMMY_VAR       -23  /**< unexpected input data */
//---------------------------------------------------------------------------------
#define ERROR_AMG_INTERP_TYPE -30  /**< unknown interpolation type */
#define ERROR_AMG_SMOOTH_TYPE -31  /**< unknown smoother type */
#define ERROR_AMG_COARSE_TYPE -32  /**< unknown coarsening type */
#define ERROR_AMG_COARSEING   -33  /**< coarsening step failed to complete */
//---------------------------------------------------------------------------------
#define ERROR_SOLVER_TYPE     -40  /**< unknown solver type */
#define ERROR_SOLVER_PRECTYPE -41  /**< unknow precond type */
#define ERROR_SOLVER_STAG     -42  /**< solver stagnates */
#define ERROR_SOLVER_SOLSTAG  -43  /**< solver's solution is too small */
#define ERROR_SOLVER_TOLSMALL -44  /**< solver's tolerance is too small */
#define ERROR_SOLVER_ILUSETUP -45  /**< ILU setup error */
#define ERROR_SOLVER_MISC     -46  /**< misc solver error during run time */
#define ERROR_SOLVER_MAXIT    -48  /**< maximal iteration number exceeded */
#define ERROR_SOLVER_EXIT     -49  /**< solver does not quit successfully */
//---------------------------------------------------------------------------------
#define ERROR_QUAD_TYPE       -60  /**< unknown quadrature type */
#define ERROR_QUAD_DIM        -61  /**< unsupported quadrature dim */
//---------------------------------------------------------------------------------
#define ERROR_LIC_TYPE        -80  /**< wrong license type */
//---------------------------------------------------------------------------------
#define RUN_FAIL              -99  /**< general failure */
#define SUCCESS                 0  /**< return from funtion successfully */

/**
 * \brief Definition of logic type  
 */
#define TRUE                    1  /**< logic TRUE */
#define FALSE                   0  /**< logic FALSE */

/** 
 * \brief Definition of switch  
 */
#define ON                      1  /**< turn on certain parameter */
#define OFF                     0  /**< turn off certain parameter */

/** 
 * \brief Print level for all subroutines -- not including DEBUG output
 */
#define PRINT_NONE              0  /**< slient: no printout at all */
#define PRINT_MIN               1  /**< quiet: minimal print, like convergence */
#define PRINT_SOME              2  /**< some: print cpu time, iteration number */
#define PRINT_MORE              4  /**< more: print more useful information */
#define PRINT_MOST              8  /**< most: maximal printouts, no disk files*/
#define PRINT_ALL              10  /**< everything: all printouts allowed */

/**
 * \brief Definition of matrix format
 **/
#define MAT_FREE                0  /**< matrix-free format: only mxv action */
#define MAT_CSR                 1  /**< compressed sparse row */
#define MAT_BSR                 2  /**< blockwise compressed sparse row */
#define MAT_STR                 3  /**< structured sparse matrix */
#define MAT_bCSR                4  /**< block matrix of CSR */
#define MAT_bBSR                5  /**< block matrix of BSR for borded systems */
#define MAT_CSRL                6  /**< modified CSR to reduce cache missing */
#define MAT_SymCSR              7  /**< symmetric CSR format */

/**
 * \brief Definition of solver types for iterative methods
 */
#define SOLVER_CG               1  /**< Conjugate Gradient */
#define SOLVER_BiCGstab         2  /**< Biconjugate Gradient Stabilized */
#define SOLVER_MinRes           3  /**< Minimal Residual */
#define SOLVER_GMRES            4  /**< Generalized Minimal Residual */
#define SOLVER_VGMRES           5  /**< Variable Restarting GMRES */
#define SOLVER_VFGMRES          6  /**< Variable Restarting Flexible GMRES */
#define SOLVER_GCG              7  /**< Generalized Conjugate Gradient */
//---------------------------------------------------------------------------------
#define SOLVER_SCG             11  /**< Conjugate Gradient with safe net */
#define SOLVER_SBiCGstab       12  /**< BiCGstab with safe net */
#define SOLVER_SMinRes         13  /**< MinRes with safe net */
#define SOLVER_SGMRES          14  /**< GMRes with safe net */
#define SOLVER_SVGMRES         15  /**< Variable-restart GMRES with safe net */
#define SOLVER_SVFGMRES        16  /**< Variable-restart FGMRES with safe net */
#define SOLVER_SGCG            17  /**< GCG with safe net */
//---------------------------------------------------------------------------------
#define SOLVER_AMG             21  /**< AMG as an iterative solver */
#define SOLVER_FMG             22  /**< Full AMG as an solver */

/** 
 * \brief Definition of solver types for direct methods (requires external libs)
 */
#define SOLVER_SUPERLU         31  /**< SuperLU Direct Solver */
#define SOLVER_UMFPACK         32  /**< UMFPack Direct Solver */
#define SOLVER_MUMPS           33  /**< MUMPS   Direct Solver */

/** 
 * \brief Definition of iterative solver stopping criteria types
 */
#define STOP_REL_RES            1  /**< relative residual ||r||/||b|| */
#define STOP_REL_PRECRES        2  /**< relative B-residual ||r||_B/||b||_B */
#define STOP_MOD_REL_RES        3  /**< modified relative residual ||r||/||x|| */

/** 
 * \brief Definition of preconditioner type for iterative methods
 */
#define PREC_NULL               0  /**< with no precond */
#define PREC_DIAG               1  /**< with diagonal precond */
#define PREC_AMG                2  /**< with AMG precond */
#define PREC_FMG                3  /**< with full AMG precond */
#define PREC_ILU                4  /**< with ILU precond */
#define PREC_SCHWARZ            5  /**< with Schwarz preconditioner */

/**
 * \brief Type of ILU methods
 */
#define ILUk                    1  /**< ILUk */
#define ILUt                    2  /**< ILUt */
#define ILUtp                   3  /**< ILUtp */

/**
 * \brief Definition of AMG types
 */
#define CLASSIC_AMG             1  /**< classic AMG */
#define SA_AMG                  2  /**< smoothed aggregation AMG */
#define UA_AMG                  3  /**< unsmoothed aggregation AMG */

/**
 * \brief Definition of aggregation types
 */
#define PAIRWISE                1  /**< pairwise aggregation AMG */
#define VMB                     2  /**< VMB aggregation AMG */

/**
 * \brief Definition of cycle types
 */
#define V_CYCLE	                1  /**< V-cycle */
#define W_CYCLE                 2  /**< W-cycle */
#define AMLI_CYCLE              3  /**< AMLI-cycle */
#define NL_AMLI_CYCLE           4  /**< Nonlinear AMLI-cycle */

/** 
 * \brief Definition of standard smoother types
 */
#define SMOOTHER_JACOBI         1  /**< Jacobi smoother */
#define SMOOTHER_GS             2  /**< Gauss-Seidel smoother */
#define SMOOTHER_SGS            3  /**< symm Gauss-Seidel smoother */
#define SMOOTHER_CG             4  /**< CG as a smoother */
#define SMOOTHER_SOR            5  /**< SOR smoother */
#define SMOOTHER_SSOR           6  /**< SSOR smoother */
#define SMOOTHER_GSOR           7  /**< GS + SOR smoother */
#define SMOOTHER_SGSOR          8  /**< SGS + SSOR smoother */
#define SMOOTHER_POLY           9  /**< Polynomial smoother */
#define SMOOTHER_L1DIAG        10  /**< L1 norm diagonal scaling smoother */

/**
 * \brief Definition of specialized smoother types
 */
#define SMOOTHER_BLKOIL        11  /**< Used in monolithic AMG for black-oil */
#define SMOOTHER_SPETEN        19  /**< Used in monolithic AMG for black-oil */

/**
 * \brief Definition of coarsening types
 */
#define COARSE_RS               1  /**< Classical coarsening */
#define COARSE_CR               3  /**< Compatible relaxation */
#define COARSE_AC               4  /**< Aggressive coarsening */

/**
 * \brief Definition of interpolation types
 */
#define INTERP_DIR              1  /**< Direct interpolation */
#define INTERP_STD              2  /**< Standard interpolation */
#define INTERP_ENG              3  /**< energy minimization interp in C */

/** 
 * \brief Type of vertices (dofs) for C/F splitting
 */
#define UNPT                   -1  /**< undetermined points */
#define FGPT                    0  /**< fine grid points  */
#define CGPT                    1  /**< coarse grid points */
#define ISPT                    2  /**< isolated points */

/** 
 * \brief Definition of smoothing order
 */
#define NO_ORDER                0  /**< Natural order smoothing */
#define CF_ORDER                1  /**< C/F order smoothing */

/** 
 * \brief Type of ordering for smoothers
 */
#define USERDEFINED             0  /**< USERDEFINED order */
#define CPFIRST                 1  /**< C-points first order */
#define FPFIRST                -1  /**< F-points first order */
#define ASCEND                 12  /**< Asscending order */
#define DESCEND                21  /**< Dsscending order */

#endif                             /* end if for __FASP_MESSAGES__ */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
