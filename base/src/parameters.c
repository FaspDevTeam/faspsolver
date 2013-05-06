/*! \file parameters.c
 *  \brief Initialize, set, or print input data and parameters
 *
 */

#include <stdio.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_param_init (char *inputfile,
 *                           input_param *inparam,
 *                           itsolver_param *itsparam,
 *                           AMG_param *amgparam,
 *                           ILU_param *iluparam,
 *                           Schwarz_param *schparam)
 *
 * \brief Initialize parameters, global variables, etc
 *
 * \param inputfile     Filename of the input file
 * \param inparam       Input parameters
 * \param itsparam      Iterative solver parameters
 * \param amgparam      AMG parameters
 * \param iluparam      ILU parameters
 * \param schparam      Schwarz parameters
 *
 * \author Chensong Zhang
 * \date   2010/08/12
 *
 * Modified by Xiaozhe Hu (01/23/2011): initialize, then set value
 * Modified by Chensong Zhang (09/12/2012): find a bug during debugging in VS08
 */
void fasp_param_init (char *inputfile,
                      input_param *inparam,
                      itsolver_param *itsparam,
                      AMG_param *amgparam,
                      ILU_param *iluparam,
                      Schwarz_param *schparam)
{
    total_alloc_mem   = 0; // initialize total memeory amount
    total_alloc_count = 0; // initialize alloc count
    
    if (itsparam) fasp_param_solver_init(itsparam);
    
    if (amgparam) fasp_param_amg_init(amgparam);
    
    if (iluparam) fasp_param_ilu_init(iluparam);
    
    if (schparam) fasp_param_schwarz_init(schparam);
    
    if (inputfile) {
        fasp_param_input(inputfile,inparam);
        if (itsparam) fasp_param_solver_set(itsparam,inparam);
        if (amgparam) fasp_param_amg_set(amgparam,inparam);
        if (iluparam) fasp_param_ilu_set(iluparam,inparam);
        if (schparam) fasp_param_schwarz_set(schparam,inparam);
    }
    else {
        printf("### WARNING: No input specified. Use default values instead!\n");
    }
}

/**
 * \fn void fasp_param_input_init (input_param *inparam)
 *
 * \brief Initialize input parameters
 *
 * \param inparam    Input parameters
 *
 * \author Chensong Zhang
 * \date   2010/03/20
 */
void fasp_param_input_init (input_param *inparam)
{
    strcpy(inparam->workdir,"data/");

    // Input/output
    inparam->print_level              = PRINT_MIN;
    inparam->output_type              = 0;
    
    // Problem information
    inparam->problem_num              = 10;
    inparam->solver_type              = SOLVER_CG;
    inparam->precond_type             = PREC_AMG;
    inparam->stop_type                = STOP_REL_RES;
    
    // Solver parameters
    inparam->itsolver_tol             = 1e-6;
    inparam->itsolver_maxit           = 500;
    inparam->restart                  = 25;
    
    // ILU method parameters
    inparam->ILU_type                 = ILUk;
    inparam->ILU_lfil                 = 0;
    inparam->ILU_droptol              = 0.001;
    inparam->ILU_relax                = 0;
    inparam->ILU_permtol              = 0.0;
    
    // Schwarz method parameters
    inparam->Schwarz_mmsize           = 200;
	inparam->Schwarz_maxlvl           = 2;
	inparam->Schwarz_type             = 1;
    
    // AMG method parameters
    inparam->AMG_type                 = CLASSIC_AMG;
    inparam->AMG_levels               = 15;
    inparam->AMG_cycle_type           = V_CYCLE;
    inparam->AMG_smoother             = SMOOTHER_GS;
    inparam->AMG_presmooth_iter       = 2;
    inparam->AMG_postsmooth_iter      = 2;
    inparam->AMG_relaxation           = 1.0;
    inparam->AMG_coarse_dof           = 500;
    inparam->AMG_tol                  = 1e-4*inparam->itsolver_tol;
    inparam->AMG_maxit                = 1;
    inparam->AMG_ILU_levels           = 0;
    inparam->AMG_schwarz_levels       = 0;
    inparam->AMG_coarse_scaling       = OFF; // Keep this off. Require investigation --Chensong
    inparam->AMG_amli_degree          = 1;
    inparam->AMG_nl_amli_krylov_type  = 2;
    
    // Classical AMG specific
    inparam->AMG_coarsening_type      = 1;
    inparam->AMG_interpolation_type   = 1;
    inparam->AMG_strong_threshold     = 0.5;
    inparam->AMG_truncation_threshold = 0.4;
    inparam->AMG_max_row_sum          = 0.9;
    inparam->AMG_aggressive_level     = 0;
    inparam->AMG_aggressive_path      = 1;
    
    // Aggregation AMG specific
    inparam->AMG_strong_coupled       = 0.08;
    inparam->AMG_max_aggregation      = 9;
    inparam->AMG_tentative_smooth     = 0.67;
    inparam->AMG_smooth_filter        = ON;
}

/**
 * \fn void fasp_param_amg_init (AMG_param *amgparam)
 *
 * \brief Initialize AMG parameters
 *
 * \param amgparam    Parameters for AMG
 *
 * \author Chensong Zhang
 * \date   2010/04/03
 */
void fasp_param_amg_init (AMG_param *amgparam)
{
    // AMG type 
    amgparam->AMG_type             = CLASSIC_AMG;
    amgparam->print_level          = PRINT_NONE;
    amgparam->maxit                = 1;
    amgparam->tol                  = 1e-8;
    
    // AMG method parameters
    amgparam->max_levels           = 15;
    amgparam->coarse_dof           = 500;
    amgparam->cycle_type           = V_CYCLE;
    amgparam->smoother             = SMOOTHER_GS;
    amgparam->smooth_order         = CF_ORDER;
    amgparam->presmooth_iter       = 2;
    amgparam->postsmooth_iter      = 2;
    amgparam->relaxation           = 1.0;
    amgparam->polynomial_degree    = 3;
    amgparam->coarse_scaling       = OFF; // Keep this off. Require investigation --Chensong
    amgparam->amli_degree          = 1;
    amgparam->amli_coef            = NULL;
    amgparam->nl_amli_krylov_type  = 2;
    
    // Classical AMG specific
    amgparam->coarsening_type      = 1;
    amgparam->interpolation_type   = 1;
    amgparam->strong_threshold     = 0.5;
    amgparam->truncation_threshold = 0.4;
    amgparam->max_row_sum          = 0.9;
    amgparam->aggressive_level     = 0;
    amgparam->aggressive_path      = 1;
    
    // Aggregation AMG specific
    amgparam->strong_coupled       = 0.08;
    amgparam->max_aggregation      = 9;
    amgparam->tentative_smooth     = 0.0;
    amgparam->smooth_filter        = OFF;
    
    // ILU smoother parameters
    amgparam->ILU_type             = ILUk;
    amgparam->ILU_levels           = 0;
    amgparam->ILU_lfil             = 0;
    amgparam->ILU_droptol          = 0.001;
    amgparam->ILU_relax            = 0;
    
    // Schwarz smoother parameters
    amgparam->schwarz_levels       = 0;
    amgparam->schwarz_mmsize       = 200;
    amgparam->schwarz_maxlvl       = 2;
    amgparam->schwarz_type         = 1;
}

/**
 * \fn void fasp_param_solver_init (itsolver_param *itsparam)
 *
 * \brief Initialize itsolver_param
 *
 * \param itsparam   Parameters for iterative solvers
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_param_solver_init (itsolver_param *itsparam)
{
    itsparam->print_level   = 0;
    itsparam->itsolver_type = SOLVER_CG;
    itsparam->precond_type  = PREC_AMG;
    itsparam->stop_type     = STOP_REL_RES;
    itsparam->maxit         = 500;
    itsparam->restart       = 25;
    itsparam->tol           = 1e-8;
}

/**
 * \fn void fasp_param_ilu_init (ILU_param *iluparam)
 *
 * \brief Initialize ILU parameters
 *
 * \param iluparam  Parameters for ILU
 *
 * \author Chensong Zhang
 * \date   2010/04/06
 */
void fasp_param_ilu_init (ILU_param *iluparam)
{
    iluparam->print_level  = 0;
    iluparam->ILU_type     = ILUk;
    iluparam->ILU_lfil     = 2;
    iluparam->ILU_droptol  = 0.001;
    iluparam->ILU_relax    = 0;
    iluparam->ILU_permtol  = 0.01;
}

/**
 * \fn void fasp_param_schwarz_init (Schwarz_param *schparam)
 *
 * \brief Initialize Schwarz parameters
 *
 * \param schparam    Parameters for Schwarz method
 *
 * \author Xiaozhe Hu
 * \date   05/22/2012
 */
void fasp_param_schwarz_init (Schwarz_param *schparam)
{
    schparam->print_level    = 0;
    schparam->schwarz_type   = 3;
    schparam->schwarz_maxlvl = 2;
    schparam->schwarz_mmsize = 200;
}

/**
 * \fn void fasp_param_amg_set (AMG_param *param, input_param *inparam)
 *
 * \brief Set AMG_param from INPUT
 *
 * \param param     Parameters for AMG
 * \param inparam   Input parameters
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_param_amg_set (AMG_param *param,
                         input_param *inparam)
{
    param->AMG_type    = inparam->AMG_type;
    param->print_level = inparam->print_level;
    
    if (inparam->solver_type == SOLVER_AMG) {
        param->maxit = inparam->itsolver_maxit;
        param->tol   = inparam->itsolver_tol;
    }
    else if (inparam->solver_type == SOLVER_FMG) {
        param->maxit = inparam->itsolver_maxit;
        param->tol   = inparam->itsolver_tol;
    }
    else {
        param->maxit = inparam->AMG_maxit;
        param->tol   = inparam->AMG_tol;
    }
    
    param->max_levels           = inparam->AMG_levels;
    param->cycle_type           = inparam->AMG_cycle_type;
    param->smoother             = inparam->AMG_smoother;
    param->relaxation           = inparam->AMG_relaxation;
    param->polynomial_degree    = inparam->AMG_polynomial_degree;
    param->presmooth_iter       = inparam->AMG_presmooth_iter;
    param->postsmooth_iter      = inparam->AMG_postsmooth_iter;
    param->coarse_dof           = inparam->AMG_coarse_dof;
    param->coarse_scaling       = inparam->AMG_coarse_scaling;
    param->amli_degree          = inparam->AMG_amli_degree;
    param->amli_coef            = NULL;
    param->nl_amli_krylov_type  = inparam->AMG_nl_amli_krylov_type;
    
    param->coarsening_type      = inparam->AMG_coarsening_type;
    param->interpolation_type   = inparam->AMG_interpolation_type;
    param->strong_threshold     = inparam->AMG_strong_threshold;
    param->truncation_threshold = inparam->AMG_truncation_threshold;
    param->max_row_sum          = inparam->AMG_max_row_sum;
    param->aggressive_level     = inparam->AMG_aggressive_level;
    param->aggressive_path      = inparam->AMG_aggressive_path;
    
    param->strong_coupled       = inparam->AMG_strong_coupled;
    param->max_aggregation      = inparam->AMG_max_aggregation;
    param->tentative_smooth     = inparam->AMG_tentative_smooth;
    param->smooth_filter        = inparam->AMG_smooth_filter;
    
    param->ILU_levels           = inparam->AMG_ILU_levels;
    param->ILU_type             = inparam->ILU_type;
    param->ILU_lfil             = inparam->ILU_lfil;
    param->ILU_droptol          = inparam->ILU_droptol;
    param->ILU_relax            = inparam->ILU_relax;
    param->ILU_permtol          = inparam->ILU_permtol;
    
    param->schwarz_levels       = inparam->AMG_schwarz_levels;
	param->schwarz_mmsize       = inparam->Schwarz_mmsize;
	param->schwarz_maxlvl       = inparam->Schwarz_maxlvl;
	param->schwarz_type         = inparam->Schwarz_type;
}

/**
 * \fn void fasp_param_ilu_set (ILU_param *iluparam, input_param *inparam)
 *
 * \brief Set ILU_param with INPUT
 *
 * \param iluparam    Parameters for ILU
 * \param inparam     Input parameters
 *
 * \author Chensong Zhang
 * \date   2010/04/03
 */
void fasp_param_ilu_set (ILU_param *iluparam,
                         input_param *inparam)
{
    iluparam->print_level = inparam->print_level;
    iluparam->ILU_type    = inparam->ILU_type;
    iluparam->ILU_lfil    = inparam->ILU_lfil;
    iluparam->ILU_droptol = inparam->ILU_droptol;
    iluparam->ILU_relax   = inparam->ILU_relax;
    iluparam->ILU_permtol = inparam->ILU_permtol;
}

/**
 * \fn void fasp_param_schwarz_set (Schwarz_param *schparam, input_param *inparam)
 *
 * \brief Set Schwarz_param with INPUT
 *
 * \param schparam    Parameters for Schwarz method
 * \param inparam     Input parameters
 *
 * \author Xiaozhe Hu
 * \date   05/22/2012
 */
void fasp_param_schwarz_set (Schwarz_param *schparam,
                             input_param *inparam)
{
    schparam->print_level    = inparam->print_level;
    schparam->schwarz_type   = inparam->Schwarz_type;
    schparam->schwarz_maxlvl = inparam->Schwarz_maxlvl;
    schparam->schwarz_mmsize = inparam->Schwarz_mmsize;
}

/**
 * \fn void fasp_param_solver_set (itsolver_param *itsparam, input_param *inparam)
 *
 * \brief Set itsolver_param with INPUT
 *
 * \param itsparam   Parameters for iterative solvers
 * \param inparam    Input parameters
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_param_solver_set (itsolver_param *itsparam,
                            input_param *inparam)
{
    itsparam->print_level    = inparam->print_level;
    itsparam->itsolver_type  = inparam->solver_type;
    itsparam->precond_type   = inparam->precond_type;
    itsparam->stop_type      = inparam->stop_type;
    itsparam->restart        = inparam->restart;
    
    if (itsparam->itsolver_type == SOLVER_AMG) {
        itsparam->tol   = inparam->AMG_tol;
        itsparam->maxit = inparam->AMG_maxit;
    }
    else {
        itsparam->tol   = inparam->itsolver_tol;
        itsparam->maxit = inparam->itsolver_maxit;
    }
}

/**
 * \fn void fasp_precond_data_null (precond_data *pcdata)
 *
 * \brief Initialize precond_data
 *
 * \param pcdata   Preconditioning data structure
 *
 * \author Chensong Zhang
 * \date   2010/03/23
 */
void fasp_precond_data_null (precond_data *pcdata)
{
    pcdata->AMG_type            = CLASSIC_AMG;
    pcdata->print_level         = PRINT_MIN;
    pcdata->maxit               = 500;
    pcdata->max_levels          = 20;
    pcdata->tol                 = 1e-8;
    pcdata->cycle_type          = V_CYCLE;
    pcdata->smoother            = SMOOTHER_GS;
    pcdata->smooth_order        = CF_ORDER;
    pcdata->presmooth_iter      = 2;
    pcdata->postsmooth_iter     = 2;
    pcdata->relaxation          = 1.1;
    pcdata->coarsening_type     = 1;
    pcdata->coarse_scaling      = ON;
    pcdata->amli_degree         = 1;
    pcdata->nl_amli_krylov_type = SOLVER_GCG;
}

/**
 * \fn void fasp_param_amg_to_prec (precond_data *pcdata, AMG_param *amgparam)
 *
 * \brief Set precond_data with AMG_param
 *
 * \param pcdata      Preconditioning data structure
 * \param amgparam    Parameters for AMG
 *
 * \author Chensong Zhang
 * \date   2011/01/10
 */
void fasp_param_amg_to_prec (precond_data *pcdata,
                             AMG_param *amgparam)
{
    pcdata->AMG_type            = amgparam->AMG_type;
    pcdata->print_level         = amgparam->print_level;
    pcdata->maxit               = amgparam->maxit;
    pcdata->max_levels          = amgparam->max_levels;
    pcdata->tol                 = amgparam->tol;
    pcdata->cycle_type          = amgparam->cycle_type;
    pcdata->smoother            = amgparam->smoother;
    pcdata->smooth_order        = amgparam->smooth_order;
    pcdata->presmooth_iter      = amgparam->presmooth_iter;
    pcdata->postsmooth_iter     = amgparam->postsmooth_iter;
    pcdata->coarsening_type     = amgparam->coarsening_type;
    pcdata->relaxation          = amgparam->relaxation;
    pcdata->polynomial_degree   = amgparam->polynomial_degree;
    pcdata->coarse_scaling      = amgparam->coarse_scaling;
    pcdata->amli_degree         = amgparam->amli_degree;
    pcdata->amli_coef           = amgparam->amli_coef;
    pcdata->nl_amli_krylov_type = amgparam->nl_amli_krylov_type;
    pcdata->tentative_smooth    = amgparam->tentative_smooth;
}

/**
 * \fn void fasp_param_prec_to_amg (AMG_param *amgparam, precond_data *pcdata)
 *
 * \brief Set AMG_param with precond_data
 *
 * \param amgparam    Parameters for AMG
 * \param pcdata      Preconditioning data structure
 *
 * \author Chensong Zhang
 * \date   2011/01/10
 */
void fasp_param_prec_to_amg (AMG_param *amgparam,
                             precond_data *pcdata)
{
    amgparam->AMG_type            = pcdata->AMG_type;
    amgparam->print_level         = pcdata->print_level;
    amgparam->cycle_type          = pcdata->cycle_type;
    amgparam->smoother            = pcdata->smoother;
    amgparam->smooth_order        = pcdata->smooth_order;
    amgparam->presmooth_iter      = pcdata->presmooth_iter;
    amgparam->postsmooth_iter     = pcdata->postsmooth_iter;
    amgparam->relaxation          = pcdata->relaxation;
    amgparam->polynomial_degree   = pcdata->polynomial_degree;
    amgparam->coarse_scaling      = pcdata->coarse_scaling;
    amgparam->amli_degree         = pcdata->amli_degree;
    amgparam->amli_coef           = pcdata->amli_coef;
    amgparam->nl_amli_krylov_type = pcdata->nl_amli_krylov_type;
    amgparam->tentative_smooth    = pcdata->tentative_smooth;
    amgparam->ILU_levels          = pcdata->mgl_data->ILU_levels;
}

/**
 * \fn void fasp_param_amg_to_prec_bsr (precond_data_bsr *pcdata, AMG_param *amgparam)
 *
 * \brief Set precond_data_bsr with AMG_param
 *
 * \param pcdata      Preconditioning data structure
 * \param amgparam    Parameters for AMG
 *
 * \author Xiaozhe Hu
 * \date   02/06/2012
 */
void fasp_param_amg_to_prec_bsr (precond_data_bsr *pcdata,
                                 AMG_param *amgparam)
{
    pcdata->AMG_type            = amgparam->AMG_type;
    pcdata->print_level         = amgparam->print_level;
    pcdata->maxit               = amgparam->maxit;
    pcdata->max_levels          = amgparam->max_levels;
    pcdata->tol                 = amgparam->tol;
    pcdata->cycle_type          = amgparam->cycle_type;
    pcdata->smoother            = amgparam->smoother;
    pcdata->smooth_order        = amgparam->smooth_order;
    pcdata->presmooth_iter      = amgparam->presmooth_iter;
    pcdata->postsmooth_iter     = amgparam->postsmooth_iter;
    pcdata->coarsening_type     = amgparam->coarsening_type;
    pcdata->relaxation          = amgparam->relaxation;
    pcdata->coarse_scaling      = amgparam->coarse_scaling;
    pcdata->amli_degree         = amgparam->amli_degree;
    pcdata->amli_coef           = amgparam->amli_coef;
    pcdata->nl_amli_krylov_type = amgparam->nl_amli_krylov_type;
    pcdata->tentative_smooth    = amgparam->tentative_smooth;
}

/**
 * \fn void fasp_param_prec_to_amg_bsr (AMG_param *amgparam, precond_data_bsr *pcdata)
 *
 * \brief Set AMG_param with precond_data
 *
 * \param amgparam    Parameters for AMG
 * \param pcdata      Preconditioning data structure
 *
 * \author Xiaozhe Hu
 * \date   02/06/2012
 */
void fasp_param_prec_to_amg_bsr (AMG_param *amgparam,
                                 precond_data_bsr *pcdata)
{
    amgparam->AMG_type            = pcdata->AMG_type;
    amgparam->print_level         = pcdata->print_level;
    amgparam->cycle_type          = pcdata->cycle_type;
    amgparam->smoother            = pcdata->smoother;
    amgparam->smooth_order        = pcdata->smooth_order;
    amgparam->presmooth_iter      = pcdata->presmooth_iter;
    amgparam->postsmooth_iter     = pcdata->postsmooth_iter;
    amgparam->relaxation          = pcdata->relaxation;
    amgparam->coarse_scaling      = pcdata->coarse_scaling;
    amgparam->amli_degree         = pcdata->amli_degree;
    amgparam->amli_coef           = pcdata->amli_coef;
    amgparam->nl_amli_krylov_type = pcdata->nl_amli_krylov_type;
    amgparam->tentative_smooth    = pcdata->tentative_smooth;
    amgparam->ILU_levels          = pcdata->mgl_data->ILU_levels;
}

/**
 * \fn void fasp_param_amg_print (AMG_param *param)
 *
 * \brief Print out AMG parameters
 *
 * \param param   Parameters for AMG
 *
 * \author Chensong Zhang
 * \date   2010/03/22
 */
void fasp_param_amg_print (AMG_param *param)
{
    
    if ( param ) {
        
        printf("\n       Parameters in AMG_param\n");
        printf("-----------------------------------------------\n");
        
        printf("AMG print level:                   %d\n", param->print_level);
        printf("AMG max num of iter:               %d\n", param->maxit);
        printf("AMG type:                          %d\n", param->AMG_type);
        printf("AMG tolerance:                     %.2e\n", param->tol);
        printf("AMG max levels:                    %d\n", param->max_levels);
        printf("AMG cycle type:                    %d\n", param->cycle_type);
        printf("AMG scaling of coarse correction:  %d\n", param->coarse_scaling);
        printf("AMG smoother type:                 %d\n", param->smoother);
        printf("AMG smoother order:                %d\n", param->smooth_order);
        printf("AMG num of presmoothing:           %d\n", param->presmooth_iter);
        printf("AMG num of postsmoothing:          %d\n", param->postsmooth_iter);
        
        if ( param->smoother == SMOOTHER_SOR  ||
             param->smoother == SMOOTHER_SSOR ||
             param->smoother == SMOOTHER_GSOR ||
             param->smoother == SMOOTHER_SGSOR ) {
            printf("AMG relax factor:                  %.4f\n", param->relaxation);
        }
        
        if ( param->smoother == SMOOTHER_POLY ) {
            printf("AMG polynomial smoother degree:    %d\n", param->polynomial_degree);
        }
        
        if ( param->cycle_type == AMLI_CYCLE ) {
            printf("AMG AMLI degree of polynomial:     %d\n", param->amli_degree);
        }
        
        if ( param->cycle_type == NL_AMLI_CYCLE ) {
            printf("AMG Nonlinear AMLI Krylov type:    %d\n", param->nl_amli_krylov_type);
        }
        
        switch (param->AMG_type) {
            case CLASSIC_AMG:
                printf("AMG coarsening type:               %d\n", param->coarsening_type);
                printf("AMG interpolation type:            %d\n", param->interpolation_type);
                printf("AMG dof on coarsest grid:          %d\n", param->coarse_dof);
                printf("AMG strong threshold:              %.4f\n", param->strong_threshold);
                printf("AMG truncation threshold:          %.4f\n", param->truncation_threshold);
                printf("AMG max row sum:                   %.4f\n", param->max_row_sum);
                printf("AMG aggressive levels:             %d\n", param->aggressive_level);
                printf("AMG aggressive path:               %d\n", param->aggressive_path);
                break;
                
            default: // SA_AMG or UA_AMG
                printf("Aggregation AMG strong coupling:   %.4f\n", param->strong_coupled);
                printf("Aggregation AMG max aggregation:   %d\n", param->max_aggregation);
                printf("Aggregation AMG tentative smooth:  %.4f\n", param->tentative_smooth);
                printf("Aggregation AMG smooth filter:     %d\n", param->smooth_filter);
                break;
        }
        
        if (param->ILU_levels>0) {
            printf("AMG ILU type:                      %d\n", param->ILU_type);
            printf("AMG ILU level:                     %d\n", param->ILU_levels);
            printf("AMG ILU level of fill-in:          %d\n", param->ILU_lfil);
            printf("AMG ILU drop tol:                  %e\n", param->ILU_droptol);
            printf("AMG ILU relaxation:                %f\n", param->ILU_relax);
        }
        
        if (param->schwarz_levels>0){
            printf("AMG Schwarz type:                  %d\n", param->schwarz_type);
            printf("AMG Schwarz level:                 %d\n", param->schwarz_levels);
            printf("AMG Schwarz level of forming block:%d\n", param->schwarz_maxlvl);
            printf("AMG Schwarz maximal block size:    %d\n", param->schwarz_mmsize);
        }
        
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    } // end if (param)
    
}

/**
 * \fn void fasp_param_ilu_print (ILU_param *param)
 *
 * \brief Print out ILU parameters
 *
 * \param param    Parameters for ILU
 *
 * \author Chensong Zhang
 * \date   2011/12/20
 */
void fasp_param_ilu_print (ILU_param *param)
{
    if ( param ) {
        
        printf("\n       Parameters in ILU_param\n");
        printf("-----------------------------------------------\n");
        printf("ILU print level:                   %d\n",   param->print_level);
        printf("ILU type:                          %d\n",   param->ILU_type);
        printf("ILU level of fill-in:              %d\n",   param->ILU_lfil);
        printf("ILU relaxation factor:             %.4f\n", param->ILU_relax);
        printf("ILU drop tolerance:                %.2e\n", param->ILU_droptol);
        printf("ILU permutation tolerance:         %.2e\n", param->ILU_permtol);
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    }
}

/**
 * \fn void fasp_param_schwarz_print (Schwarz_param *param)
 *
 * \brief Print out Schwarz parameters
 *
 * \param param    Parameters for Schwarz
 *
 * \author Xiaozhe Hu
 * \date   05/22/2012
 */
void fasp_param_schwarz_print (Schwarz_param *param)
{
    if ( param ) {
        
        printf("\n       Parameters in Schwarz_param\n");
        printf("-----------------------------------------------\n");
        printf("Schwarz print level:               %d\n",   param->print_level);
        printf("Schwarz type:                      %d\n",   param->schwarz_type);
        printf("Schwarz level of forming blokcs:   %d\n",   param->schwarz_maxlvl);
        printf("Schwarz maximal block size:        %d\n",   param->schwarz_mmsize);
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    }
}


/**
 * \fn void fasp_param_solver_print (itsolver_param *param)
 *
 * \brief Print out itsolver parameters
 *
 * \param param    Paramters for iterative solvers
 *
 * \author Chensong Zhang
 * \date   2011/12/20
 */
void fasp_param_solver_print (itsolver_param *param)
{
    if ( param ) {
        
        printf("\n       Parameters in itsolver_param\n");
        printf("-----------------------------------------------\n");
        
        printf("Solver print level:                %d\n", param->print_level);
        printf("Solver type:                       %d\n", param->itsolver_type);
        printf("Solver precond type:               %d\n", param->precond_type);
        printf("Solver max num of iter:            %d\n", param->maxit);
        printf("Solver tolerance:                  %.2e\n", param->tol);
        printf("Solver stopping type:              %d\n", param->stop_type);
        
        if (param->itsolver_type==SOLVER_GMRES ||
            param->itsolver_type==SOLVER_VGMRES) {
            printf("Solver restart number:             %d\n", param->restart);
        }
        
        printf("-----------------------------------------------\n\n");
        
    }
    else {
        printf("### WARNING: param has not been set!\n");
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
