/*! \file input.c
 *
 *  \brief Read input parameters
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn SHORT fasp_param_check (input_param *inparam)
 *
 * \brief Simple check on input parameters
 *
 * \param inparam     Input parameters
 *
 * \author Chensong Zhang
 * \date   09/29/2013
 */
SHORT fasp_param_check (input_param *inparam)
{
    SHORT status = SUCCESS;
    
    if ( inparam->problem_num<0
        || inparam->solver_type<0
        || inparam->solver_type>50
        || inparam->precond_type<0
        || inparam->itsolver_tol<=0
        || inparam->itsolver_maxit<=0
        || inparam->stop_type<=0
        || inparam->stop_type>3
        || inparam->restart<0
        || inparam->ILU_type<=0
        || inparam->ILU_type>3
        || inparam->ILU_lfil<0
        || inparam->ILU_droptol<=0
        || inparam->ILU_relax<0
        || inparam->ILU_permtol<0
        || inparam->Schwarz_mmsize<0
        || inparam->Schwarz_maxlvl<0
        || inparam->Schwarz_type<0
        || inparam->AMG_type<=0
        || inparam->AMG_type>3
        || inparam->AMG_cycle_type<=0
        || inparam->AMG_cycle_type>4
        || inparam->AMG_levels<0
        || inparam->AMG_ILU_levels<0
        || inparam->AMG_coarse_dof<=0
        || inparam->AMG_tol<0
        || inparam->AMG_maxit<0
        || inparam->AMG_coarsening_type<=0
        || inparam->AMG_coarsening_type>4
        || inparam->AMG_interpolation_type<0
        || inparam->AMG_interpolation_type>5
        || inparam->AMG_smoother<0
        || inparam->AMG_smoother>20
        || inparam->AMG_strong_threshold<0.0
        || inparam->AMG_strong_threshold>0.9999
        || inparam->AMG_truncation_threshold<0.0
        || inparam->AMG_truncation_threshold>0.9999
        || inparam->AMG_max_row_sum<0.0
        || inparam->AMG_presmooth_iter<0
        || inparam->AMG_postsmooth_iter<0
        || inparam->AMG_amli_degree<0
        || inparam->AMG_aggressive_level<0
        || inparam->AMG_aggressive_path<0
        || inparam->AMG_strong_coupled<0
        || inparam->AMG_max_aggregation<=0
        || inparam->AMG_tentative_smooth<0
        || inparam->AMG_smooth_filter<0
        ) status = ERROR_INPUT_PAR;
    
    return status;
}

/**
 * \fn void fasp_param_input (char *filenm, input_param *inparam)
 *
 * \brief Read input parameters from disk file
 *
 * \param filenm    File name for input file
 * \param inparam   Input parameters
 *
 * \author Chensong Zhang
 * \date   03/20/2010 
 *
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle
 * Modified by Chensong Zhang on 01/10/2012
 * Modified by Ludmil Zikatanov on 02/15/2013
 * Modified by Chensong Zhang on 05/10/2013: add a new input.
 */
void fasp_param_input (char *filenm, 
                       input_param *inparam)
{
    char     buffer[500]; // Note: max number of char for each line!
    int      val;
    SHORT    status = SUCCESS;
    
    // set default input parameters
    fasp_param_input_init(inparam);

    // if input file is not specified, use the default values
    if (filenm==NULL) return;
    
    FILE *fp = fopen(filenm,"r");
    if (fp==NULL) {
        printf("### ERROR: Could not open file %s...\n", filenm);
        fasp_chkerr(ERROR_OPEN_FILE, "fasp_param_input");

    }
    
    while ( status == SUCCESS ) {
        int     ibuff;
        double  dbuff;
        char    sbuff[500];
    
        val = fscanf(fp,"%s",buffer);
        if (val==EOF) break;
        if (val!=1){ status = ERROR_INPUT_PAR; break; }
        if (buffer[0]=='[' || buffer[0]=='%' || buffer[0]=='|') {
            fgets(buffer,500,fp); // skip rest of line
            continue;
        }
    
        // match keyword and scan for value
        if (strcmp(buffer,"workdir")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",sbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            strncpy(inparam->workdir,sbuff,128);
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"problem_num")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->problem_num=ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"print_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->print_level = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"output_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->output_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"solver_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->solver_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"stop_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->stop_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"precond_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->precond_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"itsolver_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->itsolver_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"itsolver_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->itsolver_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_ILU_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_ILU_levels = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_schwarz_levels")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = SUCCESS; break; }
			inparam->AMG_schwarz_levels = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
    
        else if (strcmp(buffer,"itsolver_restart")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->restart = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
    
            if ((strcmp(buffer,"C")==0)||(strcmp(buffer,"c")==0))
                inparam->AMG_type = CLASSIC_AMG;
            else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
                inparam->AMG_type = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
                inparam->AMG_type = UA_AMG;
            else
                { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_strong_coupled")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_strong_coupled = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_max_aggregation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_max_aggregation = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_tentative_smooth")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_tentative_smooth = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_smooth_filter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
    
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                inparam->AMG_smooth_filter = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                inparam->AMG_smooth_filter = OFF;
            else
                { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_coarse_scaling")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
    
            if ((strcmp(buffer,"ON")==0)||(strcmp(buffer,"on")==0)||
                (strcmp(buffer,"On")==0)||(strcmp(buffer,"oN")==0))
                inparam->AMG_coarse_scaling = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                inparam->AMG_coarse_scaling = OFF;
            else
                { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_levels = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_coarse_dof")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_coarse_dof = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_cycle_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
    
            if ((strcmp(buffer,"V")==0)||(strcmp(buffer,"v")==0))
                inparam->AMG_cycle_type = V_CYCLE;
            else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
                inparam->AMG_cycle_type = W_CYCLE;
            else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
                inparam->AMG_cycle_type = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
                inparam->AMG_cycle_type = NL_AMLI_CYCLE; 
            else
                { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_smoother")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
    
            if ((strcmp(buffer,"JACOBI")==0)||(strcmp(buffer,"jacobi")==0))
                inparam->AMG_smoother = SMOOTHER_JACOBI;
            else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
                inparam->AMG_smoother = SMOOTHER_GS;
            else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
                inparam->AMG_smoother = SMOOTHER_SGS;
            else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
                inparam->AMG_smoother = SMOOTHER_CG;    
            else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
                inparam->AMG_smoother = SMOOTHER_SOR;
            else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
                inparam->AMG_smoother = SMOOTHER_SSOR;
            else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
                inparam->AMG_smoother = SMOOTHER_GSOR;
            else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
                inparam->AMG_smoother = SMOOTHER_SGSOR;
            else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
                inparam->AMG_smoother = SMOOTHER_POLY;
            else if ((strcmp(buffer,"L1DIAG")==0)||(strcmp(buffer,"l1diag")==0))
                inparam->AMG_smoother = SMOOTHER_L1DIAG;
            else if ((strcmp(buffer,"BLKOIL")==0)||(strcmp(buffer,"blkoil")==0))
                inparam->AMG_smoother = SMOOTHER_BLKOIL;
            else
                { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_smooth_order")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%s",buffer);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            
            if ((strcmp(buffer,"NO")==0)||(strcmp(buffer,"no")==0))
                inparam->AMG_smooth_order = NO_ORDER;
            else if ((strcmp(buffer,"CF")==0)||(strcmp(buffer,"cf")==0))
                inparam->AMG_smooth_order = CF_ORDER;
            else
            { status = ERROR_INPUT_PAR; break; }
            fgets(buffer,500,fp); // skip rest of line
        }

        else if (strcmp(buffer,"AMG_coarsening_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_coarsening_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_interpolation_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_interpolation_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_aggressive_level = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_aggressive_path")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_aggressive_path = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_presmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_presmooth_iter = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_postsmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_postsmooth_iter = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_relaxation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_relaxation=dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_polynomial_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_polynomial_degree = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_strong_threshold")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_strong_threshold = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_truncation_threshold")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_truncation_threshold = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_max_row_sum")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_max_row_sum = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_amli_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_amli_degree = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_nl_amli_krylov_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->AMG_nl_amli_krylov_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->ILU_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_lfil")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->ILU_lfil = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_droptol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->ILU_droptol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_relax")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->ILU_relax = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_permtol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            inparam->ILU_permtol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"Schwarz_mmsize")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			inparam->Schwarz_mmsize = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"Schwarz_maxlvl")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) {status = ERROR_INPUT_PAR; break; }
			inparam->Schwarz_maxlvl = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
		
		else if (strcmp(buffer,"Schwarz_type")==0)
		{
			val = fscanf(fp,"%s",buffer);
			if (val!=1 || strcmp(buffer,"=")!=0) {
				status = ERROR_INPUT_PAR; break;
			}
			val = fscanf(fp,"%d",&ibuff);
			if (val!=1) { status = ERROR_INPUT_PAR; break; }
			inparam->Schwarz_type = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

        else {
            status = ERROR_INPUT_PAR;
            break;
        }    
    }
    
    fclose(fp);
    
    // sanity checks
    status = fasp_param_check(inparam);
    
#if DEBUG_MODE
    printf("### DEBUG: Reading input status = %d\n", status);
#endif
    
    // if meet unexpected input, stop the program
    fasp_chkerr(status,"fasp_param_input");
}

/**
 * \fn void fasp_param_set (int argc, const char *argv [], input_param *inparam)
 *
 * \brief read input from arguments
 *
 * \param argc       Number of arg input
 * \param argv       Input arguments
 * \param inparam    Parameters to be set
 *
 * \author Chensong Zhang
 * \date   12/29/2013
 */
void fasp_param_set (int argc,
                     const char *argv[],
                     input_param *inparam)
{
    int      arg_index   = 1;
    int      print_usage = 0;
    SHORT    status      = SUCCESS;

    // Option 1. set default input parameters
    fasp_param_input_init(inparam);
    
    while ( arg_index < argc ) {
        
        if ( strcmp(argv[arg_index], "-help") == 0 ) {
            print_usage = 1; break;
        }
        
        // Option 2. Get parameters from an ini file
        else if ( strcmp(argv[arg_index], "-ini") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Missing ini file name!\n");
                print_usage = 1; break;
            }
            strcpy(inparam->inifile, argv[arg_index]);
            fasp_param_input(inparam->inifile,inparam);
            if ( ++arg_index >= argc ) break;
        }
        
        // Option 3. Get parameters from command line input
        else if ( strcmp(argv[arg_index], "-print") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting print level (an integer between 0 and 10).\n");
                print_usage = 1; break;
            }
            inparam->print_level = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-output") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting output type (0 or 1).\n");
                print_usage = 1; break;
            }
            inparam->output_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-solver") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting solver type.\n");
                print_usage = 1; break;
            }
            inparam->solver_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-precond") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting preconditioner type.\n");
                print_usage = 1; break;
            }
            inparam->precond_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-maxit") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting max number of iterations.\n");
                print_usage = 1; break;
            }
            inparam->itsolver_maxit = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }
        
        else if ( strcmp(argv[arg_index], "-tol") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting tolerance for itsolver.\n");
                print_usage = 1; break;
            }
            inparam->itsolver_tol = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgmaxit") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting max number of iterations for AMG.\n");
                print_usage = 1; break;
            }
            inparam->AMG_maxit = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgtol") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting tolerance for AMG.\n");
                print_usage = 1; break;
            }
            inparam->AMG_tol = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgtype") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG type (1, 2, 3).\n");
                print_usage = 1; break;
            }
            inparam->AMG_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgcycle") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG cycle type (1, 2, 3).\n");
                print_usage = 1; break;
            }
            inparam->AMG_cycle_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgcoarsening") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG coarsening type.\n");
                print_usage = 1; break;
            }
            inparam->AMG_coarsening_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amginterplation") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG interpolation type.\n");
                print_usage = 1; break;
            }
            inparam->AMG_interpolation_type = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgsmoother") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG smoother type.\n");
                print_usage = 1; break;
            }
            inparam->AMG_smoother = atoi(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgsthreshold") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG strong threshold.\n");
                print_usage = 1; break;
            }
            inparam->AMG_strong_threshold = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else if ( strcmp(argv[arg_index], "-amgscouple") == 0 ) {
            arg_index++;
            if ( arg_index >= argc ) {
                printf("### ERROR: Expecting AMG strong coupled threshold.\n");
                print_usage = 1; break;
            }
            inparam->AMG_strong_coupled = atof(argv[arg_index]);
            if ( ++arg_index >= argc ) break;
        }

        else {
            print_usage = 1;
            break;
        }

    }
    
    if ( print_usage ) {
        
        printf("FASP command line options:\n");
        printf("================================================================\n");
        printf("  -ini              [CharValue] : Ini file name\n");
        printf("  -print            [IntValue]  : Print level\n");
        printf("  -output           [IntValue]  : Output to screen or a log file\n");
        printf("  -solver           [IntValue]  : Solver type\n");
        printf("  -precond          [IntValue]  : Preconditioner type\n");
        printf("  -maxit            [IntValue]  : Max number of iterations\n");
        printf("  -tol              [RealValue] : Tolerance for iterative solvers\n");
        printf("  -amgmaxit         [IntValue]  : Max number of AMG iterations\n");
        printf("  -amgtol           [RealValue] : Tolerance for AMG methods\n");
        printf("  -amgtype          [IntValue]  : AMG type\n");
        printf("  -amgcycle         [IntValue]  : AMG cycle type\n");
        printf("  -amgcoarsening    [IntValue]  : AMG coarsening type\n");
        printf("  -amginterpolation [IntValue]  : AMG interpolation type\n");
        printf("  -amgsmoother      [IntValue]  : AMG smoother type\n");
        printf("  -amgsthreshold    [RealValue] : AMG strong threshold\n");
        printf("  -amgscoupled      [RealValue] : AMG strong coupled threshold\n");
        printf("  -help                         : Brief help messages\n");

        exit(ERROR_INPUT_PAR);
        
    }
    
    // sanity checks
    status = fasp_param_check(inparam);
    
    // if meet unexpected input, stop the program
    fasp_chkerr(status,"fasp_param_set");
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
