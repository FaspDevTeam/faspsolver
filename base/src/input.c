/*! \file input.c
 *  \brief Read input parameters
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_param_input (char *filenm, input_param *Input)
 *
 * \brief Read input parameters from disk file
 *
 * \param filenm    File name for input file
 * \param Input     Input parameters
 *
 * \author Chensong Zhang
 * \date   03/20/2010 
 *
 * Modified by Xiaozhe Hu on 01/23/2011: add AMLI cycle
 * Modified by Chensong Zhang on 01/10/2012
 * Modified by Ludmil Zikatanov on 02/15/2013
 */
void fasp_param_input (char *filenm, 
                       input_param *Input)
{
    char     buffer[500]; // Note: max number of char for each line!
    int      val;
    SHORT    status = SUCCESS;
    
    // set default input parameters
    fasp_param_input_init(Input);

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
            strncpy(Input->workdir,sbuff,128);
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"problem_num")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->problem_num=ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"print_level")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->print_level = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"output_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->output_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"solver_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->solver_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"stop_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->stop_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"precond_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->precond_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"itsolver_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->itsolver_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"itsolver_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->itsolver_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_ILU_levels")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_ILU_levels = ibuff;
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
			Input->AMG_schwarz_levels = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}
    
        else if (strcmp(buffer,"itsolver_restart")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->restart = ibuff;
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
                Input->AMG_type = CLASSIC_AMG;
            else if ((strcmp(buffer,"SA")==0)||(strcmp(buffer,"sa")==0))
                Input->AMG_type = SA_AMG;
            else if ((strcmp(buffer,"UA")==0)||(strcmp(buffer,"ua")==0))
                Input->AMG_type = UA_AMG;
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
            Input->AMG_strong_coupled = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_max_aggregation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_max_aggregation = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_tentative_smooth")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_tentative_smooth = dbuff;
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
                Input->AMG_smooth_filter = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                Input->AMG_smooth_filter = OFF;
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
                Input->AMG_coarse_scaling = ON;
            else if ((strcmp(buffer,"OFF")==0)||(strcmp(buffer,"off")==0)||
                     (strcmp(buffer,"ofF")==0)||(strcmp(buffer,"oFf")==0)||
                     (strcmp(buffer,"Off")==0)||(strcmp(buffer,"oFF")==0)||
                     (strcmp(buffer,"OfF")==0)||(strcmp(buffer,"OFf")==0))
                Input->AMG_coarse_scaling = OFF;
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
            Input->AMG_levels = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_tol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_tol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_maxit")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_maxit = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_coarse_dof")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_coarse_dof = ibuff;
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
                Input->AMG_cycle_type = V_CYCLE;
            else if ((strcmp(buffer,"W")==0)||(strcmp(buffer,"w")==0))
                Input->AMG_cycle_type = W_CYCLE;
            else if ((strcmp(buffer,"A")==0)||(strcmp(buffer,"a")==0))
                Input->AMG_cycle_type = AMLI_CYCLE;
            else if ((strcmp(buffer,"NA")==0)||(strcmp(buffer,"na")==0))
                Input->AMG_cycle_type = NL_AMLI_CYCLE; 
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
                Input->AMG_smoother = SMOOTHER_JACOBI;
            else if ((strcmp(buffer,"GS")==0)||(strcmp(buffer,"gs")==0))
                Input->AMG_smoother = SMOOTHER_GS;
            else if ((strcmp(buffer,"SGS")==0)||(strcmp(buffer,"sgs")==0))
                Input->AMG_smoother = SMOOTHER_SGS;
            else if ((strcmp(buffer,"CG")==0)||(strcmp(buffer,"cg")==0))
                Input->AMG_smoother = SMOOTHER_CG;    
            else if ((strcmp(buffer,"SOR")==0)||(strcmp(buffer,"sor")==0))
                Input->AMG_smoother = SMOOTHER_SOR;
            else if ((strcmp(buffer,"SSOR")==0)||(strcmp(buffer,"ssor")==0))
                Input->AMG_smoother = SMOOTHER_SSOR;
            else if ((strcmp(buffer,"GSOR")==0)||(strcmp(buffer,"gsor")==0))
                Input->AMG_smoother = SMOOTHER_GSOR;
            else if ((strcmp(buffer,"SGSOR")==0)||(strcmp(buffer,"sgsor")==0))
                Input->AMG_smoother = SMOOTHER_SGSOR;
            else if ((strcmp(buffer,"POLY")==0)||(strcmp(buffer,"poly")==0))
                Input->AMG_smoother = SMOOTHER_POLY;
            else if ((strcmp(buffer,"L1DIAG")==0)||(strcmp(buffer,"l1diag")==0))
                Input->AMG_smoother = SMOOTHER_L1DIAG;
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
            Input->AMG_coarsening_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_interpolation_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_interpolation_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_presmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_presmooth_iter = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_postsmooth_iter")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_postsmooth_iter = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_relaxation")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_relaxation=dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_polynomial_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_polynomial_degree = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_strong_threshold")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_strong_threshold = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_truncation_threshold")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_truncation_threshold = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_max_row_sum")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_max_row_sum = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"AMG_amli_degree")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_amli_degree = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
        
        else if (strcmp(buffer,"AMG_nl_amli_krylov_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->AMG_nl_amli_krylov_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_type")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_type = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_lfil")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%d",&ibuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_lfil = ibuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_droptol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_droptol = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_relax")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_relax = dbuff;
            fgets(buffer,500,fp); // skip rest of line
        }
    
        else if (strcmp(buffer,"ILU_permtol")==0) {
            val = fscanf(fp,"%s",buffer);
            if (val!=1 || strcmp(buffer,"=")!=0) {
                status = ERROR_INPUT_PAR; break;
            }
            val = fscanf(fp,"%lf",&dbuff);
            if (val!=1) { status = ERROR_INPUT_PAR; break; }
            Input->ILU_permtol = dbuff;
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
			Input->Schwarz_mmsize = ibuff;
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
			Input->Schwarz_maxlvl = ibuff;
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
			Input->Schwarz_type = ibuff;
			fgets(buffer,500,fp); // skip rest of line
		}

        else {
            status = ERROR_INPUT_PAR;
            break;
        }    
    }
    
    fclose(fp);
    
    // sanity checks
    if ( Input->problem_num<0 
         || Input->print_level<0 
         || Input->solver_type<0
         || Input->precond_type<0 
         || Input->itsolver_tol<=0 
         || Input->itsolver_maxit<=0 
         || Input->ILU_type<=0
         || Input->ILU_lfil<0
         || Input->ILU_droptol<=0
         || Input->ILU_relax<0
         || Input->ILU_permtol<0
         || Input->AMG_type<=0 
         || Input->AMG_levels<0 
         || Input->AMG_ILU_levels<0
         || Input->AMG_tol<0 
         || Input->AMG_maxit<0 
         || Input->AMG_coarsening_type<0 
         || Input->AMG_smoother<0
         || Input->AMG_strong_threshold<0.0 
         || Input->AMG_strong_threshold>1.0 
         || Input->AMG_truncation_threshold<0.0 
         || Input->AMG_truncation_threshold>1.0 
         || Input->AMG_max_row_sum<0
         || Input->AMG_presmooth_iter<0 
         || Input->AMG_postsmooth_iter<0 
         || Input->AMG_strong_coupled<0 
         || Input->AMG_max_aggregation<=0 
         || Input->AMG_tentative_smooth<0 
         || Input->AMG_smooth_filter<0 
         || Input->stop_type<=0
         || Input->stop_type>3 
         || Input->restart<0
         ) status = ERROR_INPUT_PAR;
    
#if DEBUG_MODE
    printf("### DEBUG: Reading input status = %d\n", status);
#endif
    
    // if meet unexpected input, stop the program
    fasp_chkerr(status,"fasp_param_input");
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
