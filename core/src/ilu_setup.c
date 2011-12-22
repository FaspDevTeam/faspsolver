/*! \file ilu_setup.c
 *  \brief Interface between FASP and ILU packages: setup ILU.
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/* ilu.for */
#ifdef __cplusplus 
extern "C" {void ilus_(const int *n,double *a,int *ja,int *ia,int *lfil,const double *relax,int *nzlu,double *luval,int *ijlu,int *maxstr,int *icode);}
extern "C" {void iluk_(const int *n,double *a,int *ja,int *ia,int *lfil,double *alu,int *jlu,int *iwk,int *ierr,int *nzlu);}
extern "C" {void ilut_(const int *n,double *a,int *ja,int *ia,int *lfil,const double *droptol,double *alu,int *jlu,int *iwk,int *ierr,int *nz);}
extern "C" {void ilutp_(const int *n,double *a,int *ja,int *ia,int *lfil,const double *droptol,const double *permtol,const int *mbloc,double *alu,int *jlu,int *iwk,int *ierr,int *nz);}
#else
extern void ilus_(const int *n,double *a,int *ja,int *ia,int *lfil,const double *relax,int *nzlu,double *luval,int *ijlu,int *maxstr,int *icode);
extern void iluk_(const int *n,double *a,int *ja,int *ia,int *lfil,double *alu,int *jlu,int *iwk,int *ierr,int *nzlu);
extern void ilut_(const int *n,double *a,int *ja,int *ia,int *lfil,const double *droptol,double *alu,int *jlu,int *iwk,int *ierr,int *nz);
extern void ilutp_(const int *n,double *a,int *ja,int *ia,int *lfil,const double *droptol,const double *permtol,const int *mbloc,double *alu,int *jlu,int *iwk,int *ierr,int *nz);
#endif

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int fasp_ilu_dcsr_setup(dCSRmat *A, ILU_data *iludata, ILU_param *param)
 * \brief Get ILU decoposition of a CSR matrix A
 *
 * \param *A         pointer to CSR matrir of double type
 * \param *iludata   pointer to ILU_data
 * \param *param     pointer to ILU parameters
 *
 * \author Shiquan Zhang
 * \date 12/27/2009
 */
int fasp_ilu_dcsr_setup (dCSRmat *A, 
												 ILU_data *iludata, 
												 ILU_param *param)
{
#if FASP_USE_ILU
	const int type=param->ILU_type, print_level=param->print_level;
	const int m=A->row, n=A->col, nnz=A->nnz, mbloc=n;
	const double ILU_relax=param->ILU_relax, ILU_droptol=param->ILU_droptol;
	const double permtol=param->ILU_permtol;
	
	int lfil=param->ILU_lfil, lfilt=param->ILU_lfil;
	int ierr, iwk, nzlu, nwork, maxstr, *ijlu;
	double *luval;
	
	clock_t setup_start, setup_end;
	double setup_duration;
	int status = SUCCESS;
	
#if DEBUG_MODE
	printf("fasp_ilu_dcsr_setup ...... [Start]\n");
#endif
	
	if (param->print_level>8) printf("fasp_ilu_dcsr_setup: m=%d, n=%d, nnz=%d\n",m,n,nnz);
	
	setup_start=clock();
	
	// Expected amount of memory for ILU needed and allocate memory 
	switch (type) {
		case ILUk:
			if (lfil == 0) iwk=nnz+500;
			else iwk=(lfil+2)*nnz;
			break;
		case ILUt:
			iwk=2*nnz; 	// iwk is the maxim possible nnz for ILU	
			lfilt=floor(n*0.5)+1;
			break;
		case ILUtp:
			iwk=2*nnz; 	// iwk is the maxim possible nnz for ILU	
			lfilt=floor(n*0.5)+1;
			break;
		default: // ILUs
			iwk=(lfil+10)*nnz;	
	} 
	
	nwork  = 4*n;
	maxstr = 3*(iwk+nwork);
	
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_dcsr_setup: fill-in=%d, iwk=%d, nwork=%d\n", lfil, iwk, nwork);
#endif
	
	// setup ILU preconditioner
	iludata->row=iludata->col=n;	
	fasp_ilu_data_alloc(iwk, nwork, iludata);
	
#if CHMEM_MODE
	printf("fasp_ilu_dcsr_setup: memory usage after fasp_ilu_data_alloc: \n");
	fasp_mem_usage();
#endif
	
	// ILU decomposition
	ijlu=iludata->ijlu;
	luval=iludata->luval;
	
	switch (type) {
		case ILUk:
			iluk_(&n,A->val,A->JA,A->IA,&lfil,luval,ijlu,&iwk,&ierr,&nzlu);
			break;
		case ILUt:
			ilut_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,luval,ijlu,&iwk,&ierr,&nzlu);
			break;
		case ILUtp:
			ilutp_(&n,A->val,A->JA,A->IA,&lfilt,&ILU_droptol,&permtol, \
						 &mbloc,luval,ijlu,&iwk,&ierr,&nzlu);
			break;
		default: // ILUs
			ilus_(&n,A->val,A->JA,A->IA,&lfil,&ILU_relax,&nzlu,luval,ijlu,&maxstr,&ierr);
			break;	
	} 
	
	fasp_dcsr_shift(A, -1);
	
#if CHMEM_MODE
	printf("fasp_ilu_dcsr_setup: memory usage after ILU setup: \n");
	fasp_mem_usage();
#endif
	
	iludata->nzlu=nzlu;
	iludata->nwork=nwork;
	
#if DEBUG_MODE
	printf("[DEBUG] fasp_ilu_dcsr_setup: iwk=%d, nzlu=%d\n",iwk,nzlu);
#endif	
	
	if (ierr!=0) {
		printf("Error: ILU setup failed (ierr=%d)!\n", ierr);
		status = ERROR_SOLVER_ILUSETUP;
		goto FINISHED;
	}
	
	if (iwk<nzlu) {
		printf("Error: Need more memory for ILU %d!\n", iwk-nzlu);
		status = ERROR_SOLVER_ILUSETUP;
		goto FINISHED;
	}
	
	if (print_level>0) {
		setup_end=clock();
		setup_duration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
		
		switch (type) {
			case ILUk:
				printf("ILUk setup costs %f seconds.\n", setup_duration);	
				break;
			case ILUt:
				printf("ILUt setup costs %f seconds.\n", setup_duration);	
				break;
			case ILUtp:
				printf("ILUtp setup costs %f seconds.\n", setup_duration);	
				break;
			case ILUs:
				printf("ILUs setup costs %f seconds.\n", setup_duration);	
				break;	
		} 		
	}
	
FINISHED: 	
	
#if DEBUG_MODE
	printf("fasp_ilu_dcsr_setup ...... [Finish]\n");
#endif
	
	return status;
	
#else // WITH_ILU
	
	printf("Error: ILU is not enabled!!!\n");
	exit(ERROR_MISC);
	
#endif
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
