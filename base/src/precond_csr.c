/*! \file precond_csr.c
 *  \brief Preconditioners
 */

#include "fasp.h"
#include "fasp_functs.h"
#include "forts_ns.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn precond *fasp_precond_setup (SHORT precond_type, AMG_param *amgparam, 
 *                                  ILU_param *iluparam, dCSRmat *A)
 *
 * \brief Setup preconditioner interface for iterative methods
 *
 * \param precond_type   Preconditioner type
 * \param *amgparam      AMG parameters
 * \param *iluparam      ILU parameters
 * \param *A             Pointer to coefficient matrix
 *
 * \return               Pointer to preconditioner
 *
 * \author Feiteng Huang
 * \date   05/18/2009
 */
precond *fasp_precond_setup (SHORT precond_type, 
                             AMG_param *amgparam, 
                             ILU_param *iluparam, 
                             dCSRmat *A)
{
    precond           *pc = NULL;
    
    AMG_data         *mgl = NULL;
    precond_data  *pcdata = NULL;
    ILU_data         *ILU = NULL;
    dvector         *diag = NULL;

    INT           max_levels, nnz, m, n;
    
    switch (precond_type) {
            
    case PREC_AMG: // AMG preconditioner
            
        pc = (precond *)fasp_mem_calloc(1, sizeof(precond));
        max_levels = amgparam->max_levels;
        nnz=A->nnz, m=A->row, n=A->col;    
            
        // initialize A, b, x for mgl[0]    
        mgl=fasp_amg_data_create(max_levels);
        mgl[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
        mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n); 
            
        // setup preconditioner  
        switch (amgparam->AMG_type) {
        case SA_AMG: // Smoothed Aggregation AMG
            fasp_amg_setup_sa(mgl, amgparam); break;
        case UA_AMG: // Unsmoothed Aggregation AMG
            fasp_amg_setup_ua(mgl, amgparam); break;
        default: // Classical AMG
#ifdef _OPENMP
            fasp_amg_setup_rs_omp(mgl, amgparam); break;
#else
            fasp_amg_setup_rs(mgl, amgparam); break;
#endif
        }
            
        pcdata = (precond_data *)fasp_mem_calloc(1, sizeof(precond_data));
        fasp_param_amg_to_prec(pcdata, amgparam);
        pcdata->max_levels = mgl[0].num_levels;
        pcdata->mgl_data = mgl;
            
        pc->data = pcdata;
            
        switch (amgparam->cycle_type) {
        case AMLI_CYCLE: // AMLI cycle
            pc->fct = fasp_precond_amli; break;
        case NL_AMLI_CYCLE: // Nonlinear AMLI AMG
            pc->fct = fasp_precond_nl_amli; break;
        default: // V,W-Cycle AMG
            pc->fct = fasp_precond_amg; break;
        }
            
        break;
            
    case PREC_FMG: // FMG preconditioner
            
        pc = (precond *)fasp_mem_calloc(1, sizeof(precond));
        max_levels = amgparam->max_levels;
        nnz=A->nnz, m=A->row, n=A->col;    
            
        // initialize A, b, x for mgl[0]    
        mgl=fasp_amg_data_create(max_levels);
        mgl[0].A=fasp_dcsr_create(m,n,nnz); fasp_dcsr_cp(A,&mgl[0].A);
        mgl[0].b=fasp_dvec_create(n); mgl[0].x=fasp_dvec_create(n); 
            
        // setup preconditioner  
        switch (amgparam->AMG_type) {
        case SA_AMG: // Smoothed Aggregation AMG
            fasp_amg_setup_sa(mgl, amgparam); break;
        case UA_AMG: // Unsmoothed Aggregation AMG
            fasp_amg_setup_ua(mgl, amgparam); break;
        default: // Classical AMG
#ifdef _OPENMP
            fasp_amg_setup_rs_omp(mgl, amgparam); break;
#else
            fasp_amg_setup_rs(mgl, amgparam); break;
#endif
        }
            
        pcdata = (precond_data *)fasp_mem_calloc(1, sizeof(precond_data));
        fasp_param_amg_to_prec(pcdata, amgparam);
        pcdata->max_levels = mgl[0].num_levels;
        pcdata->mgl_data = mgl;
            
        pc->data = pcdata; pc->fct = fasp_precond_famg;
            
        break;
            
    case PREC_ILU: // ILU preconditioner
            
        pc = (precond *)fasp_mem_calloc(1, sizeof(precond));
        ILU = (ILU_data *)fasp_mem_calloc(1, sizeof(ILU_data));
        fasp_ilu_dcsr_setup(A, ILU, iluparam);
        pc->data = ILU;
        pc->fct = fasp_precond_ilu;
            
        break;
            
    case PREC_DIAG: // Diagonal preconditioner
            
        pc = (precond *)fasp_mem_calloc(1, sizeof(precond));
        diag = (dvector *)fasp_mem_calloc(1, sizeof(dvector));
        fasp_dcsr_getdiag(0, A, diag);    
            
        pc->data = diag; 
        pc->fct  = fasp_precond_diag;
            
        break;
            
    default: // No preconditioner
            
        break;
            
    }
    
    return pc;
}

/**
 * \fn void fasp_precond_diag (REAL *r, REAL *z, void *data)
 *
 * \brief Diagonal preconditioner z=inv(D)*r
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Chensong Zhang
 * \date   04/06/2010
 */
void fasp_precond_diag (REAL *r, 
                        REAL *z, 
                        void *data)
{
    dvector *diag=(dvector *)data;
    REAL *diagptr=diag->val;
    INT i, m=diag->row;    
    
    memcpy(z,r,m*sizeof(REAL));
    for (i=0;i<m;++i) {
        if (ABS(diag->val[i])>SMALLREAL) z[i]/=diagptr[i];
    }    
}

/**
 * \fn void fasp_precond_ilu (REAL *r, REAL *z, void *data)
 *
 * \brief ILU preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Shiquan Zhang
 * \date   04/06/2010
 */
void fasp_precond_ilu (REAL *r, 
                       REAL *z, 
                       void *data)
{
    ILU_data *iludata=(ILU_data *)data;
    const INT m=iludata->row, mm1=m-1, memneed=2*m;
    REAL *zz, *zr;
    
    if (iludata->nwork<memneed) goto MEMERR; // check this outside this subroutine!!
    
    zz = iludata->work; 
    zr = iludata->work+m;
    fasp_array_cp(m, r, zr);     
    
    {
        INT i, j, jj, begin_row, end_row, mm2=m-2;
        INT *ijlu=iludata->ijlu;
        REAL *lu=iludata->luval;
//	REAL tmp;
    
        // forward sweep: solve unit lower matrix equation L*zz=zr
        zz[0]=zr[0];
/*
 * #ifdef _OPENMP	
 * #pragma omp parallel for private(i,j,jj,begin_row,end_row) ordered //reduction(+:tmp) 
 * #endif
*/	
        for (i=1;i<=mm1;++i) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=begin_row;j<=end_row;++j) {
                jj=ijlu[j];
                if (jj<i) zr[i]-=lu[j]*zz[jj];
                else break;
            }
            zz[i]=zr[i]; 
        }
    
        // backward sweep: solve upper matrix equation U*z=zz
        z[mm1]=zz[mm1]*lu[mm1];
        for (i=mm2;i>=0;i--) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=end_row;j>=begin_row;j--) {
                jj=ijlu[j];
                if (jj>i) zz[i]-=lu[j]*z[jj];
                else break;
            } 
            z[i]=zz[i]*lu[i];
        }
    }
    
    return;
    
 MEMERR:
    printf("### ERROR: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
    exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_ilu_forward (REAL *r, REAL *z, void *data)
 *
 * \brief ILU preconditioner: only forwear sweep
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu, Shiquang Zhang
 * \date   04/06/2010
 */
void fasp_precond_ilu_forward (REAL *r, 
                               REAL *z, 
                               void *data)
{
    ILU_data *iludata=(ILU_data *)data;
    const INT m=iludata->row, mm1=m-1, memneed=2*m;
    REAL *zz, *zr;
    
    if (iludata->nwork<memneed) goto MEMERR; 
    
    zz = iludata->work; 
    zr = iludata->work+m;
    fasp_array_cp(m, r, zr);     
    
    {
        INT i, j, jj, begin_row, end_row;
        INT *ijlu=iludata->ijlu;
        REAL *lu=iludata->luval;
    
        // forward sweep: solve unit lower matrix equation L*z=r
        zz[0]=zr[0];
        for (i=1;i<=mm1;++i) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=begin_row;j<=end_row;++j) {
                jj=ijlu[j];
                if (jj<i) zr[i]-=lu[j]*zz[jj];
                else break;
            }
            zz[i]=zr[i];
        }
    }
    
    fasp_array_cp(m, zz, z); 
    
    return;
    
 MEMERR:
    printf("### ERROR: Need %d memory, only %d available!!!", memneed, iludata->nwork);
    exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_ilu_backward (REAL *r, REAL *z, void *data)
 *
 * \brief ILU preconditioner: only backward sweep
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu, Shiquan  Zhang
 * \date   04/06/2010
 */
void fasp_precond_ilu_backward (REAL *r, 
                                REAL *z, 
                                void *data)
{
    ILU_data *iludata=(ILU_data *)data;
    const INT m=iludata->row, mm1=m-1, memneed=2*m;
    REAL *zz;//, *zr;
    
    if (iludata->nwork<memneed) goto MEMERR; 
    
    zz = iludata->work; 
    //zr = iludata->work+m;
    fasp_array_cp(m, r, zz);     
    
    {
        INT i, j, jj, begin_row, end_row, mm2=m-2;
        INT *ijlu=iludata->ijlu;
        REAL *lu=iludata->luval;
    
        // backward sweep: solve upper matrix equation U*z=zz
        z[mm1]=zz[mm1]*lu[mm1];
        for (i=mm2;i>=0;i--) {
            begin_row=ijlu[i]; end_row=ijlu[i+1]-1;
            for (j=end_row;j>=begin_row;j--) {
                jj=ijlu[j];
                if (jj>i) zz[i]-=lu[j]*z[jj];
                else break;
            } 
            z[i]=zz[i]*lu[i];
        }
    
    }
    
    return;
    
 MEMERR:
    printf("### ERROR: Need %d memory, only %d available!!!", memneed, iludata->nwork);
    exit(ERROR_ALLOC_MEM);
}

/**
 * \fn void fasp_precond_schwarz(REAL *r, REAL *z, void *data)
 * \brief get z from r by schwarz
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date 03/22/2010
 */
void fasp_precond_schwarz(REAL *r, 
                          REAL *z, 
                          void *data)
{
	Schwarz_data *schwarz_data=(Schwarz_data *)data;
	
	INT n = schwarz_data->A.row;
	INT *ia = schwarz_data->A.IA;
	INT *ja = schwarz_data->A.JA;
	REAL *a = schwarz_data->A.val;
	
	INT nblk = schwarz_data->nblk;
	INT *iblock = schwarz_data->iblock;
	INT *jblock = schwarz_data->jblock;
	REAL *rhsloc = schwarz_data->rhsloc;
	REAL *au = schwarz_data->au;
	REAL *al = schwarz_data->al;
	
	INT schwarz_type = schwarz_data->schwarz_type;
	
	INT memt = schwarz_data->memt;
	INT *mask = schwarz_data->mask;
	INT *maxa = schwarz_data->maxa;
	
	fasp_array_set(n, z, 0.0);
	
	switch (schwarz_type)
	{
		case 2:
			bbgs2ns_(&n,ia,ja,a,z,r,&nblk,iblock,jblock,mask,maxa,au,al,rhsloc,&memt);
			break;
		case 3:
			fbgs2ns_(&n,ia,ja,a,z,r,&nblk,iblock,jblock,mask,maxa,au,al,rhsloc,&memt);
			bbgs2ns_(&n,ia,ja,a,z,r,&nblk,iblock,jblock,mask,maxa,au,al,rhsloc,&memt);
			break;
		default:
			fbgs2ns_(&n,ia,ja,a,z,r,&nblk,iblock,jblock,mask,maxa,au,al,rhsloc,&memt);
			break;
	}
    
}

/**
 * \fn void fasp_precond_amg (REAL *r, REAL *z, void *data)
 *
 * \brief AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Chensong Zhang
 * \date   04/06/2010
 */
void fasp_precond_amg (REAL *r, 
                       REAL *z, 
                       void *data)
{
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    INT i;
    
    AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
    mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) fasp_solver_mgcycle(mgl,&amgparam);

    // We can use a recurcive version of MG:
    //   for (i=0;i<maxit;++i) fasp_solver_mgrecur(mgl,&amgparam,0);
    
    fasp_array_cp(m,mgl->x.val,z);    
}

/**
 * \fn void fasp_precond_famg (REAL *r, REAL *z, void *data)
 *
 * \brief Full AMG perconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   02/27/2011
 */
void fasp_precond_famg (REAL *r, 
                        REAL *z, 
                        void *data)
{
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    INT i;
    
    AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
    mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) fasp_solver_fmgcycle(mgl,&amgparam);
    
    fasp_array_cp(m,mgl->x.val,z);    
}

/**
 * \fn void fasp_precond_amli(REAL *r, REAL *z, void *data)
 *
 * \brief AMLI AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   01/23/2011
 */
void fasp_precond_amli (REAL *r, 
                        REAL *z, 
                        void *data)
{
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    INT i;
    
    AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
    mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) fasp_solver_amli(mgl,&amgparam,0); 
    
    fasp_array_cp(m,mgl->x.val,z);    
}

/**
 * \fn void fasp_precond_nl_amli(REAL *r, REAL *z, void *data)
 *
 * \brief Nonliear AMLI AMG preconditioner
 *
 * \param r     Pointer to the vector needs preconditioning
 * \param z     Pointer to preconditioned vector
 * \param data  Pointer to precondition data
 *
 * \author Xiaozhe Hu
 * \date   04/25/2011
 */
void fasp_precond_nl_amli (REAL *r, 
                           REAL *z, 
                           void *data)
{
    
    precond_data *pcdata=(precond_data *)data;
    const INT m=pcdata->mgl_data[0].A.row;
    const INT maxit=pcdata->maxit;
    const SHORT num_levels = pcdata->max_levels;
    INT i;
    
    AMG_param amgparam; fasp_param_amg_init(&amgparam);
    fasp_param_prec_to_amg(&amgparam,pcdata);
    
    AMG_data *mgl = pcdata->mgl_data;
    mgl->b.row=m; fasp_array_cp(m,r,mgl->b.val); // residual is an input 
    mgl->x.row=m; fasp_dvec_set(m,&mgl->x,0.0);
    
    for (i=0;i<maxit;++i) fasp_solver_nl_amli(mgl,&amgparam,0, num_levels); 
    
    fasp_array_cp(m,mgl->x.val,z);    
}

/**
 * \fn void fasp_precond_free (SHORT precond_type, precond *pc) 
 *
 * \brief free preconditioner
 *
 * \param precond_type   Preconditioner type
 * \param *pc            precondition data & fct
 *
 * \return               void
 *
 * \author Feiteng Huang
 * \date   12/24/2012
 */
void fasp_precond_free (SHORT precond_type, precond *pc)
{
    switch (precond_type) {
            
    case PREC_AMG: // AMG preconditioner

        fasp_amg_data_free(((precond_data*)(pc->data))->mgl_data);
        fasp_mem_free((precond_data*)(pc->data));
        fasp_mem_free(pc);
            
        break;
            
    case PREC_FMG: // FMG preconditioner
            
        fasp_amg_data_free(((precond_data*)(pc->data))->mgl_data);
        fasp_mem_free((precond_data*)(pc->data));
        fasp_mem_free(pc);
            
        break;
            
    case PREC_ILU: // ILU preconditioner
            
        fasp_ilu_data_free((ILU_data*)(pc->data));
        fasp_mem_free(pc);
            
        break;
            
    case PREC_DIAG: // Diagonal preconditioner
            
        fasp_dvec_free((dvector*)(pc->data));
        fasp_mem_free(pc);
            
        break;
            
    default: // No preconditioner
            
        break;
            
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
