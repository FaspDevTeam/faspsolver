/*! \file rap.c
 *
 *  \brief R*A*P driver
 *
 *------------------------------------------------------
 * C-version by Ludmil Zikatanov 2010-04-08
 *                        tested 2010-04-08
 *------------------------------------------------------
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/*!
 * \fn dCSRmat fasp_blas_dcsr_rap2 (INT *ir, INT *jr, REAL *r, INT *ia, INT *ja, 
 *                                  REAL *a, INT *ipt, INT *jpt, REAL *pt, INT n, 
 *                                  INT nc, INT *maxrpout, INT *ipin, INT *jpin) 
 *
 * \brief Compute R*A*P 
 *
 * \author Ludmil Zikatanov
 * \date   04/08/2010
 *  
 * \note It uses dCSRmat only. The functions called from here are in sparse_util.c
 */
dCSRmat fasp_blas_dcsr_rap2 (INT *ir, 
                             INT *jr, 
                             REAL *r,
                             INT *ia, 
                             INT *ja, 
                             REAL *a,        
                             INT *ipt, 
                             INT *jpt, 
                             REAL *pt,    
                             INT n, 
                             INT nc,
                             INT *maxrpout,
                             INT *ipin, 
                             INT *jpin) 
{
    dCSRmat ac;
    INT n1,nc1,nnzp,maxrp;
    INT *ip=NULL,*jp=NULL;
    
    /*=========================================================*/
    /* 
       if ipin is null, this
       means that we need to do the transpose of p here; otherwise,
       these are considered to be input
    */
    maxrp=0;
    nnzp=ipt[nc]-1;
    n1=n+1;
    
    if(!ipin) {
        ip = (INT *)calloc(n1,sizeof(INT));
        jp = (INT *)calloc(nnzp,sizeof(INT));
        /* these must be null anyway, so no need to assign null
           ipin=NULL;
           jpin=NULL;
        */
    } else {
        ip=ipin;
        jp=jpin;
    }
    
    fasp_sparse_iit_(ipt,jpt,&nc,&n,ip,jp);
    
    /* tripple matrix product: R * A * transpose(P^T)=R*A*P.*/
    /* A is square n by n*/
    /* Note: to compute R*A* P the input are R, A and P^T */ 
    /* we need to transpose now the structure of P, because the input is P^T */
    /* end of transpose of the boolean corresponding to P */
    /* ic are the addresses of the rows of the output */
    nc1 = nc+1;
    ac.IA = (INT *)calloc(nc1,sizeof(INT));
    
    /*
       First call is with jc=null so that we find the number of
       nonzeroes in the result 
    */
    ac.JA=NULL;
    fasp_sparse_rapms_(ir,jr,ia,ja,ip,jp,&n,&nc,ac.IA,ac.JA,&maxrp);
    ac.nnz = ac.IA[nc]-1;
    ac.JA = (INT *)calloc(ac.nnz,sizeof(INT));
    
    /*
       second call is to fill the column indexes array jc. 
    */
    fasp_sparse_rapms_(ir,jr,ia,ja,ip,jp,&n,&nc,ac.IA,ac.JA,&maxrp);
    if(!ipin){
        if(ip) free(ip);
        if(jp) free(jp);
    }
    ac.val = (REAL *)calloc(ac.nnz,sizeof(REAL));
    /* this is the compute with the entries */
    fasp_sparse_rapcmp_(ir,jr,r,ia,ja,a,ipt,jpt,pt,&n,&nc,ac.IA,ac.JA,ac.val,&maxrp);
    ac.row=nc;
    ac.col=nc;
    
    /*=========================================================*/
    *maxrpout=maxrp;
    
    return ac;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
