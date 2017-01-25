/*! \file AuxGivens.c
 *
 *  \brief Givens transformation
 *
 *  \note This file contains Level-0 (Aux) functions
 */

#include <math.h>

#include "fasp.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_aux_givens (const REAL beta, const dCSRmat *H, dvector *y, 
 *                           REAL *tmp)
 *
 * \brief Perform Givens rotations to compute y |beta*e_1- H*y|
 *
 * \param beta   Norm of residual r_0
 * \param H      Upper Hessenberg dCSRmat matrix: (m+1)*m
 * \param y      Minimizer of |beta*e_1- H*y|
 * \param tmp    Temporary work array
 *
 * \author Xuehai Huang
 * \date   10/19/2008
 */
void fasp_aux_givens (const REAL      beta,
                      const dCSRmat  *H,
                      dvector        *y,
                      REAL           *tmp)
{
    const INT  Hsize=H->row;
    INT        i, j, istart, idiag, ip1start;
    REAL       h0,h1,r,c,s,tempi,tempip1,sum;
    
    tmp[0]=beta;
    memset(&tmp[1], 0x0, sizeof(REAL)*(Hsize-1));

    for (i=0;i<Hsize-1;++i) {
        istart=H->IA[i];
        ip1start=H->IA[i+1];
        if (i==0) idiag=istart;
        else idiag=istart+1;
    
        h0=H->val[idiag]; // h0=H[i][i]
        h1=H->val[H->IA[i+1]]; // h1=H[i+1][i]
        r=sqrt(h0*h0+h1*h1);
        c=h0/r; s=h1/r;
    
        for (j=idiag;j<ip1start;++j) {
            tempi=H->val[j];
            tempip1=H->val[ip1start+(j-idiag)];
            H->val[j]=c*tempi+s*tempip1;
            H->val[ip1start+(j-idiag)]=c*tempip1-s*tempi;
        }
    
        tempi=c*tmp[i]+s*tmp[i+1];
        tempip1=c*tmp[i+1]-s*tmp[i];
    
        tmp[i]=tempi; tmp[i+1]=tempip1;
    }
    
    for (i=Hsize-2;i>=0;--i) {
        sum=tmp[i];
        istart=H->IA[i];
        if (i==0) idiag=istart;
        else idiag=istart+1;
    
        for (j=Hsize-2;j>i;--j) sum-=H->val[idiag+j-i]*y->val[j];
    
        y->val[i]=sum/H->val[idiag];
    }
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
