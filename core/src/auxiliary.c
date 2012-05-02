/*
 *  auxiliary.c
 *  
 *------------------------------------------------------
 * C-version: by Ludmil Zikatanov 2010-04-08
 *                         tested 2010-04-08
 *------------------------------------------------------
 *
 */

// Todo: add doxygen for each function! remove unwanted functions. --Chensong

/*! \file auxiliary.c
 *  \brief Read write and outher auxiliary routines. Also routine to convert 
 *  \brief ij value format to crs.
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

void fasp_aux_sh00 (dCSRmat *X, 
                    const INT ish)
{
    INT n=X->row, nnzX=X->nnz;
    unsigned INT i;
    // shift back
    for (i=0;i<=n;i++)   { X->IA[i] += ish; }
    for (i=0;i<nnzX;i++) { X->JA[i] += ish; }
    return;
}    

void fasp_aux_ijvcrs (INT *nnzi, 
                      INT *ia, 
                      INT *ja, 
                      REAL *a, 
                      INT *n,
                      INT *nnz,
                      INT *irow, 
                      INT *jcol, 
                      REAL *aval)
{
    INT nnzo=0,ni=0,nj=0,k=0,nzk,n1,irk=0,ica=0,icb=0,iend=0,jp=0;
    /*   ---------------------------------------------------------------*/
    /*  This subroutine converts the structures of a given matrix from*/
    /*  IJV=(IROW,ICOL,AVAL) to CRS=(IA,JA,A).*/
    /**/
    /*  Input:*/
    /*    nnzi number of nonzeroes in the input A */
    /*    IROW, JCOL - the row and column indicies, each o fsize nnzi */
    /*    AVAL(k)=A(I,J), I=IROW(K), J=JCOL(K); k=1:nnzi*/
    /**/
    /*  Output: */
    /*    N      - number of rows in A*/
    /*    NNZ    - number of nonzeros of A, should be equal tp nnzi*/
    /*    IA, JA - structure of A. */
    /*    A      - numerical values of nonzeros of A.*/
    /*-----------------------------------------------------------------*/
    fprintf(stderr,"\n\n");
    nnzo = *nnzi;
    for(k = 0; k < nnzo; k++){
        if(irow[k]>ni) ni =irow[k];
        if(jcol[k]>nj) nj=jcol[k];
        //    fprintf(stdout,"%d %d %d %d\n",ni,irow[k],nj,jcol[k]);
    }
    /* no need of such warning */
    /*    if (ni != nj) 
          fprintf(stderr,    \
          "**Warning: num of rows and columns do not match:%i .neq. %i\n",    \
          ni, nj);
    */
    n1=ni+1;
    for (k = 0; k< n1; k++) {
        ia[k] = 0;
    }
    for (k = 0; k< nnzo; k++) {
        /*  irk = irow[k] + 1; */
        irk = irow[k];
        ia[irk] = ia[irk] + 1;
    }
    nzk = 1;
    for (k = 1; k< n1; k++) {
        if (ia[k] != 0) {
            nzk = k;
            break;
        }
    }
    for (k = 0; k< nzk; k++) 
        ia[k] = 1;
    for (k = nzk; k< n1; k++) {
        ia[k] = ia[k-1] + ia[k];
    }
    //  fprintf(stdout,"========\n");
    for (k = 0; k < ni; k++) {
        ica = ia[k]-1;
        icb = ia[k+1]-1;
        //    fprintf(stdout,"row: from %d  to %d inclusive\n",ica,icb-1);
        if (icb > ica) ja[icb-1] = ica;
        /*      ja(icb-1) = ia(k) */
    }
    //  fprintf(stdout,"========\n");
    //  for (k = 0; k<=ni;k++){
    //    fprintf(stdout,"ia[%d]=%d\n",k+1,ia[k]);
    //  }
    //  fprintf(stdout,"========\n");
    //  fprintf(stdout,"---------------\n");
    for (k = 0; k< nnzo;k++){
        irk     = irow[k]-1;
        iend   = ia[irk+1]-2;
        jp     = ja[iend];
        ja[jp] = jcol[k];
        a[jp]  = aval[k];
        //    fprintf(stdout,"%d  %d\n",irk,iend);
        if (iend != jp)  
            ja[iend] = ja[iend] + 1;
    }
    //  fprintf(stdout,"----------------\n");
    *nnz=nnzo;
    *n=ni;
    return;
}

void fasp_aux_uuplv0_ (REAL *u,
                       REAL *v, 
                       INT *n)
{
    /*  
        This computes y = y + x.
    */
    INT i;
    for (i=0; i < *n ; i++) {
        u[i]=u[i]+v[i];
    }
    return; 
}

void fasp_aux_rveci (FILE *inp, 
                     INT *vec, 
                     INT *nn)
/* reads a vector of integers of size *nn from a file inp*/
/* the file "inp" should be open for reading */
{
    
    INT n;
    INT *vec_end;
    INT wall;
    n = *nn;
    vec_end  =  vec + n;
    for ( ; vec < vec_end; ++vec)
        wall = fscanf(inp,"%i",vec);
    fprintf(stdout,"Read %d INTEGERS", n);
    return;
}

void fasp_aux_rvecd (FILE *inp,
                     REAL *vec, 
                     INT *nn)
/* reads a vector of REALS of size nn from a file inp*/
{
    INT n;
    REAL *vec_end;
    INT wall;
    n=*nn;
    vec_end =  vec + n;
    for ( ; vec < vec_end; ++vec)
        wall = fscanf(inp,"%lg",vec);
    fprintf(stdout,"Read %d REALS", n);
    return;
}

void fasp_aux_wveci (FILE *inp, 
                     INT *vec, 
                     INT *nn)
/* writes a vector of integers of size nn from a file inp*/
{
    
    INT n;
    INT *vec_end;
    n = *nn;
    vec_end  =  vec + n;
    for ( ; vec < vec_end; ++vec)
        fprintf(inp,"%d ",*vec);
    fprintf(inp,"\n");
    fprintf(stdout,"Wrote %d INTEGERS", n);
    return;
}

void fasp_aux_wvecd (FILE *inp,
                     REAL *vec, 
                     INT *nn)
/* writes a vector of REALS of size nn from a file inp*/
{
    INT n;
    REAL *vec_end;  
    n=*nn;
    vec_end =  vec + n;
    for ( ; vec < vec_end; ++vec)
        fprintf(inp,"%24.16lg",*vec);
    fprintf(inp,"\n");
    fprintf(stdout,"Wrote %d REALS", n);
    return;
}

void fasp_aux_norm1_ (INT *ia,
                      INT *ja, 
                      REAL *a, 
                      INT *nn, 
                      REAL *a1norm)
{
    INT n,i,jk,iaa,iab;
    REAL sum,s;
    /* computes one norm of a matrix a and stores it in the variable
       pointed to by *a1norm*/
    n=*nn;
    s = 0e+00;
    for (i=0; i < n ; i++) {
        //fprintf(stdout,"\nia[%d] = %d",i,ia[i]);
        iaa = ia[i];
        iab = ia[i+1];
        sum = 0e+00;
        for (jk = iaa; jk < iab; jk++) {
            //   fprintf(stdout,"\nja[%d] = %d",jk,ja[jk]);
            //   fprintf(stdout,"a[%d] = %e",jk,a[jk]);
            sum += fabs(a[jk]);
        }
        if ( sum > s) s = sum; 
    }
    *a1norm=s;
}

void fasp_aux_auv_ (INT *ia,
                    INT *ja, 
                    REAL *a, 
                    REAL *u, 
                    REAL *v,
                    INT *nn, 
                    REAL *aauv)
{
    /* Calculation a(u,v)=(Au,v) */
    INT n,i,j,ij,iaa,iab;
    REAL sum,s;
    n=*nn;
    s = 0e+00;
    for (i=0; i < n ; i++) {
        iaa = ia[i];
        iab = ia[i+1];
        sum = 0e+00;
        for (ij = iaa; ij < iab; ij++) {
            j=ja[ij]-1;
            sum += a[ij]*u[j];
        }
        s=s+v[i]*sum;
    }
    *aauv=s;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
