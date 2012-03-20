/*! \file smoother_poly.c
 *  \brief Smoothers for sparse matrix in CSR format using poly. approx. to A^{-1} 
 */

#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include "fasp.h"
#include "fasp_functs.h"

static void bminax(REAL *b,INT *ia,INT *ja, REAL *a, REAL *x,INT *nn, REAL *res);

/*---------------------------------*/
/*--      Public Function        --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dcsr_poly (dCSRmat *Amat, dvector *brhs, dvector *usol, INT n, INT ndeg, INT L)
 *
 * \brief poly approx to A^{-1} as MG smoother: JK&LTZ2010
 *
 * \param Amat  Pointer to stiffness matrix
 * \param brhs  Pointer to right hand side
 * \param usol  Pointer to solution 
 * \param n     Problem size 
 * \param ndeg  Degree of poly 
 * \param L     Number of iterations
 *
 * \author James Brannick and Ludmil T Zikatanov
 * \date   06/28/2010
 */
void fasp_smoother_dcsr_poly (dCSRmat *Amat, 
                              dvector *brhs, 
                              dvector *usol, 
                              INT n, 
                              INT ndeg, 
                              INT L)
{
    INT  *ia=Amat->IA,*ja=Amat->JA;
    REAL *a=Amat->val, *b=brhs->val, *u=usol->val;
	
    INT   i,j,k,it,jk,iaa,iab,ndeg0;  // id and ij for scaling of A
    REAL *v,*v0,*r,*vsave;  // one can get away without r as well;
    REAL  smaxa,smina,delinv,s,smu0,smu1,skappa,th,th1,sq; 
    REAL  ri,ari,vj,ravj,snj,sm,sm01,smsqrt,delta,delta2,chi;
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_smoother_dcsr_poly ...... [Start]\n");
#endif	
	
    /* WORKING MEM */
    v     = (REAL *) fasp_mem_calloc(n,sizeof(REAL));  
    v0    = (REAL *) fasp_mem_calloc(n,sizeof(REAL));  
    vsave = (REAL *) fasp_mem_calloc(n,sizeof(REAL));  
    r     = (REAL *) fasp_mem_calloc(n,sizeof(REAL));  
	
    /* COMPUTE PARAMS*/
    // min INT for approx -- could be done upfront 
    // i.e., only once per level... only norm1 ...
    fasp_aux_norm1_(ia,ja,a,&n,&smaxa);
    smina=smaxa/8;
    delinv=(smaxa+smina)/(smaxa-smina);
    th=delinv+sqrt(delinv*delinv-1e+00);
    th1=1e+00/th;
    sq=(th-th1)*(th-th1);
    //
    ndeg0=floor(log(2*(2e0+th+th1)/sq)/log(th)+1e0);
    if(ndeg0 < ndeg) ndeg0=ndeg;
    //
    smu0=1e+00/smaxa;
    smu1=1e+00/smina;
    skappa=sqrt(smaxa/smina);
    delta=(skappa-1e+00)/(skappa+1);
    delta2=delta*delta;
    s=sqrt(smu0)+sqrt(smu1);
    s=s*s;
    smsqrt=0.5e+00*s;
    chi=4e+00*smu0*smu1/s;
    sm=0.5e+00*(smu0+smu1);
    sm01=smu0*smu1;
	
#if DEBUG_MODE
    printf("### DEBUG: the degrees of polysmoothing are: %d %d\n",ndeg0,ndeg);
#endif		
	
    /* BEGIN POLY ITS */
	
    /* auv_(ia,ja,a,u,u,&n,&err0); NA: u = 0 */
    //bminax(b,ia,ja,a,u,&n,r);
    //for (i=0; i < n ; ++i){res0 += r[i]*r[i];}
    //res0=sqrt(res0);
    
    for (it = 0 ; it < L; it++){ 
        bminax(b,ia,ja,a,u,&n,r);
        for (i=0; i < n ; ++i){
            iaa = ia[i];
            iab = ia[i+1];
            ari=0e+00; /* ari is (A*r)[i] */ 
            if(iab > iaa) {
                for (jk = iaa; jk < iab; jk++) {
					j=ja[jk];
					ari += a[jk] * r[j];
                } 
            }
            ri=r[i];
            v0[i]=sm*ri; 
            v[i]=smsqrt*ri-sm01*ari; 
        }
        for (i=1; i < ndeg0; ++i){
            for (j=0; j < n ; ++j) vsave[j]=v[j];
            for (j=0; j < n ; ++j){
                /* ravj = (r- A*v)[j] */
                ravj= r[j];
                iaa = ia[j];
                iab = ia[j+1];
                if(iab > iaa) {
                    for (jk = iaa; jk < iab; jk++) {
						k=ja[jk];
						ravj -= a[jk] * vsave[k];
					}
                }
                vj=v[j];
                snj = chi*ravj+delta2*(vj-v0[j]);
                v0[j]=vj;   
                v[j]=vj+snj;
            }
        }
        fasp_aux_uuplv0_(u,v,&n);
        //bminax(b,ia,ja,a,u,&n,r);
        //for (i=0; i < n ; ++i)
        //resk += r[i]*r[i];
        //resk=sqrt(resk);
        //fprintf(stdout,"\nres0=%12.5g\n",res0);
        //fprintf(stdout,"\nresk=%12.5g\n",resk);
        //res0=resk;
        //resk=0.0e0;
    }
	
    fasp_mem_free(v);
    fasp_mem_free(v0);
    fasp_mem_free(r);
    fasp_mem_free(vsave);
	
#if DEBUG_MODE
    printf("### DEBUG: fasp_smoother_dcsr_poly ...... [FINISH]\n");
#endif	
	
    return; 
}

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void bminax(REAL *b,INT *ia,INT *ja, REAL *a, REAL *x,INT *nn, REAL *res)
 *
 * \brief ???
 *
 * \param b ???
 * \param ia ???
 * \param ja ???
 * \param a ???
 * \param x ???
 * \param nn ???
 * \param res ???
 *
 * \author James Brannick and Ludmil T Zikatanov
 * \date 06/28/2010
 */
static void bminax(REAL *b,INT *ia,INT *ja, REAL *a, REAL *x,INT *nn, REAL *res)
{
    /* Computes b-A*x */
	
    INT i,j,jk,iaa,iab;
    INT n;
    REAL u;
    n=*nn;
    for (i=0; i < n ; ++i) {
        iaa = ia[i];
        iab = ia[i+1];
        u = b[i];
        if(iab > iaa) 
            for (jk = iaa; jk < iab; jk++) {
				j=ja[jk];
				//	u = u - a[jk] * x[j];
				u -= a[jk] * x[j];
            }
        res[i] = u;
    }
    return; 
}

/*
 sqdinv= (REAL *)fasp_mem_calloc(n,sizeof(REAL));
 
 // find the diagonals: note f numbering assumption
 for (i=0; i < n; ++i) {
 iaa=ia[i]-1;
 iab=ia[i+1]-1;
 id=-1;
 for(ij=iaa; ij<iab; ++ij){
 j=ja[ij]-1;
 id=ij;
 if(j == i) break;
 }
 sqdinv[i]=1e+00/sqrt(a[id]);
 }
 for (i=0; i < n; ++i) {
 iaa=ia[i]-1;
 iab=ia[i+1]-1;
 for(ij=iaa; ij<iab; ++ij){
 j=ja[ij]-1;
 a[ij]=sqdinv[i]*sqdinv[j]*a[ij];
 //fprintf(stdout,"\na[%d,%d]=%25.16lg",i+1,j+1,a[ij]);
 }
 }
 */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
