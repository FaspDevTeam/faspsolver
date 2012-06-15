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
static void Diaginv(dCSRmat *Amat, REAL *Dinv);
static REAL DinvAnorminf(dCSRmat *Amat, REAL *Dinv);
static void Diagx(REAL *Dinv, INT n, REAL *x, REAL *b);
static void Rr(dCSRmat *Amat, REAL *Dinv, REAL *r, REAL *rbar, REAL *v0, REAL *v1, REAL *vnew, REAL *k, INT m);


/*---------------------------------*/
/*--      Public Function        --*/
/*---------------------------------*/

/**
 * \fn void fasp_smoother_dcsr_poly (dCSRmat *Amat, dvector *brhs, dvector *usol, INT n, INT ndeg, INT L)
 *
 * \brief poly approx to A^{-1} as MG smoother 
 *
 * \param Amat  Pointer to stiffness matrix, consider square matrix.
 * \param brhs  Pointer to right hand side
 * \param usol  Pointer to solution 
 * \param n     Problem size 
 * \param ndeg  Degree of poly 
 * \param L     Number of iterations
 *
 * \author Fei Cao, Xiaozhe Hu
 * \date   05/24/2012
 */
void fasp_smoother_dcsr_poly (dCSRmat *Amat, 
                              dvector *brhs, 
                              dvector *usol, 
                              INT n, 
                              INT ndeg, 
                              INT L)
{
    // local variables
    unsigned INT i;
    REAL *b = brhs->val, *u = usol->val;
    REAL *Dinv = NULL, *r = NULL, *rbar = NULL, *v0 = NULL, *v1 = NULL, *error = NULL, *k = NULL;
    REAL mu0, mu1, smu0, smu1;
    
    /* allocate memory */
    Dinv  = (REAL *) fasp_mem_calloc(n,sizeof(REAL)); 
    r     = (REAL *) fasp_mem_calloc(n,sizeof(REAL)); 
    rbar  = (REAL *) fasp_mem_calloc(n,sizeof(REAL)); 
    v0    = (REAL *) fasp_mem_calloc(n,sizeof(REAL)); 
    v1    = (REAL *) fasp_mem_calloc(n,sizeof(REAL)); 
    error = (REAL *) fasp_mem_calloc(n,sizeof(REAL));
    k     = (REAL *) fasp_mem_calloc(6,sizeof(REAL));// coefficients for calculation
    
    
    // get the inverse of the diagonal of A
    Diaginv(Amat, Dinv);
    
    // set up parameter
    mu0 = DinvAnorminf(Amat, Dinv); // get the inf norm of Dinv*A;
    
    
    mu0 = 1.0/mu0;
    mu1 = 4.0*mu0;// default set 8;
    smu0 =  sqrt(mu0); smu1 = sqrt(mu1);
    
    k[1] = (mu0+mu1)/2.0; 
    k[2] = (smu0 + smu1)*(smu0 + smu1)/2.0;
    k[3] = mu0 * mu1;
    k[4] = 2.0*k[3]/k[2];// 4.0*mu0*mu1/(sqrt(mu0)+sqrt(mu1))/(sqrt(mu0)+sqrt(mu1));
    k[5] = (mu1-2.0*smu0*smu1+mu0)/(mu1+2.0*smu0*smu1+mu0);// square of (sqrt(kappa)-1)/(sqrt(kappa)+1);
    
    
#if DEBUG_MODE
    printf("### DEBUG: fasp_smoother_dcsr_poly ...... [Start]\n");
#endif    
    
    // Update 
    for(i=0; i<L;i++ )
    {
        // get residual
        fasp_blas_dcsr_mxv(Amat, u, r);// r= Amat*u;
        fasp_blas_array_axpyz(n, -1, r, b, r);// r= -r+b;
        
        // Get correction error = R*r
        Rr(Amat, Dinv, r, rbar, v0, v1, error, k, ndeg);
        
        // update solution
        fasp_blas_array_axpy(n, 1, error, u);
        
    }
    
#if DEBUG_MODE
    printf("### DEBUG: the degrees of polysmoothing are: %d\n",ndeg);
#endif   
    
    // free memory
    fasp_mem_free(Dinv);
    fasp_mem_free(r);
    fasp_mem_free(rbar);
    fasp_mem_free(v0);
    fasp_mem_free(v1);
    fasp_mem_free(error);
    fasp_mem_free(k);
    
    return;
    
}

/**
 * \fn void fasp_smoother_dcsr_poly_old (dCSRmat *Amat, dvector *brhs, 
 *                                       dvector *usol, INT n, INT ndeg, INT L)
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
void fasp_smoother_dcsr_poly_old (dCSRmat *Amat, 
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
    //for (i=0; i < n; ++i) {res0 += r[i]*r[i];}
    //res0=sqrt(res0);
    
    for (it = 0 ; it < L; it++) { 
        bminax(b,ia,ja,a,u,&n,r);
        for (i=0; i < n ; ++i) {
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
        for (i=1; i < ndeg0; ++i) {
            for (j=0; j < n ; ++j) vsave[j]=v[j];
            for (j=0; j < n ; ++j) {
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
                //    u = u - a[jk] * x[j];
                u -= a[jk] * x[j];
            }
        res[i] = u;
    }
    return; 
}

/**
 * \fn    static void Diaginv(dCSRmat *Amat, REAL *Dinv)
 *
 * \brief find the inverse of the diagonal of A
 * 
 * \param Amat sparse matrix in CSR format
 * \param Dinv vector to store the inverse of diagonal of A
 *
 * \author Fei Cao, Xiaozhe Hu
 * \date 05/24/2012
 */
static void Diaginv(dCSRmat *Amat, REAL *Dinv)
{
    const INT   n  = Amat->row;
    const INT  *ia = Amat->IA, *ja = Amat->JA;
    const REAL *a  = Amat->val;
    unsigned INT i,j;
    
    for (i=0; i<n; i++)
    {
        for(j=ia[i]; j<ia[i+1]; j++)
        {
            if(i==ja[j]) // find the diagonal 
                break;
        }
        Dinv[i] = 1.0/a[j];
    }
    return;
}

/**
 * \fn    static REAL DinvAnorminf(dCSRmat *Amat, REAL *Dinv)
 *
 * \brief Get the inf norm of Dinv*A
 * 
 * \param Amat sparse matrix in CSR format
 * \param Dinv vector to store the inverse of diagonal of A
 *
 * \return Inf norm of Dinv*A
 *
 * \author Fei Cao, Xiaozhe Hu
 * \date   05/24/2012
 */
static REAL DinvAnorminf(dCSRmat *Amat, REAL *Dinv)
{
    //local variable
    const INT   n  = Amat->row;
    const INT  *ia = Amat->IA, *ja = Amat->JA;
    const REAL *a  = Amat->val;
    unsigned INT i,j;
    REAL norm, temp;
    
    // get the infinity norm of Dinv*A
    for (norm=0,i=0; i<n; i++) {

        for (temp=0,j=ia[i]; j<ia[i+1]; j++) {
            temp += ABS(a[ja[j]]); // Missing ja --Chensong 06/14/2012
        }
        
        temp *= Dinv[i]; // temp is the L1 norm of the ith row of Dinv*A;
        norm  = MAX(norm, temp);
    }
    
    return norm;
}

/**
 * \fn    static void Diagx(REAL *Dinv, INT n, REAL *x, REAL *b)
 *
 * \brief Compute b = Dinv * x;
 * 
 * \param Dinv vector to represent Diagonal matrix
 * \param n   problem size
 * \param b   vector
 * \param x   vector 
 *
 * \author Fei Cao, Xiaozhe Hu
 * \date 05/24/2012
 */
static void Diagx(REAL *Dinv, INT n, REAL *x, REAL *b)
{
    unsigned INT i;
    for (i=0; i<n; i++)
    {
        b[i] = Dinv[i] * x[i];
    }
    return;
}

/**
 * \fn static void Rr(dCSRmat *Amat, REAL *Dinv, REAL *r, REAL *rbar, REAL *v0, REAL *v1, REAL *vnew, REAL *k, INT m)
 *
 * \brief Compute action R*r, where R =  q_m(Dinv*A)*Dinv;
 * 
 * \param Amat sparse matrix A in CSR format
 * \param Dinv inverse of the diagoal of A
 * \param r  initial residual
 * \param rbar intermediate memory
 * \param v0 0 degree action R*r
 * \param v1 1 degree action R*r
 * \param vnew store final action vnew = R*r;
 * \param k  coefficients in equations
 * \param m  degree of polynomial
 *
 * \author Fei Cao, Xiaozhe Hu
 * \date 05/24/2012
 */
static void Rr(dCSRmat *Amat, REAL *Dinv, REAL *r, REAL *rbar, REAL *v0, REAL *v1, REAL *vnew, REAL *k, INT m)
{
    // local variables
    const INT   n  = Amat->row;
    INT i,j;
    
    
    //1 set up rbar
    Diagx(Dinv, n, r, rbar);// rbar = Dinv *r;
    
    //2 set up v0, v1;
    fasp_blas_dcsr_mxv(Amat, rbar, v1);//v1= A*rbar;
    Diagx(Dinv, n, v1, v1); // v1=Dinv *v1;
    
    for(i=0;i<n;i++)
    {
        v0[i] = k[1] * rbar[i];
        v1[i] = k[2] * rbar[i] - k[3] * v1[i];
    }
    
    //3 iterate to get v_(j+1)
    
    for (j=1;j<m;j++)
    {
        fasp_blas_dcsr_mxv(Amat, v1, rbar);//rbar= A*v_(j);
        
        for(i=0;i<n;i++)
        {
            rbar[i] = (r[i] - rbar[i])*Dinv[i];// indeed rbar=Dinv*(r-A*v_(j));
            vnew[i] = v1[i] + k[5] *(v1[i] - v0[i]) + k[4] * rbar[i];// compute v_(j+1)
            
            // prepare for next cycle
            v0[i]=v1[i]; 
            v1[i]=vnew[i];
            
        }
    }
    
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/