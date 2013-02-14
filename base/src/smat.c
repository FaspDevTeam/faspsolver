/*! \file smat.c
 *  \brief Simple operations for *small* full matrices in row-major format
 *
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_blas_smat_inv_nc2 (REAL *a)
 *
 * \brief Compute the inverse matrix of a 2*2 full matrix A (in place)
 *
 * \param a   Pointer to the REAL array which stands a 2*2 matrix 
 * 
 * \author Xiaozhe Hu
 * \date   18/11/2011
 */
void fasp_blas_smat_inv_nc2 (REAL *a)
{ 
    const REAL a0 = a[0], a1 = a[1];
    const REAL a2 = a[2], a3 = a[3];
    
    const REAL det = a0*a3 - a1*a2;
    REAL det_inv;
    
    if (ABS(det)<1e-22) {
        printf("### WARNING: Matrix is nearly singular! det = %e\n", det);
        /*
         printf("##----------------------------------------------\n");
         printf("## %12.5e %12.5e \n", a0, a1);
         printf("## %12.5e %12.5e \n", a2, a3);
         printf("##----------------------------------------------\n");
         getchar(); 
         */
    }
    
    det_inv = 1.0/det;
    
    a[0] = a3*det_inv; a[1] = - a1*det_inv; 
    
    a[2] = - a2*det_inv; a[3] = a0*det_inv;     
}

/**
 * \fn void fasp_blas_smat_inv_nc3 (REAL *a)
 *
 * \brief Compute the inverse matrix of a 3*3 full matrix A (in place)
 *
 * \param a  Pointer to the REAL array which stands a 3*3 matrix 
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   05/01/2010
 */
void fasp_blas_smat_inv_nc3 (REAL *a)
{ 
    const REAL a0 = a[0], a1 = a[1], a2 = a[2];
    const REAL a3 = a[3], a4 = a[4], a5 = a[5];
    const REAL a6 = a[6], a7 = a[7], a8 = a[8];
    
    const REAL M0 = a4*a8-a5*a7, M3 = a2*a7-a1*a8, M6 = a1*a5-a2*a4;
    const REAL M1 = a5*a6-a3*a8, M4 = a0*a8-a2*a6, M7 = a2*a3-a0*a5;
    const REAL M2 = a3*a7-a4*a6, M5 = a1*a6-a0*a7, M8 = a0*a4-a1*a3;
    const REAL det = a0*M0+a3*M3+a6*M6;
    REAL det_inv;
    
    if (ABS(det)<1e-22) {
        printf("### WARNING: Matrix is nearly singular! det = %e\n", det);
        /*
         printf("##----------------------------------------------\n");
         printf("## %12.5e %12.5e %12.5e \n", a0, a1, a2);
         printf("## %12.5e %12.5e %12.5e \n", a3, a4, a5);
         printf("## %12.5e %12.5e %12.5e \n", a5, a6, a7);
         printf("##----------------------------------------------\n");
         getchar(); 
         */
    }
    
    det_inv = 1.0/det;
    
    a[0] = M0*det_inv; a[1] = M3*det_inv; a[2] = M6*det_inv;
    
    a[3] = M1*det_inv; a[4] = M4*det_inv; a[5] = M7*det_inv;
    
    a[6] = M2*det_inv; a[7] = M5*det_inv; a[8] = M8*det_inv;    
}

/**
 * \fn void fasp_blas_smat_inv_nc4 (REAL *a)
 *
 * \brief Compute the inverse matrix of a 4*4 full matrix A (in place)
 *
 * \param a  Pointer to the REAL array which stands a 4*4 matrix
 *
 * \author Xiaozhe Hu
 * \date   01/12/2013
 */
void fasp_blas_smat_inv_nc4 (REAL *a)
{
    const REAL a11 = a[0],  a12 = a[1],  a13 = a[2],  a14 = a[3];
    const REAL a21 = a[4],  a22 = a[5],  a23 = a[6],  a24 = a[7];
    const REAL a31 = a[8],  a32 = a[9],  a33 = a[10], a34 = a[11];
    const REAL a41 = a[12], a42 = a[13], a43 = a[14], a44 = a[15];

    const REAL M11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
    const REAL M12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
    const REAL M13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
    const REAL M14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
    const REAL M21 = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
    const REAL M22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;
    const REAL M23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a41 - a13*a24*a41 - a14*a21*a43;
    const REAL M24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;
    const REAL M31 = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
    const REAL M32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
    const REAL M33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
    const REAL M34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;
    const REAL M41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
    const REAL M42 = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
    const REAL M43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
    const REAL M44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;
    
    const REAL det = a11*M11 + a12*M21 + a13*M31 + a14*M41;
    REAL det_inv;
    
    if (ABS(det)<1e-22) {
        printf("### WARNING: Matrix is nearly singular! det = %e\n", det);
        /*
         printf("##----------------------------------------------\n");
         printf("## %12.5e %12.5e %12.5e \n", a0, a1, a2);
         printf("## %12.5e %12.5e %12.5e \n", a3, a4, a5);
         printf("## %12.5e %12.5e %12.5e \n", a5, a6, a7);
         printf("##----------------------------------------------\n");
         getchar();
         */
    }
    
    det_inv = 1.0/det;
    
    a[0] = M11*det_inv;  a[1] = M12*det_inv;  a[2] = M13*det_inv;  a[3] = M14*det_inv;
    a[4] = M21*det_inv;  a[5] = M22*det_inv;  a[6] = M23*det_inv;  a[7] = M24*det_inv;
    a[8] = M31*det_inv;  a[9] = M32*det_inv;  a[10] = M33*det_inv; a[11] = M34*det_inv;
    a[12] = M41*det_inv; a[13] = M42*det_inv; a[14] = M43*det_inv; a[15] = M44*det_inv;
    
}

/**
 * \fn void fasp_blas_smat_inv_nc5 (REAL *a)
 *
 * \brief Compute the inverse matrix of a 5*5 full matrix A (in place)
 *
 * \param a  Pointer to the REAL array which stands a 5*5 matrix 
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   05/01/2010
 */
void fasp_blas_smat_inv_nc5 (REAL *a)
{ 
    const REAL a0=a[0],   a1=a[1],   a2=a[2],   a3=a[3],   a4=a[4];
    const REAL a5=a[5],   a6=a[6],   a7=a[7],   a8=a[8],   a9=a[9];
    const REAL a10=a[10], a11=a[11], a12=a[12], a13=a[13], a14=a[14];
    const REAL a15=a[15], a16=a[16], a17=a[17], a18=a[18], a19=a[19];
    const REAL a20=a[20], a21=a[21], a22=a[22], a23=a[23], a24=a[24];
    
    REAL det0, det1, det2, det3, det4, det, det_inv;
    
    det0  =  a6  * ( a12 * (a18*a24-a19*a23) + a17 * (a14*a23-a13*a24) + a22 * (a13*a19 - a14*a18) );
    det0 +=  a11 * ( a7  * (a19*a23-a18*a24) + a17 * (a8*a24 -a9*a23 ) + a22 * (a9*a18  - a8*a19)  );
    det0 +=  a16 * ( a7  * (a13*a24-a14*a23) + a12 * (a9*a23 -a8*a24 ) + a22 * (a8*a14  - a9*a13)  ); 
    det0 +=  a21 * ( a17 * (a9*a13 -a8*a14 ) + a7  * (a14*a18-a13*a19) + a12 * (a8*a19  - a9*a18)  );
    
    det1  =  a1  * ( a22 * (a14*a18-a13*a19) + a12 * (a19*a23-a18*a24) + a17 * (a13*a24 - a14*a23) );
    det1 +=  a11 * ( a17 * (a4*a23 - a3*a24) + a2  * (a18*a24-a19*a23) + a22 * (a3*a19  - a4*a18)  ); 
    det1 +=  a16 * ( a12 * (a3*a24 - a4*a23) + a2  * (a14*a23-a13*a24) + a22 * (a4*a13  - a3*a14)  );
    det1 +=  a21 * ( a2  * (a13*a19-a14*a18) + a12 * (a4*a18 -a3*a19 ) + a17 * (a3*a14  - a4*a13)  ); 
    
    det2  =  a1  * ( a7 * (a18*a24-a19*a23) + a17 * (a9*a23-a8*a24) + a22 * (a8*a19 - a9*a18) );
    det2 +=  a6  * ( a2 * (a19*a23-a18*a24) + a17 * (a3*a24-a4*a23) + a22 * (a4*a18 - a3*a19) );    
    det2 +=  a16 * ( a2 * (a8*a24 -a9*a23 ) + a7  * (a4*a23-a3*a24) + a22 * (a3*a9  - a4*a8)  );
    det2 +=  a21 * ( a7 * (a3*a19 -a4*a18 ) + a2  * (a9*a18-a8*a19) + a17 * (a4*a8  - a3*a9)  );
    
    det3  =  a1  * ( a12* (a8*a24 -a9*a23)  + a7  * (a14*a23-a13*a24) + a22 * (a9*a13 - a8*a14) );
    det3 +=  a6  * ( a2 * (a13*a24-a14*a23) + a12 * (a4*a23 -a3*a24 ) + a22 * (a3*a14 - a4*a13) );    
    det3 +=  a11 * ( a7 * (a3*a24 -a4*a23)  + a2  * (a9*a23 -a8*a24 ) + a22 * (a4*a8  - a3*a9)  );
    det3 +=  a21 * ( a2 * (a8*a14 -a9*a13)  + a7  * (a4*a13 -a3*a14 ) + a12 * (a3*a9  - a4*a8)  );
    
    det4  =  a1  * ( a7 * (a13*a19-a14*a18) + a12 * (a9*a18 -a8*a19 ) + a17 * (a8*a14 - a9*a13) );
    det4 +=  a6  * ( a12* (a3*a19 -a4*a18 ) + a17 * (a4*a13 -a3*a14 ) + a2  * (a14*a18- a13*a19));
    det4 +=  a11 * ( a2 * (a8*a19 -a9*a18 ) + a7  * (a4*a18 -a3*a19 ) + a17 * (a3*a9  - a4*a8)  );
    det4 +=  a16 * ( a7 * (a3*a14 -a4*a13 ) + a2  * (a9*a13 -a8*a14 ) + a12 * (a4*a8  - a3*a9)  );
    
    det = det0*a0 + det1*a5+ det2*a10 + det3*a15 + det4*a20;
    
    if (ABS(det)<1e-22) {
        printf("### WARNING: Matrix is nearly singular! det = %e\n", det);
        /*
         printf("##----------------------------------------------\n");
         printf("## %12.5e %12.5e %12.5e %12.5e %12.5e\n", a0,  a1,  a2,  a3,  a4);
         printf("## %12.5e %12.5e %12.5e %12.5e %12.5e\n", a5,  a6,  a7,  a8,  a9);
         printf("## %12.5e %12.5e %12.5e %12.5e %12.5e\n", a10, a11, a12, a13, a14);
         printf("## %12.5e %12.5e %12.5e %12.5e %12.5e\n", a15, a16, a17, a18, a19);
         printf("## %12.5e %12.5e %12.5e %12.5e %12.5e\n", a20, a21, a22, a23, a24);
         printf("##----------------------------------------------\n");
         getchar(); 
         */
    }
    
    det_inv = 1/det;
    
    a[0]  = a6  * (a12*a18*a24-a12*a19*a23-a17*a13*a24+a17*a14*a23+a22*a13*a19-a22*a14*a18);
    a[0] += a11 * (a7*a19*a23 -a7*a18*a24 +a17*a8*a24 -a17*a9*a23 -a22*a8*a19 +a22*a9*a18);
    a[0] += a16 * (a7*a13*a24 -a7*a14*a23 -a12*a8*a24 +a12*a9*a23 +a22*a8*a14 -a22*a9*a13);
    a[0] += a21 * (a7*a14*a18 -a7*a13*a19 +a12*a8*a19 -a12*a9*a18 -a17*a8*a14 +a17*a9*a13);
    a[0] *= det_inv;
    
    a[1]  = a1  * (a12*a19*a23-a12*a18*a24+a22*a14*a18-a17*a14*a23-a22*a13*a19+a17*a13*a24);
    a[1] += a11 * (a22*a3*a19 +a2*a18*a24 -a17*a3*a24 -a22*a4*a18 -a2*a19*a23 +a17*a4*a23 );
    a[1] += a16 * (a12*a3*a24 -a12*a4*a23 -a22*a3*a14 +a2*a14*a23 +a22*a4*a13 -a2*a13*a24 );
    a[1] += a21 * (a12*a4*a18 -a12*a3*a19 -a2*a14*a18 -a17*a4*a13 +a2*a13*a19 +a17*a3*a14 );
    a[1] *= det_inv;
    
    a[2]  = a1  * (a7*a18*a24-a7*a19*a23-a17*a8*a24+a17*a9*a23+a22*a8*a19-a22*a9*a18);
    a[2] += a6  * (a2*a19*a23-a2*a18*a24+a17*a3*a24-a17*a4*a23-a22*a3*a19+a22*a4*a18);
    a[2] += a16 * (a2*a8*a24 -a2*a9*a23 -a7*a3*a24 +a7*a4*a23 +a22*a3*a9 -a22*a4*a8);
    a[2] += a21 * (a2*a9*a18 -a2*a8*a19 +a7*a3*a19 -a7*a4*a18 -a17*a3*a9 +a17*a4*a8);
    a[2] *= det_inv;
    
    a[3]  = a1  * (a12*a8*a24-a12*a9*a23+a7*a14*a23-a7*a13*a24+a22*a9*a13-a22*a8*a14);
    a[3] += a6  * (a12*a4*a23-a12*a3*a24+a22*a3*a14-a22*a4*a13+a2*a13*a24-a2*a14*a23);
    a[3] += a11 * (a7*a3*a24 -a7*a4*a23 +a22*a4*a8 -a22*a3*a9 +a2*a9*a23 -a2*a8*a24);
    a[3] += a21 * (a12*a3*a9 -a12*a4*a8 +a2*a8*a14 -a2*a9*a13 +a7*a4*a13 -a7*a3*a14 );
    a[3] *= det_inv;
    
    a[4]  = a1  * (a7*a13*a19-a7*a14*a18-a12*a8*a19+a12*a9*a18+a17*a8*a14-a17*a9*a13);
    a[4] += a6  * (a2*a14*a18-a2*a13*a19+a12*a3*a19-a12*a4*a18-a17*a3*a14+a17*a4*a13);
    a[4] += a11 * (a2*a8*a19 -a2*a9*a18 -a7*a3*a19 +a7*a4*a18 +a17*a3*a9 -a17*a4*a8 );
    a[4] += a16 * (a2*a9*a13 -a2*a8*a14 +a7*a3*a14 -a7*a4*a13 -a12*a3*a9 +a12*a4*a8 );
    a[4] *= det_inv;
    
    a[5]  = a5  * (a12*a19*a23-a12*a18*a24+a22*a14*a18-a22*a13*a19+a17*a13*a24-a17*a14*a23);
    a[5] += a20 * (a12*a9*a18 -a12*a8*a19 +a7*a13*a19 -a18*a7*a14 +a17*a8*a14 -a9*a17*a13 );
    a[5] += a15 * (a22*a9*a13 -a12*a9*a23 +a12*a24*a8 +a7*a14*a23 -a24*a7*a13 -a22*a14*a8 );
    a[5] += a10 * (a18*a7*a24 -a18*a22*a9 -a17*a8*a24 +a17*a9*a23 +a22*a8*a19 -a19*a23*a7 );
    a[5] *= det_inv;
    
    a[6]  = a2  * (a19*a23*a10-a14*a23*a15    -a18*a24*a10+a18*a14*a20-a13*a19*a20+a24*a13*a15);
    a[6] += a12 * (a18*a0*a24 -a18*a20*a4 +a3*a19*a20 -a19*a23*a0 +a4*a23*a15 -a24*a15*a3 );
    a[6] += a17 * (a4*a13*a20 -a13*a24*a0 +a14*a23*a0 -a3*a14*a20 +a24*a3*a10 -a4*a23*a10 );
    a[6] += a22 * (a14*a15*a3 -a18*a14*a0 +a18*a4*a10 -a4*a13*a15 +a13*a19*a0 -a3*a19*a10 );
    a[6] *= det_inv;
    
    a[7]  = a0  * (a18*a9*a22 -a18*a24*a7 +a19*a23*a7 -a9*a23*a17 +a24*a8*a17 -a8*a19*a22 );
    a[7] += a5  * (a2*a18*a24 -a2*a19*a23 +a17*a4*a23 -a17*a3*a24 +a22*a3*a19 -a22*a4*a18 );
    a[7] += a15 * (a4*a8*a22  -a3*a9*a22  -a24*a8*a2  +a9*a23*a2  -a4*a23*a7  +a24*a3*a7  );
    a[7] += a20 * (a18*a4*a7  -a18*a9*a2  +a9*a3*a17  -a4*a8*a17  +a8*a19*a2  -a3*a19*a7  );
    a[7] *= det_inv;
    
    a[8]  = a0  * (a12*a9*a23 -a12*a24*a8 +a22*a14*a8 -a7*a14*a23 +a24*a7*a13 -a9*a22*a13 );
    a[8] += a5  * (a12*a3*a24 -a12*a4*a23 -a22*a3*a14 +a2*a14*a23 -a2*a13*a24 +a22*a4*a13 );
    a[8] += a10 * (a22*a9*a3  -a4*a22*a8  +a4*a7*a23  -a2*a9*a23  +a24*a2*a8  -a7*a24*a3  );
    a[8] += a20 * (a7*a14*a3  -a4*a7*a13  +a9*a2*a13  +a12*a4*a8  -a12*a9*a3  -a2*a14*a8  );
    a[8] *= det_inv;
    
    a[9]  = a0  * (a12*a8*a19-a12*a18*a9+a18*a7*a14-a8*a17*a14+a17*a13*a9-a7*a13*a19);
    a[9] += a5  * (a2*a13*a19-a2*a14*a18-a12*a3*a19+a12*a4*a18+a17*a3*a14-a17*a4*a13);
    a[9] += a10 * (a18*a2*a9-a18*a7*a4+a3*a7*a19-a2*a8*a19+a17*a8*a4-a3*a17*a9);
    a[9] += a15 * (a8*a2*a14-a12*a8*a4+a12*a3*a9-a3*a7*a14+a7*a13*a4-a2*a13*a9);
    a[9] *= det_inv;
    
    a[10] = a5  * (a18*a24*a11-a24*a13*a16+a14*a23*a16-a19*a23*a11+a13*a19*a21-a18*a14*a21);
    a[10]+= a10 * (a19*a23*a6-a9*a23*a16+a24*a8*a16-a8*a19*a21+a18*a9*a21-a18*a24*a6);
    a[10]+= a15 * (a24*a13*a6-a14*a23*a6-a24*a8*a11+a9*a23*a11+a14*a8*a21-a13*a9*a21);
    a[10]+= a20 * (a18*a14*a6-a18*a9*a11+a8*a19*a11-a13*a19*a6+a9*a13*a16-a14*a8*a16);
    a[10]*= det_inv;
    
    a[11] = a4  * (a21*a13*a15-a11*a23*a15+a16*a23*a10-a13*a16*a20+a18*a11*a20-a18*a21*a10);
    a[11]+= a14 * (a18*a0*a21-a1*a18*a20+a16*a3*a20-a23*a0*a16+a1*a23*a15-a21*a3*a15);
    a[11]+= a19 * (a1*a13*a20-a1*a23*a10+a23*a0*a11+a21*a3*a10-a11*a3*a20-a13*a0*a21);
    a[11]+= a24 * (a13*a0*a16-a18*a0*a11+a11*a3*a15+a1*a18*a10-a1*a13*a15-a16*a3*a10);
    a[11]*= det_inv;
    
    a[12] = a4  * (a5*a21*a18-a18*a20*a6+a20*a16*a8-a5*a16*a23+a15*a6*a23-a21*a15*a8);
    a[12]+= a9  * (a1*a20*a18-a1*a15*a23+a0*a16*a23-a18*a0*a21-a20*a16*a3+a15*a21*a3);
    a[12]+= a19 * (a20*a6*a3-a5*a21*a3+a0*a21*a8-a23*a0*a6+a1*a5*a23-a1*a20*a8);
    a[12]+= a24 * (a1*a15*a8-a0*a16*a8+a18*a0*a6-a1*a5*a18+a5*a16*a3-a6*a15*a3);
    a[12]*= det_inv;
    
    a[13] = a0  * (a24*a11*a8-a6*a24*a13+a21*a9*a13-a11*a9*a23+a14*a6*a23-a14*a21*a8);
    a[13]+= a1  * (a5*a13*a24-a5*a14*a23+a14*a20*a8+a10*a9*a23-a24*a10*a8-a20*a9*a13);
    a[13]+= a3  * (a6*a10*a24-a10*a9*a21+a5*a14*a21-a5*a24*a11+a20*a9*a11-a14*a6*a20);
    a[13]+= a4  * (a5*a11*a23-a5*a21*a13+a21*a10*a8-a6*a10*a23+a20*a6*a13-a11*a20*a8);
    a[13]*= det_inv;
    
    a[14] = a0  * (a13*a19*a6-a14*a18*a6-a11*a19*a8+a14*a16*a8+a11*a18*a9-a13*a16*a9);
    a[14]+= a1  * (a14*a18*a5-a13*a19*a5+a10*a19*a8-a14*a15*a8-a10*a18*a9+a13*a15*a9);
    a[14]+= a3  * (a11*a19*a5-a11*a15*a9+a10*a16*a9-a10*a19*a6+a14*a15*a6-a14*a16*a5);
    a[14]+= a4  * (a11*a15*a8-a11*a18*a5+a13*a16*a5-a13*a15*a6+a10*a18*a6-a10*a16*a8); 
    a[14]*= det_inv;
    
    a[15] = a5  * (a19*a22*a11-a24*a17*a11+a12*a24*a16-a22*a14*a16-a12*a19*a21+a17*a14*a21);
    a[15]+= a10 * (a24*a17*a6-a19*a22*a6-a24*a7*a16+a22*a9*a16+a19*a7*a21-a17*a9*a21);
    a[15]+= a15 * (a22*a14*a6-a9*a22*a11+a24*a7*a11-a12*a24*a6-a7*a14*a21+a12*a9*a21);
    a[15]+= a20 * (a12*a19*a6-a17*a14*a6-a19*a7*a11+a9*a17*a11+a7*a14*a16-a12*a9*a16);
    a[15]*= det_inv;
    
    a[16] = a0  * (a11*a17*a24-a11*a19*a22-a12*a16*a24+a12*a19*a21+a14*a16*a22-a14*a17*a21);
    a[16]+= a1  * (a10*a19*a22-a10*a17*a24+a12*a15*a24-a12*a19*a20-a14*a15*a22+a14*a17*a20);
    a[16]+= a2  * (a10*a16*a24-a10*a19*a21-a11*a15*a24+a11*a19*a20+a14*a15*a21-a14*a16*a20);
    a[16]+= a4  * (a10*a17*a21*+a11*a15*a22-a11*a17*a20-a12*a15*a21+a12*a16*a20-a10*a16*a22);
    a[16]*= det_inv;    
    
    a[17] = a0  * (a21*a9*a17-a6*a24*a17+a19*a6*a22-a0*a16*a9*a22+a24*a16*a7-a19*a21*a7);
    a[17]+= a1  * (a5*a24*a17-a5*a19*a22+a19*a20*a7-a20*a9*a17+a15*a9*a22-a24*a15*a7);
    a[17]+= a2  * (a5*a19*a21-a19*a6*a20-a5*a24*a16+a24*a6*a15-a15*a9*a21+a20*a9*a16);
    a[17]+= a4  * (a16*a5*a22-a6*a15*a22+a20*a6*a17-a5*a21*a17-a6*a15*a22+a21*a15*a7-a16*a20*a7);
    a[17]*= det_inv;    
    
    a[18] = a0  * (a12*a24*a6 - a14*a22*a6 - a11*a24*a7 + a14*a21*a7 + a11*a22*a9 - a12*a21*a9);
    a[18]+= a1  * (a14*a22*a5 - a12*a24*a5 + a10*a24*a7 - a14*a20*a7 - a10*a22*a9 + a12*a20*a9);
    a[18]+= a2  * (a11*a24*a5 - a11*a20*a9 + a14*a20*a6 - a14*a21*a5 + a10*a21*a9  - a10*a24*a6);  
    a[18]+= a4  * (a11*a20*a7 - a11*a22*a5 + a12*a21*a5 + a10*a22*a6 - a12*a20*a6 - a10*a21*a7);
    a[18]*= det_inv;    
    
    a[19] = a0  * (a12*a16*a9-a6*a12*a19+a6*a17*a14-a17*a11*a9+a11*a7*a19-a16*a7*a14);
    a[19]+= a1  * (a5*a12*a19-a5*a17*a14-a12*a15*a9+a17*a10*a9+a15*a7*a14-a10*a7*a19);
    a[19]+= a2  * (a11*a15*a9-a5*a11*a19+a5*a16*a14-a6*a15*a14+a6*a10*a19-a16*a10*a9);
    a[19]+= a4  * (a5*a17*a11-a5*a12*a16+a12*a6*a15+a10*a7*a16-a17*a6*a10-a15*a7*a11);
    a[19]*= det_inv;    
    
    a[20] = a5  * (a12*a18*a21-a12*a23*a16+a22*a13*a16-a18*a22*a11+a23*a17*a11-a17*a13*a21);
    a[20]+= a15 * (a12*a23*a6-a12*a8*a21+a8*a22*a11-a23*a7*a11+a7*a13*a21-a22*a13*a6);
    a[20]+= a20 * (a12*a8*a16-a12*a18*a6+a18*a7*a11-a8*a17*a11+a17*a13*a6-a7*a13*a16);
    a[20]+= a10 * (a17*a8*a21-a22*a8*a16-a18*a7*a21+a18*a22*a6+a23*a7*a16-a23*a17*a6);
    a[20]*= det_inv;    
    
    a[21] = a0  * (a12*a23*a16-a12*a18*a21+a17*a13*a21+a18*a22*a11-a23*a17*a11-a22*a13*a16);
    a[21]+= a1  * (a12*a18*a20-a12*a23*a15+a22*a13*a15+a23*a17*a10-a17*a13*a20-a18*a22*a10);
    a[21]+= a2  * (a18*a21*a10-a18*a11*a20-a21*a13*a15+a16*a13*a20-a23*a16*a10+a23*a11*a15);
    a[21]+= a3  * (a17*a11*a20-a12*a16*a20+a12*a21*a15-a21*a17*a10-a22*a11*a15+a16*a22*a10);
    a[21]*= det_inv;    
    
    a[22] = a0  * (a18*a21*a7-a18*a6*a22+a23*a6*a17+a16*a8*a22-a21*a8*a17-a23*a16*a7);
    a[22]+= a1  * (a5*a18*a22-a5*a23*a17-a15*a8*a22+a20*a8*a17-a18*a20*a7+a23*a15*a7);
    a[22]+= a3  * (a16*a20*a7+a6*a15*a22-a6*a20*a17-a5*a16*a22+a5*a21*a17-a21*a15*a7);
    a[22]+= a2  * (a5*a23*a16-a5*a18*a21+a18*a6*a20+a15*a8*a21-a20*a8*a16-a23*a6*a15);
    a[22]*= det_inv;    
    
    a[23] = a0  * (a12*a21*a8-a22*a11*a8+a11*a7*a23-a6*a12*a23-a21*a7*a13+a6*a22*a13);
    a[23]+= a1  * (a5*a12*a23-a5*a22*a13-a10*a7*a23+a20*a7*a13+a22*a10*a8-a12*a20*a8);
    a[23]+= a2  * (a5*a21*a13+a11*a20*a8+a6*a10*a23-a5*a11*a23-a21*a10*a8-a6*a20*a13);
    a[23]+= a3  * (a5*a22*a11-a5*a12*a21+a10*a7*a21-a22*a6*a10-a20*a7*a11+a12*a6*a20);
    a[23]*= det_inv;    
    
    a[24] = a0  * (a17*a11*a8-a11*a7*a18+a6*a12*a18-a12*a16*a8+a16*a7*a13-a6*a17*a13);
    a[24]+= a1  * (a5*a17*a13-a5*a12*a18+a10*a7*a18+a12*a15*a8-a17*a10*a8-a15*a7*a13);
    a[24]+= a2  * (a5*a11*a18-a5*a16*a13+a16*a10*a8+a6*a15*a13-a11*a15*a8-a6*a10*a18);
    a[24]+= a3  * (a5*a12*a16+a17*a6*a10-a5*a17*a11-a12*a6*a15-a10*a7*a16+a15*a7*a11);
    a[24]*= det_inv;    
}

/**
 * \fn void fasp_blas_smat_inv_nc7 (REAL *a)
 *
 * \brief Compute the inverse matrix of a 7*7 matrix a
 *
 * \param a   Pointer to the REAL array which stands a 7*7 matrix 
 * 
 * \note This is NOT implemented yet!
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   05/01/2010
 */
void fasp_blas_smat_inv_nc7 (REAL *a)
{    
    fasp_blas_smat_inv(a,7);
}

/**
 * \fn INT fasp_blas_smat_inv (REAL *a, const INT n)
 *
 * \brief Compute the inverse matrix of a small full matrix a
 *
 * \param a   Pointer to the REAL array which stands a n*n matrix
 * \param n   Dimension of the matrix
 * 
 * \author Xiaozhe Hu, Shiquan Zhang
 * \date   04/21/2010  
 */
INT fasp_blas_smat_inv (REAL *a, 
                        const INT n)
{ 
    
    switch (n) {
            
        case 2:
            fasp_blas_smat_inv_nc2(a);
            break;
            
        case 3:
            fasp_blas_smat_inv_nc3(a);
            break;
            
        case 4:
            fasp_blas_smat_inv_nc4(a);
            break;
            
        case 5: 
            fasp_blas_smat_inv_nc5(a);
            break;    
            
        default: { 
            INT i,j,k,l,u,kn,in;
            REAL alinv;
            
            for (k=0; k<n; ++k) {
                
                kn = k*n;
                l  = kn+k;
                
                if (ABS(a[l]) < SMALLREAL) {
                    printf("### WARNING: Diagonal entry is close to zero! diag_%d=%e\n",k,a[l]);
                    a[l] = SMALLREAL;
                }
                alinv = 1.0/a[l];
                a[l] = alinv;
                
                for (j=0; j<k; ++j) {
                    u = kn+j; a[u] *= alinv;
                }
                
                for (j=k+1; j<n; ++j) {
                    u = kn+j; a[u] *= alinv;
                }
                
                for (i=0; i<k; ++i) {
                    in = i*n;
                    for (j=0; j<n; ++j)
                        if (j!=k) { 
                            u = in+j; a[u] -= a[in+k]*a[kn+j];
                        } // end if (j!=k)
                }
                
                for (i=k+1; i<n; ++i) {
                    in = i*n;
                    for (j=0; j<n; ++j)
                        if (j!=k) { 
                            u = in+j; a[u] -= a[in+k]*a[kn+j];
                        } // end if (j!=k)
                }
                
                for (i=0; i<k; ++i) {
                    u=i*n+k; a[u] *= -alinv;
                }
                
                for (i=k+1; i<n; ++i) {
                    u=i*n+k; a[u] *= -alinv;
                }    
                
            } // end for (k=0; k<n; ++k)
            break;
        }
            
    }
    
    return SUCCESS;
}

/**
 * \fn void fasp_iden_free (idenmat *A)
 *
 * \brief Free idenmat sparse matrix data memeory space
 *
 * \param A      Pointer to the idenmat matrix
 *
 * \author Chensong Zhang
 * \date   2010/04/03 
 */
void fasp_iden_free (idenmat *A)
{    
    INT i;
    
    if (A==NULL) return;
    
    for (i=0;i<A->row;++i) fasp_mem_free(A->val[i]);
    fasp_mem_free(A->val); A->val = NULL; A->row = 0;
}

/**
 * \fn void fasp_smat_identity_nc2 (REAL *a)
 *
 * \brief Set a 2*2 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 2*2 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2011/11/18
 */
void fasp_smat_identity_nc2 (REAL *a)
{
    memset(a, 0X0, 4*sizeof(REAL));
    
    a[0] = 1.0; a[3] = 1.0; 
}    

/**
 * \fn void fasp_smat_identity_nc3 (REAL *a)
 *
 * \brief Set a 3*3 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 3*3 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity_nc3 (REAL *a)
{
    memset(a, 0X0, 9*sizeof(REAL));
    
    a[0] = 1.0; a[4] = 1.0; a[8] = 1.0;
}    

/**
 * \fn void fasp_smat_identity_nc5 (REAL *a)
 *
 * \brief Set a 5*5 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 5*5 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity_nc5 (REAL *a)
{
    memset(a, 0X0, 25*sizeof(REAL));
    
    a[0]  = 1.0; 
    a[6]  = 1.0; 
    a[12] = 1.0; 
    a[18] = 1.0; 
    a[24] = 1.0;
}    

/**
 * \fn void fasp_smat_identity_nc7 (REAL *a)
 *
 * \brief Set a 7*7 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 7*7 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity_nc7 (REAL *a)
{
    memset(a, 0X0, 49*sizeof(REAL));
    
    a[0]  = 1.0; 
    a[8]  = 1.0; 
    a[16] = 1.0; 
    a[24] = 1.0;
    a[32] = 1.0; 
    a[40] = 1.0; 
    a[48] = 1.0;
}    

/**
 * \fn void fasp_smat_identity (REAL *a, INT n, INT n2)
 *
 * \brief Set a n*n full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a n*n full matrix 
 * \param n      Size of full matrix
 * \param n2     Length of the REAL vector which stores the n*n full matrix
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity (REAL *a, 
                         INT n, 
                         INT n2)
{
    memset(a, 0X0, n2*sizeof(REAL));
    
    switch (n) {

    case 2: {
        a[0] = 1.0;
        a[3] = 1.0;
    }
        break;
    
    case 3: {
        a[0] = 1.0;
        a[4] = 1.0;
        a[8] = 1.0;
    }
        break;
    
    case 4: {
        a[0] = 1.0;
        a[5] = 1.0;
        a[10] = 1.0;
        a[15] = 1.0;
    }
        break;
    
    case 5: {
        a[0] = 1.0;
        a[6] = 1.0;
        a[12] = 1.0;
        a[18] = 1.0;
        a[24] = 1.0;
    }
        break;
    
    case 6: {
        a[0] = 1.0;
        a[7] = 1.0;
        a[14] = 1.0;
        a[21] = 1.0;
        a[28] = 1.0;
        a[35] = 1.0;
    }
        break;
    
    case 7: {
        a[0] = 1.0;
        a[8] = 1.0;
        a[16] = 1.0;
        a[24] = 1.0;
        a[32] = 1.0;
        a[40] = 1.0;
        a[48] = 1.0;
    }
        break;
    
    default: {
        INT l;
        for (l = 0; l < n; l ++) a[l*n+l] = 1.0;
    }
        break;
    }
    
}    

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
