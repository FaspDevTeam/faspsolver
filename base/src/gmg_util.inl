/*! \file  gmg_util.inl
 *
 *  \brief Routines for GMG solvers
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn void compute_r_1d(REAL *u, REAL *b, REAL *r, INT k, INT *level)
 * \brief compute residue vector r of 1D problem
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param r            Pointer to the residue vector
 * \param k            Level k
 * \param level        Pointer to the start position of each level
 *
 * \author Ziteng Wang
 * \date   2013-06-07
 */
static void compute_r_1d(REAL *u,
                         REAL *b,
                         REAL *r,
                         INT   k,
                         INT *level)
{
    INT i,n;
    const INT levelk = level[k];
    
    n = level[k+1]-level[k];
    
    for (i = 1; i < n-1; i++) {
        r[levelk+i] = b[levelk+i]-2*u[levelk+i]+u[levelk+i+1]+u[levelk+i-1];
    }
}

/**
 * \fn void compute_r_2d(REAL *u, REAL *b, REAL *r, INT k, INT *level
 *                       INT *nxk, INT *nyk)
 * \brief compute residue vector r of 2D problem
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param r            Pointer to the residue vector
 * \param k            Level k
 * \param level        Pointer to the start position of each level
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 *
 * \author Ziteng Wang
 * \date   2013-06-07
 */
static void compute_r_2d(REAL *u,
                         REAL *b,
                         REAL *r,
                         INT k,
                         INT *level,
                         INT *nxk,
                         INT *nyk)
{
    INT i,j,k1;
    const INT nykk = nyk[k];
    const INT nxkk = nxk[k];
    const INT levelk = level[k];
    
    for (i = 1; i < nykk-1; i++) {
        k1 = levelk+i*nxkk;
        for (j = 1; j < nxkk-1; j++) {
            r[k1+j] = b[k1+j]-4*u[k1+j]+u[k1+j+1]+u[k1+j-1]+u[k1+nxkk+j]+u[k1-nxkk+j];
        }
    }
}

/**
 * \fn void compute_r_3d(REAL *u, REAL *b, REAL *r, INT k, INT *level,
 *                       INT *nxk, INT *nyk, INT *nzk)
 * \brief compute residue vector r of 3D problem
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param r            Pointer to the residue vector
 * \param k            Level k
 * \param level        Pointer to the start position of each level
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   2013-06-07
 */
static void compute_r_3d(REAL *u,
                         REAL *b,
                         REAL *r,
                         INT k,
                         INT *level,
                         INT *nxk,
                         INT *nyk,
                         INT *nzk)
{
    int i,j,h;
    int i0,j0,j1,j2,k0,k1,k2,k3,k4,k5,k6;
    const INT levelk = level[k];
    const INT nxkk = nxk[k];
    const INT nykk = nyk[k];
    const INT nzkk = nzk[k];
    
    // middle part of the cubic
    for (i = 1; i < nzkk-1; i++) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j++) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h++) {
                k0 = j0+h;
                k1 = j0+h-1;
                k2 = j0+h+1;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxkk*nykk;
                k6 = k0-nxkk*nykk;
                r[k0]=b[k0]-6*u[k0]+u[k1]+u[k2]+u[k3]+u[k4]+u[k5]+u[k6];
            }
        }
    }
}

/**
 * \fn REAL computenorm(REAL *r, INT *level, INT k)
 * \brief compute L2 norm of vector r
 *
 * \param r     Pointer to the residue vector
 * \param level Pointer to the start position of each level
 * \param k     level k
 *
 * \return L2 Norm
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static REAL computenorm(REAL *r,
                        INT *level,
                        INT k)
{
    INT i,n;
    REAL squarnorm;
    
    squarnorm = 0.0;
    n = level[k+1]-level[k];
    for (i = 1; i < n; i++) {
        squarnorm = squarnorm + r[level[k]+i]*r[level[k]+i];
    }
    squarnorm = sqrt(squarnorm);
    return squarnorm;
}

/**
 * \fn void xequaly(REAL *x, REAL *y, INT *level, INT k)
 * \brief x = y
 *
 * \param x     vector x
 * \param y     vector x
 * \param level Pointer to the start position of each level
 * \param k     Level k
 *
 * \author Ziteng Wang
 * \date 06/07/2013
 */
static void xequaly(REAL *x,
                    REAL *y,
                    INT *level,
                    INT k)
{
    INT i,n;
    const INT levelk = level[k];
    
    n = level[k+1] - level[k];
    
    for (i = 0; i < n; i++) {
        x[levelk+i] = y[levelk+i];
    }
}

/**
 * \fn void xequalypcz(REAL *x, REAL *y, REAL c, REAL *z, INT *level, INT k)
 * \brief x = y+c*z
 *
 * \param x     vector x
 * \param y     vector x
 * \param c     coefficient c
 * \param z     vector z
 * \param level Pointer to the start position of each level
 * \param k     Level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void xequalypcz(REAL *x,
                       REAL *y,
                       REAL c,
                       REAL *z,
                       INT *level,
                       INT k)
{
    INT i, n;
    const INT levelk = level[k];
    
    n = level[k+1]-level[k];
    
    for (i = 0; i < n; i++) {
        x[levelk+i] = y[levelk+i] + c*z[levelk+i];
    }
}

/**
 * \fn void xequalay_1d(REAL *x, REAL *y, INT *level, INT k);
 * \brief x = Ay
 *
 * \param x     vector x
 * \param y     vector x
 * \param level Pointer to the start position of each level
 * \param k     Level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void xequalay_1d(REAL *x,
                        REAL *y,
                        INT *level,
                        INT k)
{
    INT i,n;
    REAL *btemp = (REAL *)malloc(level[k+1]*sizeof(REAL));
    
    for (i = 0; i < level[k+1]; i++) {
        btemp[i] = 0.0;
    }
    
    // compute (-x)
    compute_r_1d(y, btemp, x, k, level);
    n = level[k+1]-level[k];
    for (i = 0; i < n; i++) {
        x[level[k]+i] = (-1)*x[level[k]+i];
    }
    
    free(btemp);
}

/**
 * \fn void xequalay_2d(REAL *x, REAL *y, INT *level, INT k,
 *                      INt *nxk, INT *nyk)
 * \brief x = Ay
 *
 * \param x     vector x
 * \param y     vector x
 * \param level Pointer to the start position of each level
 * \param k     Level k
 * \param nxk   Number of grids in x direction in level k
 * \param nyk   Number of grids in y direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void xequalay_2d(REAL *x,
                        REAL *y,
                        INT *level,
                        INT k,
                        INT *nxk,
                        INT *nyk)
{
    INT i,n;
    REAL *btemp = (REAL *)malloc(level[k+1]*sizeof(REAL));
    
    for (i = 0; i < level[k+1]; i++) {
        btemp[i] = 0.0;
    }
    
    // compute (-x)
    compute_r_2d(y, btemp, x, k, level, nxk, nyk);
    n = level[k+1]-level[k];
    for (i = 0; i < n; i++) {
        x[level[k]+i] = (-1)*x[level[k]+i];
    }
    
    free(btemp);
}

/**
 * \fn void xequalay_3d(REAL *x, REAL *y, INT *level, INT k,
 *                      INt *nxk, INT *nyk)
 * \brief x = Ay 3D
 *
 * \param x     vector x
 * \param y     vector x
 * \param level Pointer to the start position of each level
 * \param k     Level k
 * \param nxk   Number of grids in x direction in level k
 * \param nyk   Number of grids in y direction in level k
 * \param nzk   Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void xequalay_3d(REAL *x,
                        REAL *y,
                        INT *level,
                        INT k,
                        INT *nxk,
                        INT *nyk,
                        INT *nzk)
{
    const INT levelk = level[k];
    const INT levelk1 = level[k+1];
    
    INT i,n;
    REAL *btemp = (REAL *)malloc(levelk1*sizeof(REAL));
    
    for (i = levelk; i < levelk1; i++) {
        btemp[i] = 0.0;
    }
    // compute (-x)
    compute_r_3d(y, btemp, x, k, level, nxk, nyk, nzk);
    n = level[k+1]-level[k];
    for (i = 0; i < n; i++) {
        x[levelk+i] = (-1)*x[levelk+i];
    }
    
    free(btemp);
}

/**
 * \fn REAl innerproductxy(REAL *x, REAL *y, INT *level, INT k);
 * \brief <x,y>
 *
 * \param x     vector x
 * \param y     vector x
 * \param level Pointer to the start position of each level
 * \param k     Level k
 *
 * \return Innerproduct of x,y
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static REAL innerproductxy(REAL *x,
                           REAL *y,
                           INT *level,
                           INT k)
{
    INT i,n;
    REAL innerproduct;
    const INT levelk = level[k];
    const INT levelk1 = level[k+1];
    
    innerproduct = 0.0;
    n = levelk1 - levelk;
    
    for	(i = 0; i < n; i++) {
        innerproduct = innerproduct + x[levelk+i]*y[levelk+i];
    }
    
    return innerproduct;
}

/**
 * \fn REAL energynormu1d(REAL *u, INT nx);
 * \brief Energy norm u^t*A*u 1D
 *
 * \param u     vector u
 * \param nx    number of grids in x direction
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static REAL energynormu1d(REAL *u,
                          INT nx)
{
    REAL *utemp;
    REAL  error;
    INT   level[2];
    
    level[0] = 0; level[1] = nx+1;
    utemp = (REAL *)malloc(level[1]*sizeof(REAL));
    fasp_array_set(level[1], utemp, 0.0);
    xequalay_1d(utemp, u, level, 0);
    error = innerproductxy(utemp, u, level, 0);
    
    free(utemp);
    
    return error;
}

/**
 * \fn REAL energynormu2d(REAL *u, INT nx, INT ny);
 * \brief Energy norm u^t*A*u 2D
 *
 * \param u     vector u
 * \param nx    number of grids in x direction
 * \param ny    number of grids in y direction
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static REAL energynormu2d(REAL *u,
                          INT nx,
                          INT ny)
{
    REAL *utemp;
    REAL error;
    INT nxk[1],nyk[1],level[2];
    
    level[0] = 0; level[1] = (nx+1)*(ny+1);
    nxk[0] = nx+1; nyk[0] = ny+1;
    utemp = (REAL *)malloc(level[1]*sizeof(REAL));
    fasp_array_set(level[1], utemp, 0.0);
    xequalay_2d(utemp, u, level, 0, nxk, nyk);
    error = innerproductxy(utemp, u, level, 0);
    
    free(utemp);
    
    return error;
}

/**
 * \fn REAL energynormu3d(REAL *u, INT nx, INT ny, INT nz);
 * \brief Energy norm u^t*A*u 3D
 *
 * \param u     vector u
 * \param nx    number of grids in x direction
 * \param ny    number of grids in y direction
 * \param nz    number of grids in z direction
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static REAL energynormu3d(REAL *u,
                          INT nx,
                          INT ny,
                          INT nz)
{
    REAL *utemp;
    REAL error;
    INT nxk[1],nyk[1],nzk[1],level[2];
    
    level[0] = 0; level[1] = (nx+1)*(ny+1)*(nz+1);
    nxk[0] = nx+1; nyk[0] = ny+1; nzk[0] = nz+1;
    utemp = (REAL *)malloc(level[1]*sizeof(REAL));
    fasp_array_set(level[1], utemp, 0.0);
    xequalay_3d(utemp, u, level, 0, nxk, nyk, nzk);
    error = innerproductxy(utemp, u, level, 0);
    
    free(utemp);
    
    return error;
}

/**
 * \fn void coarsergrid7pointrestriction2d(REAL *b, REAL *r, INT *level, const INT k,
 *                                         INT *nxk, INT *nyk)
 * \brief Restriction function in multigrid of 2D
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 *
 * Modified by Chensong on 05/05/2015: Fix array out-of-bound bug.
 */
static void coarsergrid7pointrestriction2d(REAL *b,
                                           REAL *r,
                                           INT *level,
                                           const INT k,
                                           INT *nxk,
                                           INT *nyk)
{
    INT i,j;
    int const nxkk1 = nxk[k+1];
    int const nykk1 = nyk[k+1];
    int const nxkk = nxk[k];
    int const nykk = nyk[k];
    const INT levelk1 = level[k+1];
    const INT levelk = level[k];
    int k1,k11,k2;
    
    for (i = 1; i < nykk1-1; i++) {
        k11 = levelk1+i*nxkk1;
        k1 = levelk+2*i*nxkk;
        for (j = 1; j < nxk[k+1]-1; j++) {
            k2 = k1+2*j;
            b[k11+j] = (r[k2]*2+r[k2+1]+r[k2-1]+r[k2+nxkk]+r[k2-nxkk]+r[k2+nxkk+1]+r[k2-nxkk-1])/2;
        }
        b[levelk1+(i+1)*nxkk1-1] = (r[levelk+(2*i+1)*nxkk-2]+r[levelk+2*i*nxkk-2])/2;
    }
    k11 = levelk1+(nykk1-1)*nxkk1;
    k1 = levelk+(nykk-1)*nxkk-1;
    for (j = 1; j < nxk[k+1]-1; j++) {
        b[k11+j] = (r[k1+2*j]+r[k1+2*j-nxkk])/2;
    }
    b[levelk1+nxkk1*nykk1-1] = r[levelk+(nykk-1)*nxkk-2]/2;
}

/**
 * \fn void coarsergrid7pointrestriction3d(REAL *b, REAL *r, INT *level, INT k,
 *                                         INT *nxk, INT *nyk, INT *nzk)
 * \brief Restriction function in multigrid of 3D
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void coarsergrid7pointrestriction3d(REAL *b,
                                           REAL *r,
                                           INT *level,
                                           INT k,
                                           INT *nxk,
                                           INT *nyk,
                                           INT *nzk)
{
    INT i,j,h;
    int i0,j0,k0,k1,k2,k3,k4,k5,k6,i01,j01;
    const INT levelk = level[k];
    const INT levelk1 = level[k+1];
    const INT nxkk = nxk[k];
    const INT nxkk1 = nxk[k+1];
    const INT nykk = nyk[k];
    const INT nykk1 = nyk[k+1];
    const INT nzkk1 = nzk[k+1];
    
    int nxyk = nxkk*nykk;
    int nxyk1 = nxkk1*nykk1;
    
    for (i = 1; i < nzkk1-1; i++) {
        i0 = levelk+2*i*nxyk;
        i01 = levelk1+i*nxyk1;
        for (j = 1; j < nykk1-1; j++) {
            j0 = i0+2*j*nxkk;
            j01 = i01+j*nxkk1;
            for (h = 1; h < nxkk1-1; h++) {
                k0 = j0+2*h;
                k1 = k0-nxkk;
                k2 = k0+nxkk;
                k3 = k0-nxyk;
                k4 = k0-nxyk-nxkk;
                k5 = k0+nxyk;
                k6 = k0+nxyk+nxkk;
                b[j01+h] = (r[k0]*2+r[k0-1]+r[k0+1]+r[k1]+r[k1-1]+r[k2]+r[k2+1]+r[k3]
                            +r[k3-1]+r[k4]+r[k4-1]+r[k5]+r[k5+1]+r[k6]+r[k6+1])/4;
            }
        }
    }
}

/**
 * \fn void finergridinterpolation2d(REAL *u, INT *level, INT k,
 *                                   INT *nxk, INT *nyk)
 * \brief Interpolation function in multigrid of 2D
 *
 * \param u            Pointer to the vector of DOFs
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param nxk          Number of grids in x direction
 * \param nyk          Number of grids in y direction
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void finergridinterpolation2d(REAL *u,
                                     INT *level,
                                     INT k,
                                     INT *nxk,
                                     INT *nyk)
{
    INT i,j;
    
    for (i = 0; i < nyk[k]-1; i++) {
        for (j = 0; j < nxk[k]-1; j++) {
            u[level[k-1]+2*(j+i*nxk[k-1])] += u[level[k]+j+i*nxk[k]];
            u[level[k-1]+2*(j+i*nxk[k-1])+1] += (u[level[k]+j+i*nxk[k]]+
                                                 u[level[k]+j+1+i*nxk[k]])/2;
            u[level[k-1]+(2*i+1)*nxk[k-1]+2*j] +=
            (u[level[k]+i*nxk[k]+j]+u[level[k]+(i+1)*nxk[k]+j])/2;
            u[level[k-1]+(2*i+1)*nxk[k-1]+2*j+1] +=
            (u[level[k]+i*nxk[k]+j]+u[level[k]+(i+1)*nxk[k]+1+j])/2;
        }
    }
}

/**
 * \fn void finergridinterpolation2d(REAL *u, INT *level, INT k,
 *                                   INT *nxk, INT *nyk)
 * \brief Interpolation function in multigrid of 2D
 *
 * \param u            Pointer to the vector of DOFs
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          NUmber of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void finergridinterpolation3d(REAL *u,
                                     INT *level,
                                     INT k,
                                     INT *nxk,
                                     INT *nyk,
                                     INT *nzk)
{
    INT i, j, h;
    int i0,j0,j1,j2,j3,k0,k1,k2,k3;
    int i01,j01,j11,j21,j31,k01,k11,k21,k31;
    const INT levelk = level[k];
    const INT levelk1 = level[k-1];
    const INT nxkk = nxk[k];
    const INT nxkk1 = nxk[k-1];
    const INT nykk = nyk[k];
    const INT nykk1 = nyk[k-1];
    const INT nzkk = nzk[k];
    
    int nxyk = nxkk*nykk;
    int nxyk1 = nxkk1*nykk1;
    for (i = 0; i < nzkk-1; i++) {
        i0 = levelk+i*nxkk*nykk;
        i01 = levelk1+2*i*nxkk1*nykk1; 
        for (j = 0; j < nykk-1; j++) {
            j0 = i0+j*nxkk;
            j01 = i01 + 2*j*nxkk1;
            j1 = j0+nxkk;
            j11 = j01 + nxkk1;
            j2 = j0 + nxyk;
            j21 = j01+nxyk1;
            j3 = j2+nxkk;
            j31 = j21+nxkk1;
            
            for (h = 0; h < nxkk-1; h++) {
                k01 = j01+2*h;
                k0 = j0+h;
                u[k01] += u[k0];
                u[k01+1] += (u[k0]+u[k0+1])/2;         
                k11 = j11+2*h;
                k1 = j1+h;
                u[k11] +=(u[k0]+u[k1])/2; 
                u[k11+1] +=(u[k0]+u[k1+1])/2;
                k21 = j21+2*h;
                k2 = j2+h;
                u[k21] += (u[k0]+u[k2])/2;
                u[k21+1] += (u[k0]+u[k2+1])/2;            
                k31 = j31+2*h;
                k3 = j3+h;
                u[k31] += (u[k0]+u[k3])/2;
                u[k31+1] += (u[k0]+u[k3+1])/2;
            }
        }
    }
}

/**
 * \fn void gsiteration_2color_2d(REAL *u, REAL *b, INT *level, INT k,
 *                                INT maxlevel, INT *nxk, INT nyk)
 * \brief 2 color G-S iteration of 2D problem
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void gsiteration_2color_2d(REAL *u,
                                  REAL *b,
                                  INT *level,
                                  INT k,
                                  INT maxlevel,
                                  INT *nxk,
                                  INT *nyk)
{
    INT h,i;
    const INT nxkk = nxk[k];
    const INT nykk = nyk[k];
    const INT levelk = level[k];
    int k1;
    
    // red
    for (h = 1; h < nykk-1; h = h+2) {
        k1 = levelk+nxkk*h;
        
        for (i = 1; i < nxkk-1; i=i+2) {
            u[k1+i] = (b[k1+i]+u[k1+i+1]+u[k1+i-1]+u[k1+nxkk+i]+u[k1-nxkk+i])/4;
        }
    }
    for (h = 2; h < nykk-1; h = h+2) {
        k1 = levelk+nxkk*h;
        
        for (i = 2; i < nxkk-1; i=i+2) {
            u[k1+i] = (b[k1+i]+u[k1+i+1]+u[k1+i-1]+u[k1+nxkk+i]+u[k1-nxkk+i])/4;
        }
    }
    // black
    for (h = 1; h < nykk-1; h = h+2) {
        k1 = levelk+nxkk*h;
        
        for (i = 2; i < nxkk-1; i=i+2) {
            u[k1+i] = (b[k1+i]+u[k1+i+1]+u[k1+i-1]+u[k1+nxkk+i]+u[k1-nxkk+i])/4;
        }
    }
    for (h = 2; h < nykk-1; h = h+2) {
        k1 = levelk+nxkk*h;
        
        for (i = 1; i < nxkk-1; i=i+2) {
            u[k1+i] = (b[k1+i]+u[k1+i+1]+u[k1+i-1]+u[k1+nxkk+i]+u[k1-nxkk+i])/4;
        }
    }
}

/**
 * \fn void gsiteration_2color_3d(REAL *u, REAL *b, INT *level, INT k,
 *                           INT maxlevel, INT *nxk, INT *nyk, INT *nzk)
 * \brief 2 color G-S iteration of 3D problem
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void gsiteration_2color_3d(REAL *u,
                                  REAL *b,
                                  INT *level,
                                  INT k,
                                  INT maxlevel,
                                  INT *nxk,
                                  INT *nyk,
                                  INT *nzk)
{
    INT i,j,h;
    INT i0,j0,j1,j2,k0,k3,k4,k5,k6;
    const INT levelk = level[k];
    const INT nxkk = nxk[k];
    const INT nykk = nyk[k];
    const INT nzkk = nzk[k];
    const INT nxyk = nxkk*nykk;
    
    // red point of 2*i,2*j,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // red point of 2*i,2*j+1,2*h+1
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // points of (2*i+1,2*j+1,2h)
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // points of 2*i+1, 2*j, 2*h+1
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // Black points
    // 2*i,2*j,2*h+1
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i,2*j+1,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }	
    // 2*i+1,2*j,2*h
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i+1,2*j+1,2*h+1
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
}

/**
 * \fn void gsiteration3dpre(REAL *u, REAL *b, INT *level, INT k,
 *                        INT maxlevel, INT *nxk, INT *nyk, INT *nzk)
 * \brief colored G-S iteration of 3D problem, pre-smoothing
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void gsiteration3dpre(REAL *u,
                             REAL *b,
                             INT *level,
                             INT k,
                             INT maxlevel,
                             INT *nxk,
                             INT *nyk,
                             INT *nzk)
{
    INT i,j,h;
    INT i0,j0,j1,j2,k0,k3,k4,k5,k6;
    const INT levelk = level[k];
    const INT nxkk = nxk[k];
    const INT nykk = nyk[k];
    const INT nzkk = nzk[k];
    const INT nxyk = nxkk*nykk;
    
    // red point of 2*i,2*j,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i+1,2*j,2*h
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i,2*j+1,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // points of (2*i+1,2*j+1,2h)
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // Black points
    // 2*i,2*j,2*h+1
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // points of 2*i+1, 2*j, 2*h+1
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // red point of 2*i,2*j+1,2*h+1
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i+1,2*j+1,2*h+1
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0] = (b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
}

/**
 * \fn static void gsiteration3dpro (REAL *u, REAL *b, INT *level, INT k,
 *                                   INT maxlevel, INT *nxk, INT *nyk, INT *nzk)
 *
 * \brief Colored G-S iteration of 3D problem, pro-smoothing
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level
 * \param k            Level k
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/201
 */
static void gsiteration3dpro (REAL *u,
                              REAL *b,
                              INT *level,
                              INT k,
                              INT maxlevel,
                              INT *nxk,
                              INT *nyk,
                              INT *nzk)
{
    INT i,j,h;
    INT i0,j0,j1,j2,k0,k3,k4,k5,k6;
    const int levelk = level[k];
    const int nxkk = nxk[k];
    const int nykk = nyk[k];
    const int nzkk = nzk[k];
    const int nxyk = nxkk*nykk;
    
    // 2*i+1,2*j+1,2*h+1
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // 2*i,2*j+1,2*h+1
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // points of 2*i+1, 2*j, 2*h+1
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // 2*i,2*j,2*h+1
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 1; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // points of (2*i+1,2*j+1,2h)
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i,2*j+1,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 1; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    // 2*i+1,2*j,2*h
    for (i = 1; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0-1]+u[k0+1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
    
    // red point of 2*i,2*j,2*h
    for (i = 2; i < nzkk-1; i = i+2) {
        i0 = levelk+i*nxkk*nykk;
        for (j = 2; j < nykk-1; j = j+2) {
            j0 = i0+j*nxkk;
            j1 = i0+(j+1)*nxkk;
            j2 = i0+(j-1)*nxkk;
            for (h = 2; h < nxkk-1; h = h+2) {
                k0 = j0+h;
                k3 = j1+h;
                k4 = j2+h;
                k5 = k0+nxyk;
                k6 = k0-nxyk;
                u[k0]=(b[k0]+u[k0+1]+u[k0-1]+u[k3]+u[k4]+u[k5]+u[k6])/6;
            }
        }
    }
}

/**
 * \fn void multigriditeration1d(REAL *u, REAL *b, INT *level,
 *                               INT startlevel, INT maxlevel)
 * \brief V cycle starting from level k, where 0 is the finest level
 *        of 1D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param startlevel   Starting level of V-cycle
 * \param maxlevel     maxlevel of multigrids
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void multigriditeration1d(REAL *u,
                                 REAL *b,
                                 INT *level,
						         INT startlevel,
                                 INT maxlevel)
{
    INT n,i,j,k,levelk,levelk1;
    
    REAL *r = (REAL *)malloc(level[maxlevel+1]*sizeof(REAL));

    // forward sweep 
    for (k = startlevel; k < maxlevel-1; k++) {
        // initial
		levelk = level[k]; levelk1 = level[k+1];
		n = level[k+1] - level[k];
		for (i = 0; i < n; i++) {
			r[i+levelk] = 0.0;
		}
		if (k>startlevel) {
			for (i = 0; i < (level[k+1]-level[k]);i++) {
				u[level[k]+i] = 0.0;
			}
		}
        // pre-smoothing, G-S as smoother
        for (i = 0; i < 3; i++) {
            for (j = 1; j < n-1; j++) {
                u[levelk+j] = (b[levelk+j]+u[levelk-1+j]+u[levelk+1+j])/2;          
            }
   		}
		compute_r_1d(u, b, r, k, level);
		computenorm(r, level, k);

        // restriction on coarser grids
		n = level[k+2]-level[k+1];
        for (j = 1; j < n-1; j++) {
			b[levelk1+j] = (2*r[levelk+2*j]+r[levelk+2*j-1]+r[levelk+2*j+1]);
        }   
    } 
    
    // coarsest grid
    if (k==maxlevel-1) {
        i = level[maxlevel-1];
        u[i+1] = b[i+1]/2;
    }

    // back sweep
    for (k = maxlevel-1; k > startlevel; k--) {
        n = level[k+1] - level[k];
        // interpolation to finer grids
        for (i = 0; i < n-1; i++) {
			u[level[k-1]+2*i] += u[level[k]+i];
            u[level[k-1]+2*i+1] += (u[level[k]+i]+u[level[k]+i+1])/2;
        }
        // post-smoothing, G-S as smoother
        n = level[k] - level[k-1];
        for (i = 0; i < 3; i++) {
            for (j = n-2; j > 0; j--) {
                u[level[k-1]+j] = (b[level[k-1]+j]+u[level[k-1]-1+j]+u[level[k-1]+1+j])/2;    
            }
        }
    }
    free(r);
}

/**
 * \fn void multigriditeration2d(REAL *u, REAL *b, INT *level,
 *                               INT startlevel, INT maxlevel,
 *                               INT *nxk, INT *nyk)
 * \brief V cycle starting from level k, where 0 is the finest level
 *        of 2D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param startlevel   Starting level of V-cycle
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Pointer to the number of grids of x direction in level k
 * \param nyk          Pointer to the number of grids of y direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void multigriditeration2d(REAL *u,
                                 REAL *b,
                                 INT *level,
						         INT startlevel,
                                 INT maxlevel,
                                 INT *nxk,
                                 INT *nyk)
{
    INT   i,k,i1;
	INT  *prosm = (INT *)malloc(maxlevel*sizeof(INT));
	INT  *presm = (INT *)malloc(maxlevel*sizeof(INT));
    REAL *r = (REAL *)malloc(level[maxlevel+1]*sizeof(REAL));

	for (i = 0; i < maxlevel; i++) {
		prosm[i] = 3;
		presm[i] = 3;
	}

	prosm[0] = 1;
	presm[0] = 1;
	prosm[1] = 2;
	presm[1] = 2;

	k=startlevel;
        
    // pre-smoothing, GS as smoother
	// initial some vectors
	for (i = 0; i < (level[k+1]-level[k]); i++) {
		r[level[k]+i] = 0.0;
	}
	if (k>startlevel) {
		for (i = 0; i < (level[k+1]-level[k]); i++) u[level[k]+i] = 0.0;
	}

	// Gauss-Seidel 2 colors
	for (i = 0; i < presm[k]; i++) {
		gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);
	}           
	compute_r_2d(u, b, r, k, level, nxk, nyk);
    // restriction on coarser grids
    coarsergrid7pointrestriction2d(b, r, level, k, nxk, nyk);

	// forward sweep
    for (k = startlevel+1; k < maxlevel-1; k++) {
        
        // pre-smoothing, GS as smoother
		// initial some vectors
		for (i = 0; i < (level[k+1]-level[k]); i++) {
			r[level[k]+i] = 0.0;
			u[level[k]+i] = 0.0;
		}
		// Gauss-Seidel 2 colors
		for (i = 0; i < presm[k]; i++) {
			gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);
		}           
		compute_r_2d(u, b, r, k, level, nxk, nyk);
        // restriction on coarser grids
        coarsergrid7pointrestriction2d(b, r, level, k, nxk, nyk);
	}

    // coarsest grid
    if (k==maxlevel-1) {
        u[level[k]+4] = b[level[k]+4]/4;
		level[k+1] = level[k]+nxk[k]*nyk[k];
    }

    // back sweep
        
    // interpolation on finer grids
    for (k = maxlevel-1; k > startlevel; k--) {
		finergridinterpolation2d(u, level, k, nxk, nyk);
        // post-smoothing
		k = k-1;
		for (i1 = 0; i1 <prosm[k]; i1++) {
			// Gauss-Seidel 2 colors
			gsiteration_2color_2d(u, b, level, k, maxlevel, nxk, nyk);	
		}
		compute_r_2d(u, b, r, k, level, nxk, nyk);
		k = k+1;
    }
    
    free(r);
    free(prosm);
    free(presm);
}

/**
 * \fn void multigriditeration2d(REAL *u, REAL *b, INT *level,
 *                               INT startlevel, INT maxlevel,
 *                               INT *nxk, INT *nyk, INT *nzk)
 * \brief V cycle starting from level k, where 0 is the finest level
 *        of 3D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param startlevel   Starting level of V-cycle
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Pointer to the number of grids of x direction in level k
 * \param nyk          Pointer to the number of grids of y direction in level k
 * \param nzk          Pointer to the number of grids of z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void multigriditeration3d(REAL *u,
		                         REAL *b,
                                 INT *level,
		                         INT startlevel,
                                 INT maxlevel,
                                 INT *nxk,
                                 INT *nyk,
                                 INT *nzk)
{
    INT i,k;
    REAL *r;
    INT *presmoothtime,*prosmoothtime;

	presmoothtime = (INT *)malloc((maxlevel+1)*sizeof(INT));
	prosmoothtime = (INT *)malloc((maxlevel+1)*sizeof(INT));
	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));

	// set times of post and pre smoothing
	for (k = 0; k < maxlevel; k++) {
		presmoothtime[k] = 3;
	}
	presmoothtime[0] = 1;
	presmoothtime[1] = 2;
	presmoothtime[2] = 3;
	for (k = 0; k < maxlevel; k++) {
		prosmoothtime[k] = 3;
	}
	prosmoothtime[0] = 1;
	prosmoothtime[1] = 2;
	prosmoothtime[2] = 3;

#ifdef _OPENMP
		prosmoothtime[0] = 1;
#endif
    // forward sweep
	for	(i = 0; i < presmoothtime[0]; i++) {
#ifdef _OPENMP
		gsiteration_2color_3d(u, b, level, startlevel, maxlevel, nxk, nyk, nzk);
#else
		//gsiteration3dpre(u, b, level, startlevel, maxlevel, nxk, nyk, nzk);
		gsiteration_2color_3d(u, b, level, startlevel, maxlevel, nxk, nyk, nzk);
#endif
	}
	compute_r_3d(u, b, r, startlevel, level, nxk, nyk, nzk);
	// restriction on coarser grid
	coarsergrid7pointrestriction3d(b, r, level, startlevel, nxk, nyk, nzk);

	// coarser grids
    for (k = startlevel+1; k < maxlevel-1; k++) {

		// initial vectors
		for (i = 0; i < level[k+1]-level[k]; i++) {
			u[level[k]+i] = 0.0;
		}

        // pre-smoothing 
		for	(i = 0; i < presmoothtime[k]; i++) {
#ifdef _OPENMP
			gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
#else
			//gsiteration3dpre(u, b, level, startlevel, maxlevel, nxk, nyk, nzk);
			gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
#endif
		}
		compute_r_3d(u, b, r, k, level, nxk, nyk, nzk);
		// restriction on coarser grid
		coarsergrid7pointrestriction3d(b, r, level, k, nxk, nyk, nzk);
	}
	free(r);

	// coarsest level
	u[level[maxlevel-1]+13] = b[level[maxlevel-1]+13]/6;

	// back sweep
	for (k = maxlevel-1; k > 0; k--) {
		// interpolation from coarser grid		
		finergridinterpolation3d(u, level, k, nxk, nyk, nzk);
		// post smoothing
		k = k-1;
		for (i = 0; i < prosmoothtime[k]; i++) {
#ifdef _OPENMP
			gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
#else
			//gsiteration3dpro(u, b, level, k, maxlevel, nxk, nyk, nzk);
			gsiteration_2color_3d(u, b, level, k, maxlevel, nxk, nyk, nzk);
#endif
		}
		k = k+1;
	}
    
    free(presmoothtime);
    free(prosmoothtime);
}

/**
 * \fn void fullmultigrid_1d(REAL *u, REAL *b, INT *level, INT maxlevel,
 *                           INT nx)
 * \brief Full multigrid method of 1D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param maxlevel     maxlevel of multigrids
 * \param nx           Number of grids in x direction
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void fullmultigrid_1d(REAL *u,
                             REAL *b,
                             INT *level,
                             INT maxlevel,
                             INT nx)
{
    INT n,i,j,k,levelk;
    REAL *r;

	// initial
	r = (REAL *)malloc(level[maxlevel+1]*sizeof(REAL));
    n = level[1] - level[0];

    for (k = 0; k < maxlevel-1; k++) {

		// initial some vectors
		levelk = level[k];
		for (i = 0; i < (level[k+1]-level[k]); i++) {
			r[levelk+i] = 0.0;
		}
		if (k>0) {
			for (i = 0; i < (level[k+1]-level[k]); i++) {
				u[levelk+i] = 0.0;
			}
		}

		compute_r_1d(u, b, r, k, level);
        // restriction on coarser grids
		n = level[k+2]-level[k+1];
        for (j = 1; j < n-1; j++) {
			b[level[k+1]+j] = (2*r[level[k]+2*j]+r[level[k]+2*j-1]+r[level[k]+2*j+1]);      
        }   
    }

    // coarsest grid
    if (k==maxlevel-1) {
        i = level[maxlevel-1];
        u[i+1] = b[i+1]/2;
    }

    // FULL multigrid
    while (k>0) {
        n = level[k+1] - level[k];
        // interpolation to finer grids
        for (i = 0; i < n-1; i++) {
			u[level[k-1]+2*i] += u[level[k]+i];
            u[level[k-1]+2*i+1] += (u[level[k]+i]+u[level[k]+i+1])/2;
        }
		k = k-1;
        for(i=0;i<3;i++) {
			multigriditeration1d(u, b, level, k, maxlevel);
        }
    }
    
    free(r);
}

/**
 * \fn void fullmultigrid_2d(REAL *u, REAL *b, INT *level, INT maxlevel,
 *                           INT nxk, INT nyk)
 * \brief Full multigrid method of 2D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void fullmultigrid_2d(REAL *u,
                             REAL *b,
                             INT *level,
                             INT maxlevel,
                             INT *nxk,
                             INT *nyk)
{
    INT i,k;
    REAL *r;

	// initial
	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel],r,0.0);
    
    compute_r_2d(u, b, r, 0, level, nxk, nyk);
    // restriction on coarser grids
    coarsergrid7pointrestriction2d(b, r, level, 0, nxk, nyk);
    
    for (k = 1; k < maxlevel-1; k++) {
		// initial u of coarser grids
		if (k>0) {
			for (i = 0; i < (level[k+1]-level[k]); i++) {
				u[level[k]+i] = 0.0;
			}
		}
		compute_r_2d(u, b, r, k, level, nxk, nyk);
        // restriction on coarser grids
        coarsergrid7pointrestriction2d(b, r, level, k, nxk, nyk);
    }

    // coarsest grid
    if (k==maxlevel-1) {
        u[level[k]+4] = b[level[k]+4]/4;
		//level[k+1] = level[k]+nxk[k]*nyk[k];
    }

    // FULL multigrid
    while (k>0) {
        // interpolation from coarser grid
        finergridinterpolation2d(u, level, k, nxk, nyk);
		k = k-1;
        for(i=0;i<3;i++) {
			multigriditeration2d(u, b, level, k, maxlevel, nxk, nyk);
			compute_r_2d(u, b, r, k, level, nxk, nyk);
			computenorm(r, level, k);
        }
    }
    
    free(r);
}

/**
 * \fn static void fullmultigrid_3d (REAL *u, REAL *b, INT *level, INT maxlevel,
 *                                   INT *nxk, INT *nyk, INT *nzk)
 * \brief Full multigrid method of 3D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param maxlevel     maxlevel of multigrids
 * \param nxk          Number of grids in x direction in level k
 * \param nyk          Number of grids in y direction in level k
 * \param nzk          Number of grids in z direction in level k
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static void fullmultigrid_3d (REAL *u,
                              REAL *b,
                              INT *level,
                              INT maxlevel,
                              INT *nxk,
                              INT *nyk,
                              INT *nzk)
{
    INT i,k,levelk1;
    REAL *r;

	// initial
	r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));

    fasp_array_set(level[maxlevel], r, 0.0);

    for (k = 0; k < maxlevel-1; k++) {
		compute_r_3d(u, b, r, k, level, nxk, nyk, nzk);
        coarsergrid7pointrestriction3d(b, r, level, k, nxk, nyk, nzk);
    }

    // coarsest grid
	u[level[maxlevel-1]+13] = b[level[maxlevel-1]+13]/6;

    // FULL multigrid
    for (k = maxlevel-1; k > 1; k--) {
        // interpolation from coarser grid
		levelk1 = level[k-1];
		for (i=0;i<level[k]-level[k-1]+1;i++) u[levelk1+i] = 0.0;
        finergridinterpolation3d(u, level, k, nxk, nyk, nzk);
        for(i = 0;i<2;i++) {
			multigriditeration3d(u, b, level, k-1, maxlevel, nxk, nyk, nzk);
			compute_r_3d(u, b, r, k-1, level, nxk, nyk, nzk);
        }
    }
    for (i=0;i<level[1]-level[0]+1;i++) u[i] = 0.0;
    finergridinterpolation3d(u, level, 1, nxk, nyk, nzk);
	multigriditeration3d(u, b, level, 0, maxlevel, nxk, nyk, nzk);
    free(r);
}

/**
 * \fn static INT pcg_1d (REAL *u, REAL *b, INT *level, INT maxlevel,
 *                        INT nx, INT rtol, INT maxiteration, const SHORT prtlvl)
 *
 * \brief Preconditioned CG method of 1D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param maxlevel     Max level of multigrids
 * \param nx           Number of grids in x direction
 * \param rtol         Relative Tolerance judging convergence
 * \param maxiteration Number of maximum iteration number of CG method
 * \param prtlvl       Print level of iterative method
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static INT pcg_1d (REAL *u,
                   REAL *b,
                   INT *level,
		           INT maxlevel,
                   INT nx,
                   REAL rtol,
		           INT maxiteration,
				   const SHORT prtlvl)
{
    INT k;
    REAL *p, *r, *z, *q;
    REAL rh0, rh1, rh2, alfa, beta, normb, normr, resid, normr1, factor;

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    p = (REAL *)malloc(level[1]*sizeof(REAL));
    r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	z = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	q = (REAL *)malloc(level[1]*sizeof(REAL));

    k = 0;

    // initial residue and other vector
    fasp_array_set(level[maxlevel], z, 0.0);
    fasp_array_set(level[maxlevel], r, 0.0);
    fasp_array_set(level[1], p, 0.0);
    fasp_array_set(level[1], q, 0.0);    

	compute_r_1d(u, b, r, 0, level);
	normr = computenorm(r, level, 0);
	normb = computenorm(b, level, 0);
	normr1 = normr;
	if (normb==0.0) normb=1.0;
	if ((resid = normr / normb) <= rtol) goto FINISHED; 

	if ( prtlvl > PRINT_SOME ){
		printf("-----------------------------------------------------------\n");
		printf("It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
		printf("-----------------------------------------------------------\n");
	}

    multigriditeration1d(z, r, level, 0, maxlevel);
    rh0 = innerproductxy(r, z, level, 0);
    xequaly(p, z, level, 0);

    while (k < maxiteration) {
		// init z
		fasp_array_set(level[1], z, 0.0);

        // calculating alpha
        xequalay_1d(q, p, level, 0);
		rh2 = innerproductxy(q, p, level, 0);
		alfa = rh0/rh2;

        // update vector u, r
		xequalypcz(u, u, alfa, p, level, 0);
		xequalypcz(r, r, (-alfa), q, level, 0);
		normr = computenorm(r, level, 0);
		resid = normr / normb;
		factor = normr / normr1;
		if ( prtlvl > PRINT_SOME ){
			printf("%6d | %13.6e   | %13.6e  | %10.4f\n",k+1,resid,normr,factor);
		}
		normr1 = normr;
		if ((resid) <= rtol) break;

        // update z and beta
        multigriditeration1d(z, r, level, 0, maxlevel);
		//xequaly(z, r, level, 0);	
        rh1 = innerproductxy(r, z, level, 0);
        beta = rh1 / rh0;

        // update p
        xequalypcz(p, z, beta, p, level, 0);

		rh0 = rh1;
		k++;
	}

	if ( prtlvl > PRINT_NONE ){
		if (k >= maxiteration) {
			printf("### WARNING: V-cycle failed to converge.\n");
		}
		else {
			printf("Num of Iter's: %d, Relative Residual = %e.\n", k+1, normr);
		}
	} 

FINISHED:    
    free(r);
    free(q);
    free(p);
	free(z);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif
    
	return k+1;
}

/**
 * \fn static INT pcg_2d (REAL *u, REAL *b, INT *level, INT maxlevel, INT *nxk, 
 *                        INT *nyk, INT rtol, INT maxiteration, const SHORT prtlvl)
 *
 * \brief Preconditioned CG method of 2D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param maxlevel     Max level of multigrids
 * \param nxk          Number of grids in x direction
 * \param nyk          Number of grids in y direction
 * \param rtol         Relative Tolerance judging convergence
 * \param maxiteration Number of maximum iteration number of CG method
 * \param prtlvl       Print level of iterative method
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static INT pcg_2d(REAL *u,
                  REAL *b,
                  INT *level,
		          INT maxlevel,
                  INT *nxk,
                  INT *nyk,
                  REAL rtol,
		          INT maxiteration,
		          const SHORT prtlvl)
{
    INT k = 0;
    REAL *p, *r, *z, *q;
    REAL rh0, rh1, rh2, alfa, beta, normb, normr, resid, normr1, factor;

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    p = (REAL *)malloc(level[1]*sizeof(REAL));
    r = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	z = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	q = (REAL *)malloc(level[1]*sizeof(REAL));

    // initial residue and other vector
    fasp_array_set(level[maxlevel], z, 0.0);
    fasp_array_set(level[maxlevel], r, 0.0);
    fasp_array_set(level[1], p, 0.0);
    fasp_array_set(level[1], q, 0.0);    

	compute_r_2d(u, b, r, 0, level, nxk, nyk);
	normr = computenorm(r, level, 0);
	normb = computenorm(b, level, 0);
	normr1 = normr;
	if (normb==0.0) normb=1.0;
	if ((resid = normr / normb) <= rtol) goto FINISHED;

	if ( prtlvl > PRINT_SOME ){
		printf("-----------------------------------------------------------\n");
		printf("It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
		printf("-----------------------------------------------------------\n");
	}

    multigriditeration2d(z, r, level, 0, maxlevel, nxk, nyk);
    rh0 = innerproductxy(r, z, level, 0);
    xequaly(p, z, level, 0);

    while (k < maxiteration) {
		// init z
		fasp_array_set(level[1], z, 0.0);

        // calculating alpha
        xequalay_2d(q, p, level, 0, nxk, nyk);
		rh2 = innerproductxy(q, p, level, 0);
		alfa = rh0/rh2;

        // update vector u, r
		xequalypcz(u, u, alfa, p, level, 0);
		xequalypcz(r, r, (-alfa), q, level, 0);
		normr = computenorm(r, level, 0);
		resid = normr / normb;
		factor = normr / normr1;
		if ( prtlvl > PRINT_SOME ){
			printf("%6d | %13.6e   | %13.6e  | %10.4f\n",k+1,resid,normr,factor);
		}
		normr1 = normr;
		if (resid <= rtol) break;     
		
        // update z and beta
        multigriditeration2d(z, r, level, 0, maxlevel, nxk, nyk);
		//xequaly(z, r, level, 0);	
        rh1 = innerproductxy(r, z, level, 0);
        beta = rh1 / rh0;

        // update p
        xequalypcz(p, z, beta, p, level, 0);

		rh0 = rh1;
		k++;
	}

	if ( prtlvl > PRINT_NONE ){
		if (k >= maxiteration) {
			printf("### WARNING: V-cycle failed to converge.\n");
		}
		else {
			printf("Num of Iter's: %d, Relative Residual = %e.\n", k+1, normr);
		}
	} 

FINISHED:
    free(r);
    free(q);
    free(p);
	free(z);

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

	return k+1;
}

/**
 * \fn static INT pcg_3d (REAL *u, REAL *b, INT *level, INT maxlevel, INT *nxk,
 *                        INT *nyk, INT *nzk, INT rtol, INT maxiteration, 
 *                        const SHORT prtlvl)
 *
 * \brief Preconditioned CG method of 2D poisson equation
 *
 * \param u            Pointer to the vector of DOFs
 * \param b            Pointer to the right hand vector
 * \param level        Pointer to the start position of each level 
 * \param maxlevel     Max level of multigrids
 * \param nxk          Number of grids in x direction 
 * \param nyk          Number of grids in y direction
 * \param nzk          Number of grids in z direction
 * \param rtol         Relative Tolerance judging convergence
 * \param maxiteration Number of maximum iteration number of CG method
 * \param prtlvl       Print level of iterative method
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
static INT pcg_3d (REAL *u,
                   REAL *b,
                   INT *level,
		           INT maxlevel,
                   INT *nxk,
                   INT *nyk,
                   INT *nzk, 
                   REAL rtol,
		           INT maxiteration,
				   const SHORT prtlvl)
{
    INT i, k = 0, done = 0;
    REAL *p, *r, *z, *q;
    REAL rh0, rh1, rh2, alfa, beta, normb, normr, resid, normr1, factor;
	const INT level1 = level[1];
	const INT levelmax = level[maxlevel];

#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Start]\n", __FUNCTION__);
#endif

    p = (REAL *)malloc(level1*sizeof(REAL));
    r = (REAL *)malloc(levelmax*sizeof(REAL));
	z = (REAL *)malloc(levelmax*sizeof(REAL));
	q = (REAL *)malloc(level1*sizeof(REAL));

    // initial residue and other vector
    fasp_array_set(levelmax, z, 0.0);
    fasp_array_set(levelmax, r, 0.0);
    fasp_array_set(level[1], p, 0.0);
    fasp_array_set(level[1], q, 0.0);
	compute_r_3d(u, b, r, 0, level, nxk, nyk, nzk);
	normr = computenorm(r, level, 0);
	normb = computenorm(b, level, 0);
	normr1 = normr;
	if (normb==0.0) normb=1.0;
	if ((resid = normr / normb) <= rtol) goto FINISHED;

	if ( prtlvl > PRINT_SOME ){
		printf("-----------------------------------------------------------\n");
		printf("It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
		printf("-----------------------------------------------------------\n");
	}

	multigriditeration3d(z, r, level, 0, maxlevel, nxk, nyk, nzk);
    rh0 = innerproductxy(r, z, level, 0);
    xequaly(p, z, level, 0);

    while (!done && k < maxiteration) {
		// init z
        for (i = 0; i < level[1]; i++) z[i] = 0;

        // calculating alpha
        xequalay_3d(q, p, level, 0, nxk, nyk, nzk);
		rh2 = innerproductxy(q, p, level, 0);
		alfa = rh0/rh2;

        // update vector u, r
		xequalypcz(u, u, alfa, p, level, 0);
		xequalypcz(r, r, (-alfa), q, level, 0);
		normr = computenorm(r, level, 0);
		resid = normr / normb;
		factor = normr / normr1;
		if ( prtlvl > PRINT_SOME ){
			printf("%6d | %13.6e   | %13.6e  | %10.4f\n",k+1,resid,normr,factor);
		}
		normr1 = normr;
		if (resid <= rtol) break;

        // update z and beta
        multigriditeration3d(z, r, level, 0, maxlevel, nxk, nyk, nzk);
		//xequaly(z, r, level, 0);	
        rh1 = innerproductxy(r, z, level, 0);
        beta = rh1 / rh0;

        // update p
        xequalypcz(p, z, beta, p, level, 0);

		rh0 = rh1;
		k++;
	}

	if ( prtlvl > PRINT_NONE ){
		if (k >= maxiteration) {
			printf("### WARNING: V-cycle failed to converge.\n");
		}
		else {
			printf("Num of Iter's: %d, Relative Residual = %e.\n", k+1, normr);
		}
	} 

FINISHED:
    free(r);
    free(q);
    free(p);
	free(z);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

	return k+1;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
