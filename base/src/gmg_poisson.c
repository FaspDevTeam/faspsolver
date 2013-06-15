/*! \file  fasp_poisson_gmg.c
 *  \brief GMG method as an iterative solver for Poisson Problem
 */

#include <time.h>
#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "gmg_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_poisson_gmg_1D (REAL *u, REAL *b, INT nx, INT maxlevel,
 *                               REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 1D equation by Geometric Multigrid Method
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param maxlevel  Maximum levels of the multigrid
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_gmg_1D (REAL *u,
                          REAL *b,
                          INT nx,
                          INT maxlevel,
                          REAL rtol)
{
    const REAL atol = 1.0E-15;
    const INT  max_itr_num = 100;
    
    REAL      *u0, *r0, *b0;
    REAL       norm_r, norm_r0, error;
    INT        i, *level, count;
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = nx+1;
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(level[i]-level[i-1]+1)/2;
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0, r0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    
    fasp_array_cp(nx, u, u0);
    fasp_array_cp(nx, b, b0);
    
    // compute initial l2 norm of residue
    fasp_array_set(level[1], r0, 0.0);
    compute_r_1d(u0, b0, r0, 0, level);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // GMG solver of V-cycle
    count = 0;
    while (count < max_itr_num) {
        count++;
        multigriditeration1d(u0, b0, level, 0, maxlevel);
        compute_r_1d(u0, b0, r0, 0, level);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
        if (error < rtol || norm_r < atol) break;
    }
    
    if (count >= max_itr_num) {
        printf("### WARNING: V-cycle failed to converge.\n");
    }
    else {
        printf("Num of V-cycle's: %d, Relative Residual = %e.\n", count, error);
    }
    
    // Update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_gmg_2D (REAL *u, REAL *b, INT nx, INT nx,
 *                               INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 2D equation by Geometric Multigrid Method
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param ny        Number of grids in y direction
 * \param maxlevel  Maximum levels of the multigrid
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_gmg_2D (REAL *u,
                          REAL *b,
                          INT nx,
                          INT ny,
                          INT maxlevel,
                          REAL rtol)
{
    const REAL atol = 1.0E-15;
    const INT  max_itr_num = 100;

    REAL *u0,*r0,*b0;
    REAL norm_r,norm_r0,error;
    INT i, k, count, *nxk, *nyk, *level;
    
    // set nxk, nyk
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
	nyk = (INT *)malloc(maxlevel*sizeof(INT));
	nxk[0] = nx+1; nyk[0] = ny+1;
    for(k=1;k<maxlevel;k++){
		nxk[k] = (int) (nxk[k-1]+1)/2;
		nyk[k] = (int) (nyk[k-1]+1)/2;
	}
    
    // set level
    level = (INT *)malloc((maxlevel+1)*sizeof(REAL));
    level[0] = 0; level[1] = (nx+1)*(ny+1);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2.0,i)+1)*(ny/pow(2.0,i)+1);
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(level[1], u, u0);
    fasp_array_cp(level[1], b, b0);
    
    // compute initial l2 norm of residue
    compute_r_2d(u0, b0, r0, 0, level, nxk, nyk);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // GMG solver of V-cycle
    count = 0;
    while (count < max_itr_num) {
        count++;
        multigriditeration2d(u0, b0, level, 0, maxlevel, nxk, nyk);
        compute_r_2d(u0, b0, r0, 0, level, nxk, nyk);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
        if (error < rtol || norm_r < atol) break;
    }
    
    if (count >= max_itr_num) {
        printf("### WARNING: V-cycle failed to converge.\n");
    }
    else {
        printf("Num of V-cycle's: %d, Relative Residual = %e.\n", count, error);
    }
    
	// update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(nxk);
    free(nyk);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_gmg_3D (REAL *u, REAL *b, INT nx, INT nx, INT nz,
 *                               INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 3D equation by Geometric Multigrid Method
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param ny        Number of grids in y direction
 * \param nz        Number of grids in z direction
 * \param maxlevel  Maximum levels of the multigrid
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_gmg_3D (REAL *u,
                          REAL *b,
                          INT nx,
                          INT ny,
                          INT nz,
                          INT maxlevel,
                          REAL rtol)
{
    const REAL atol = 1.0E-15;
    const INT  max_itr_num = 100;

    REAL *u0,*r0,*b0;
    REAL norm_r,norm_r0,error;
    int i, k, count, *nxk, *nyk, *nzk, *level;
    
    // set nxk, nyk, nzk
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
    nyk = (INT *)malloc(maxlevel*sizeof(INT));
    nzk = (INT *)malloc(maxlevel*sizeof(INT));
	nxk[0] = nx+1; nyk[0] = ny+1; nzk[0] = nz+1;
    for(k=1;k<maxlevel;k++){
		nxk[k] = (int) (nxk[k-1]+1)/2;
		nyk[k] = (int) (nyk[k-1]+1)/2;
        nzk[k] = (int) (nyk[k-1]+1)/2;
	}
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = (nx+1)*(ny+1)*(nz+1);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2.0,i)+1)*(ny/pow(2.0,i)+1)*(nz/pow(2.0,i)+1);
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0, r0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(level[1], u, u0);
    fasp_array_cp(level[1], b, b0);
    
    // compute initial l2 norm of residue
    compute_r_3d(u0, b0, r0, 0, level, nxk, nyk, nzk);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // GMG solver of V-cycle
    count = 0;
    while (count < max_itr_num) {
        count++;
        multigriditeration3d(u0, b0, level, 0, maxlevel, nxk, nyk, nzk);
        compute_r_3d(u0, b0, r0, 0, level, nxk, nyk, nzk);
        norm_r = computenorm(r0, level, 0);
        error = norm_r / norm_r0;
        if (error < rtol || norm_r < atol) break;
    }
    
    if (count >= max_itr_num) {
        printf("### WARNING: V-cycle failed to converge.\n");
    }
    else {
        printf("Num of V-cycle's: %d, Relative Residual = %e.\n", count, error);
    }
    
	// update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(nxk);
    free(nyk);
    free(nzk);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_fgmg_1D (REAL *u, REAL *b, INT nx,
 *                                INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 1D equation by Geometric Multigrid Method
 *        (Full Multigrid)
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param maxlevel  Maximum levels of the multigrid 
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_fgmg_1D (REAL *u,
                           REAL *b,
                           INT nx,
                           INT maxlevel,
                           REAL rtol)
{
    const REAL atol = 1.0E-15;
    REAL *u0,*r0,*b0;
    REAL norm_r0;
    int i, *level;
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = nx+1;
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(level[i]-level[i-1]+1)/2;
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0, r0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(nx, u, u0);
    fasp_array_cp(nx, b, b0);    
    
    // compute initial l2 norm of residue
    fasp_array_set(level[1], r0, 0.0);
    compute_r_1d(u0, b0, r0, 0, level);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    //  Full GMG solver 
    fullmultigrid_1d(u0, b0, level, maxlevel, nx);
    
    // update u
    fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_fgmg_2D (REAL *u, REAL *b, INT nx, INT ny,
 *                                INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 2D equation by Geometric Multigrid Method
 *        (Full Multigrid)
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param ny        Number of grids in Y direction
 * \param maxlevel  Maximum levels of the multigrid 
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_fgmg_2D (REAL *u,
                           REAL *b,
                           INT nx,
                           INT ny,
                           INT maxlevel,
                           REAL rtol)
{
    const REAL atol = 1.0E-15;
    REAL *u0,*r0,*b0;
    REAL norm_r0;
    int i, k, *nxk, *nyk, *level;
    
    // set nxk, nyk
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
	nyk = (INT *)malloc(maxlevel*sizeof(INT));
    
	nxk[0] = nx+1; nyk[0] = ny+1;
    for(k=1;k<maxlevel;k++){
		nxk[k] = (int) (nxk[k-1]+1)/2;
		nyk[k] = (int) (nyk[k-1]+1)/2;
	}
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = (nx+1)*(ny+1);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2.0,i)+1)*(ny/pow(2.0,i)+1);
    }
	level[maxlevel+1] = level[maxlevel] + 1;
    
    // set u0, b0, r0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(level[1], u, u0);
    fasp_array_cp(level[1], b, b0);    
    
    // compute initial l2 norm of residue
    fasp_array_set(level[1], r0, 0.0);
    compute_r_2d(u0, b0, r0, 0, level, nxk, nyk);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // FMG solver
    fullmultigrid_2d(u0, b0, level, maxlevel, nxk, nyk);
    
	// update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(nxk);
    free(nyk);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_fgmg_3D (REAL *u, REAL *b, INT nx, INT ny, INT nz,
 *                                INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 3D equation by Geometric Multigrid Method
 *        (Full Multigrid)
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param ny        NUmber of grids in y direction
 * \param nz        NUmber of grids in z direction
 * \param maxlevel  Maximum levels of the multigrid 
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_fgmg_3D (REAL *u,
                           REAL *b,
                           INT nx,
                           INT ny,
                           INT nz,
                           INT maxlevel,
                           REAL rtol)
{
    const REAL atol = 1.0E-15;
    REAL *u0,*r0,*b0;
    REAL norm_r0;
    int i, k, *nxk, *nyk, *nzk, *level;
    
    // set nxk, nyk, nzk
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
    nyk = (INT *)malloc(maxlevel*sizeof(INT));
    nzk = (INT *)malloc(maxlevel*sizeof(INT));
    
	nxk[0] = nx+1; nyk[0] = ny+1; nzk[0] = nz+1;
    for(k=1;k<maxlevel;k++){
		nxk[k] = (int) (nxk[k-1]+1)/2;
		nyk[k] = (int) (nyk[k-1]+1)/2;
        nzk[k] = (int) (nyk[k-1]+1)/2;     
	}
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = (nx+1)*(ny+1)*(nz+1);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2.0,i)+1)*(ny/pow(2.0,i)+1)*(nz/pow(2.0,i)+1);
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0, r0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(level[1], u, u0);
    fasp_array_cp(level[1], b, b0);
    
    // compute initial l2 norm of residue
    compute_r_3d(u0, b0, r0, 0, level, nxk, nyk, nzk);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // FMG
    fullmultigrid_3d(u0, b0, level, maxlevel, nxk, nyk, nzk);
    
	// update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(nxk);
    free(nyk);
    free(nzk);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_pcg_gmg_1D (REAL *u, REAL *b, INT nx,
 *                                   INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 1D equation by Geometric Multigrid Method
 *        (GMG preconditioned Conjugate Gradient method)
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param maxlevel  Maximum levels of the multigrid 
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_pcg_gmg_1D (REAL *u,
                              REAL *b,
                              INT nx,
                              INT maxlevel,
                              REAL rtol)
{
    const REAL atol = 1.0E-15;
    const INT  max_itr_num = 100;

    REAL *u0,*r0,*b0;
    REAL norm_r0;
    int i, *level;
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = nx+1;
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(level[i]-level[i-1]+1)/2;
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(nx, u, u0);
    fasp_array_cp(nx, b, b0);    
    
    // compute initial l2 norm of residue
    fasp_array_set(level[1], r0, 0.0);
    compute_r_1d(u, b, r0, 0, level);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;

    // Preconditioned CG method
    pcg_1d(u0, b0, level, maxlevel, nx, rtol, max_itr_num);
    
    // Update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}    

/**
 * \fn void fasp_poisson_pcg_gmg_2D (REAL *u, REAL *b, INT nx, INT ny, 
 *                                   INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 2D equation by Geometric Multigrid Method
 *        (GMG preconditioned Conjugate Gradient method)
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param ny        Number of grids in y direction
 * \param maxlevel  Maximum levels of the multigrid 
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_pcg_gmg_2D (REAL *u,
                              REAL *b,
                              INT nx,
                              INT ny,
                              INT maxlevel,
                              REAL rtol)
{
    const REAL atol = 1.0E-15;
    const INT  max_itr_num = 100;

    REAL *u0,*r0,*b0;
    REAL norm_r0;
    int i, k, *nxk, *nyk, *level;
    
    // set nxk, nyk
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
	nyk = (INT *)malloc(maxlevel*sizeof(INT));
    
	nxk[0] = nx+1; nyk[0] = ny+1; 
    for(k=1;k<maxlevel;k++){
		nxk[k] = (int) (nxk[k-1]+1)/2;
		nyk[k] = (int) (nyk[k-1]+1)/2;
	}
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = (nx+1)*(ny+1);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2.0,i)+1)*(ny/pow(2.0,i)+1);
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0, r0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(level[1], u, u0);
    fasp_array_cp(level[1], b, b0);    
    
    // compute initial l2 norm of residue
    fasp_array_set(level[1], r0, 0.0);
    compute_r_2d(u0, b0, r0, 0, level, nxk, nyk);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // Preconditioned CG method
    pcg_2d(u0, b0, level, maxlevel, nxk, nyk, rtol, max_itr_num);
    
	// update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(nxk);
    free(nyk);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/**
 * \fn void fasp_poisson_pcg_gmg_3D (REAL *u, REAL *b, INT nx,
 *                                   INT maxlevel, REAL rtol)
 *
 * \brief Solve Ax=b of Poisson 3D equation by Geometric Multigrid Method
 *        (GMG preconditioned Conjugate Gradient method)
 *
 * \param u         Pointer to the vector of dofs
 * \param b         Pointer to the vector of right hand side
 * \param nx        Number of grids in x direction
 * \param ny        Number of grids in y direction
 * \param nz        Number of grids in z direction
 * \param maxlevel  Maximum levels of the multigrid 
 * \param rtol      Relative tolerance to judge convergence
 *
 * \author Ziteng Wang
 * \date   06/07/2013
 */
void fasp_poisson_pcg_gmg_3D (REAL *u,
                              REAL *b,
                              INT nx,
                              INT ny,
                              INT nz,
                              INT maxlevel,
                              REAL rtol)
{
    const REAL atol = 1.0E-15;
    const INT  max_itr_num = 100;

    REAL *u0,*r0,*b0;
    REAL norm_r0;
    int i, k, *nxk, *nyk, *nzk, *level;
    
    // set nxk, nyk, nzk
    nxk = (INT *)malloc(maxlevel*sizeof(INT));
    nyk = (INT *)malloc(maxlevel*sizeof(INT));
    nzk = (INT *)malloc(maxlevel*sizeof(INT));
    
	nxk[0] = nx+1; nyk[0] = ny+1; nzk[0] = nz+1;
    for(k=1;k<maxlevel;k++){
		nxk[k] = (int) (nxk[k-1]+1)/2;
		nyk[k] = (int) (nyk[k-1]+1)/2;
        nzk[k] = (int) (nyk[k-1]+1)/2;     
	}
    
    // set level
    level = (INT *)malloc((maxlevel+2)*sizeof(REAL));
    level[0] = 0; level[1] = (nx+1)*(ny+1)*(nz+1);
    for (i = 1; i < maxlevel; i++) {
        level[i+1] = level[i]+(nx/pow(2.0,i)+1)*(ny/pow(2.0,i)+1)*(nz/pow(2.0,i)+1);
    }
	level[maxlevel+1] = level[maxlevel]+1;
    
    // set u0, b0
    u0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    b0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
	r0 = (REAL *)malloc(level[maxlevel]*sizeof(REAL));
    fasp_array_set(level[maxlevel], u0, 0.0);
    fasp_array_set(level[maxlevel], b0, 0.0);
	fasp_array_set(level[maxlevel], r0, 0.0);
    fasp_array_cp(level[1], u, u0);
    fasp_array_cp(level[1], b, b0);
    
    // compute initial l2 norm of residue
    compute_r_3d(u0, b0, r0, 0, level, nxk, nyk, nzk);
    norm_r0 = computenorm(r0, level, 0);
    if (norm_r0 < atol) goto FINISHED;
    
    // Preconditioned CG method
    pcg_3d(u0, b0, level, maxlevel, nxk, nyk, nzk, rtol, max_itr_num);
    
	// update u
	fasp_array_cp(level[1], u0, u);
    
FINISHED:
    free(level);
    free(nxk);
    free(nyk);
    free(nzk);
    free(r0);
    free(u0);
    free(b0);
    
    return;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
