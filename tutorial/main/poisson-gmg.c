/*! \file  poisson-gmg.c
 *
 *  \brief A test example for FASP: Using GMG method to solve the discrete 
 *         Poisson equation (five-point stencil) in 2D.
 *  
 *  \note  GMG test example for FASP: C version
 *
 *  Solving the Poisson Equation with GMG method.
 */

#include <time.h>
#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

const REAL pi = 3.14159265;

/**
 * \fn static REAL f2d (INT i, INT j, INT nx, INT ny)
 *
 * \brief The right-hand side function f in the Poisson equation, where
 *        f = sin(pi x)*sin(pi y)          
 *  
 * \param i      Index in x direction
 * \param j      Index in y direction
 * \param nx     Number of grids in x direction
 * \param ny     Number of grids in y direction
 *
 * \author Ziteng Wang
 * \date   07/01/2013
 */
static REAL f2d(INT i,
                INT j,
                INT nx,
                INT ny)
{
    return sin(pi *(((REAL) j)/((REAL) nx)))*sin(pi *(((REAL) i)/((REAL) nx)));
}

/**
 * \fn static REAL L2NormError2d (REAL *u, INT nx, INT ny)
 *
 * \brief Computing Discretization Error, where exact solution
 *        u = sin(pi x)*sin(pi y)/(2*pi*pi)
 *
 * \param u      Vector of DOFs
 * \param nx     Number of grids in x direction
 * \param ny     Number of grids in y direction
 *
 * \author Ziteng Wang
 * \date   07/01/2013
 */
static REAL L2NormError2d (REAL *u,
                           INT nx,
                           INT ny)
{
    const REAL h = 1.0/nx;
    REAL l2norm  = 0.0, uexact;
    
    INT i, j;
    for ( i = 1; i < ny; i++ ) {
        for ( j = 1; j < nx; j++ ) {
            uexact  = sin(pi*i*h)*sin(pi*j*h)/(pi*pi*2.0);
            l2norm += pow((u[i*(nx+1)+j] - uexact), 2);
        }
    }
    l2norm = sqrt(l2norm*h*h);
    
    return l2norm;
}

/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief  An example of GMG V-cycle for the Poisson equation.
 *
 * \author Ziteng Wang
 * \date   07/02/2013
 *
 * \note   Number of grids of nx, ny should be all equal to 2^maxlevel.
 * 
 * Modified by Chensong Zhang on 07/02/2013: Format the output.
 */
int main (int argc, const char *argv[])
{
    const REAL rtol = 1.0e-8;
    const INT  print_level = 0;

    INT        maxlevel, iter = 0, nx, ny;
    INT        i, j;
    REAL       h, error0, *u, *b;
    REAL       GMG_start, GMG_end;
   
    printf("=================================================\n");
    printf("            FASP GMG V-cycle in 2D \n");
    printf("=================================================\n");
    printf("  Level  Iter    CPU time      L2 Error\n");
    
    // Testing different Number of DOFs nx = ny = 2^2 to 2^10
    for ( maxlevel = 2; maxlevel <= 10; maxlevel++ ) {
    
        // Step 1: Set size of the grid
        nx = ny = pow(2.0, maxlevel); h = 1.0/((REAL) nx);

        // Step 2: Solving Poisson equation with GMG solver
        // Set initial guess 0, right hand side as given f                  
        u = (REAL *)malloc((nx+1)*(ny+1)*sizeof(REAL));
        fasp_array_set((nx+1)*(ny+1), u, 0.0);
        
        b = (REAL *)malloc((nx+1)*(ny+1)*sizeof(REAL));
        for ( i = 0; i <= nx; i++ ) {
            for ( j = 0; j <= ny; j++ ) {
                b[j*(nx+1)+i] = h*h*f2d(i, j, nx, ny);
            }
        }
        
        // Step 3: Solve equation with V-cycle
        fasp_gettime(&GMG_start);
        iter = fasp_poisson_gmg_2D(u, b, nx, ny, maxlevel, rtol, print_level);
        fasp_gettime(&GMG_end);

        error0 = L2NormError2d(u, nx, ny);
    
        printf("%5d  %5d %12.6f %16.5e\n", maxlevel, iter, GMG_end-GMG_start, error0);
      
        // Clean up memory  
        free(u);
        free(b);
    }
    
    return FASP_SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
