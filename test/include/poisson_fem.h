/*! \file poisson_fem.h
 *  \brief Main header file for the 2D Finite Element Method
 */

#ifndef _SETUP_POISSON_H_
#define _SETUP_POISSON_H_

extern double f(double x, double y);
extern double u(double x, double y);

int setup_poisson (dCSRmat *ptr_A, 
                   dvector *ptr_b, 
                   int levelNum, 
                   const char *meshIn, 
                   const char *meshOut, 
                   int mo, 
                   const char *assemble_option, 
                   int num_qp_rhs, 
                   int num_qp_mat);

#endif
