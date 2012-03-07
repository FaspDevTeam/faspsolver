/**
 *       @file  setup_poisson.h
 *      @brief  setup p1 fem for poisson
 *
 *     @author  Feiteng Huang
 *
 *   @internal
 *     created  02/21/2012
 *
 * =====================================================================================
 */

#ifndef _SETUP_POISSON_H_
#define _SETUP_POISSON_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

extern double f(double x, double y);
extern double u(double x, double y);

int setup_poisson (dCSRmat *ptr_A, dvector *ptr_b, int levelNum, const char *meshIn, 
                   const char *meshOut, int mo, const char *assemble_option, int num_qp_rhs, int num_qp_mat);

#endif
