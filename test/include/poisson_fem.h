/*! \file poisson_fem.h
 *  \brief Main header file for the 2D Finite Element Method
 */

#ifndef _SETUP_POISSON_H_
#define _SETUP_POISSON_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "misc.h"
#include "mesh.h"
#include "fasp.h"
#include "fasp_functs.h"

extern double f(double *p);
extern double u(double *p);

int setup_poisson (dCSRmat *A, dvector *b, Mesh *mesh, Mesh_aux *mesh_aux, param_test *pt, dvector *ptr_uh, ivector *dof);

double get_l2_error_poisson(ddenmat *node,
														idenmat *elem,
														dvector *uh,
														int num_qp);


#endif
