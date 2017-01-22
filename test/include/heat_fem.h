/*! \file heat_fem.h
 *  \brief Main header file for the 2D Finite Element Method
 */

#ifndef _HEAT_FEM_H_
#define _HEAT_FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "misc.h"
#include "mesh.h"
#include "assemble.h"
#include "fasp.h"
#include "fasp_functs.h"

/** 
 * \struct bd_apply_info 
 * \brief info to apply boundary condition
 */ 
typedef struct Bd_apply_info
{
    ivector dof, bd, idx;
    // dof = interiori node
    // bd = boundary node
    // idx = mapping from node_index to dof_idx/bd_idx, node_flag = 0, if is a interiori node, node_flag = -1, if else
} Bd_apply_info;

extern double f(double *p);
extern double u(double *p);

int setup_heat (dCSRmat *A_heat,
                dCSRmat *Mass,
                dvector *b_heat,
                Mesh *mesh,
                Mesh_aux *mesh_aux,
                FEM_param *pt,
                dvector *uh_heat,
                Bd_apply_info *bdinfo,
                double dt);

double get_l2_error_heat (ddenmat *node,
                          idenmat *elem,
                          dvector *uh,
                          int num_qp,
                          double t);

#endif
