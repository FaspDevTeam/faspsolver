/*! \file misc.h
 *  \brief miscellaneous for test.
 */

#ifndef _MISC_H_
#define _MISC_H_

#include <string.h>

/** 
 * \struct param_test
 * \brief param for test
 */ 
typedef struct param_test{
	
  char meshIn[128];
	char meshOut[128];
  char option[128];

	int refine_lvl;
	int nt;
	double T;
	int mesh_out;
	int num_qp_rhs;
	int num_qp_mat;
	
} param_test;


void param_init(param_test *pt);

int param_set(int argc, const char *argv [], param_test * pt);

#endif
