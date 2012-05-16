/*! \file assemble.h
 *  \brief auxiliary func for assembling
 */

#ifndef _MESH_H_
#define _MESH_H_

#include "fasp.h"
#include "fasp_functs.h"

#define DIRICHLET 1
#define INTERIORI 0

void ivec_output ( ivector *t, ivector *s);
void dvec_output ( dvector *t, dvector *s);
void icsr_output ( iCSRmat *t, iCSRmat *s);
void dcsr_output ( dCSRmat *t, dCSRmat *s);
void dden_output ( ddenmat *t, ddenmat *s);

void extractNondirichletMatrix (dCSRmat *A, 
                                dvector *b, 
                                dCSRmat *A11, 
                                dvector *b1, 
                                ivector *isInNode, 
                                ivector *dirichlet, 
                                ivector *nondirichlet, 
                                ivector *index, 
                                dvector *uh);

#endif