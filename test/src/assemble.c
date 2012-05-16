/*! \file assemble.c
 *  \brief Auxiliary functions for assembling
 */

#include "assemble.h"

/** 
 * \fn void ivec_output ( ivector *t, ivector *s)
 *
 * \brief output ivector *t from *s.
 *
 * \param *t           ivector pointer for output 
 * \param *s           source pointer
 *
 * \author Feiteng Huang
 * \date   03/30/2012
 */
void ivec_output ( ivector *t, ivector *s)
{
    t->row = s->row;
    t->val = s->val;
}

/** 
 * \fn void dvec_output ( dvector *t, dvector *s)
 *
 * \brief output dvector *t from *s.
 *
 * \param *t           dvector pointer for output 
 * \param *s           source pointer
 *
 * \author Feiteng Huang
 * \date   03/30/2012
 */
void dvec_output ( dvector *t, dvector *s)
{
    t->row = s->row;
    t->val = s->val;
}

/** 
 * \fn void icsr_output ( icsrmat *t, icsrmat *s)
 *
 * \brief output iCSRmat *t from *s.
 *
 * \param *t           iCSRmat pointer for output 
 * \param *s           source pointer
 *
 * \author Feiteng Huang
 * \date   03/30/2012
 */
void icsr_output ( iCSRmat *t, iCSRmat *s)
{
    t->row = s->row;
    t->col = s->col;
    t->nnz = s->nnz;
    t->IA = s->IA;
    t->JA = s->JA;
    t->val = s->val;
}

/** 
 * \fn void dcsr_output ( dcsrmat *t, dcsrmat *s)
 *
 * \brief output dCSRmat *t from *s.
 *
 * \param *t           dCSRmat pointer for output 
 * \param *s           source pointer
 *
 * \author Feiteng Huang
 * \date   03/30/2012
 */
void dcsr_output ( dCSRmat *t, dCSRmat *s)
{
    t->row = s->row;
    t->col = s->col;
    t->nnz = s->nnz;
    t->IA = s->IA;
    t->JA = s->JA;
    t->val = s->val;
}

/** 
 * \fn void dden_output ( ddenmat *t, ddenmat *s)
 *
 * \brief output ddenmat *t from *s.
 *
 * \param *t           ddenmat pointer for output 
 * \param *s           source pointer
 *
 * \author Feiteng Huang
 * \date   03/30/2012
 */
void dden_output ( ddenmat *t, ddenmat *s)
{
    t->row = s->row;
    t->col = s->col;
    t->val = s->val;
    t->val[0] = s->val[0];
}

/**
 * \fn void extractNondirichletMatrix (dCSRmat *A, 
 *                                     dvector *b, 
 *                                     dCSRmat *A11, 
 *                                     dvector *b1, 
 *                                     ivector *isInNode, 
 *                                     ivector *dirichlet, 
 *                                     ivector *nondirichlet, 
 *                                     ivector *index, 
 *                                     dvector *uh)
 *
 * \brief Extract stiffness matrix by removing the corresponding dirichlet boundary condition  
 *
 * \param *A            pointer to the stiffness matrix with dirichelt boundary condition(without removed)
 * \param *b            pointer to the right hand side with dirichelt boundary condition(without removed)
 * \param *A11          pointer to the stiffness matrix without dirichelt boundary condition(removed)
 * \param *b1           pointer to the right hand side without dirichelt boundary condition(removed)
 * \param *isInNode     if the node is interior node, it will be 0
 *                      if the node is on the boundary, it will be -1
 *                      if it is dirichlet boundary it will be 1
 * \param *dirichlet    pointer to the indicator of the dirichlet boundary
 * \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
 * \param *index        pointer to the transpose of *dirichlet and *nondirichlet
 * \param *uh           pointer to the dirichlet boundary value
 *
 * \author Xuehai Huang, Kai Yang, and Feiteng Huang
 * \date   03/29/2009
 * 
 * Modified by Feiteng Huang on 04/12/2012
 */
void extractNondirichletMatrix (dCSRmat *A, 
                                dvector *b, 
                                dCSRmat *A11, 
                                dvector *b1, 
                                ivector *isInNode, 
                                ivector *dirichlet, 
                                ivector *nondirichlet, 
                                ivector *index, 
                                dvector *uh)
{
    // achiveve A11 due to dirichlet boundary condition
    int i,j,k,l,i1,j1;
    int count;
    
    A11->col=A11->row=A->row-dirichlet->row;
    A11->IA=(int*)fasp_mem_calloc(A11->row+1, sizeof(int));
    
    A11->JA=NULL;
    A11->val=NULL;
    b1->row=A11->row;
    b1->val=(double*)fasp_mem_calloc(b1->row, sizeof(double));
    
    // form A11->IA
    for (i=0;i<A11->row;++i) {
        l=nondirichlet->val[i];
        for (k=A->IA[l];k<A->IA[l+1];++k) {
            j=A->JA[k];
            if (isInNode->val[j]==INTERIORI) A11->IA[i+1]++;
        }
    }
    
    for (i=0;i<A11->row;++i) A11->IA[i+1]+=A11->IA[i];
    
    // form A11->JA
    A11->JA=(int*)fasp_mem_calloc(A11->IA[A11->row]+1, sizeof(int));
    
    count=0;
    for (i=0;i<A11->row;++i) {
        l=nondirichlet->val[i];
        for (k=A->IA[l];k<A->IA[l+1];++k) {
            j=A->JA[k];
            if (isInNode->val[j]==INTERIORI) {
                A11->JA[count]=index->val[j];
                ++count;
            }
        }
    }
    
    // form A11->val
    A11->val=(double*)fasp_mem_calloc(A11->IA[A11->row]+1, sizeof(double));
    
    for (i1=0;i1<A11->row;++i1) {
        i=nondirichlet->val[i1];
        for (k=A11->IA[i1];k<A11->IA[i1+1];++k) {
            j1=A11->JA[k];
            j=nondirichlet->val[j1];
            
            for (l=A->IA[i];l<A->IA[i+1];l++) {
                if (A->JA[l]==j) {
                    A11->val[k]=A->val[l];
                    break;
                }
            }
        }
    }
    A11->nnz=A11->IA[A11->row]-A11->IA[0];
    
    // yangkai: changed the loop order of previous code and reduced time.     
    for (i=0;i<A->row;++i) {
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            if (isInNode->val[A->JA[k]]==DIRICHLET) b->val[i]-=A->val[k]*uh->val[A->JA[k]]; // dirichlet point
        }
    }
    
    for (i1=0;i1<b1->row;++i1) {
        i=nondirichlet->val[i1];
        b1->val[i1]=b->val[i];
    }
    
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/