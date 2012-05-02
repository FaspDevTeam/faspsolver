/*! \file smoother_bsr.c
 *  \brief Smoothers for sparse matrix in BSR format
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

////////////////////////////////////////  JACOBI  //////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_jacobi (dBSRmat *A, dvector *b, dvector *u)
 *
 * \brief Jacobi relaxation
 *
 * \param A   Pointer to coefficient matrix
 * \param b   Pointer to right hand side vector
 * \param u   Initial guess (in) and new approximation after one iteration 
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25
 */
void fasp_smoother_dbsr_jacobi (dBSRmat *A, 
                                dvector *b, 
                                dvector *u)
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb2;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    const REAL   *val = A->val;
    
    // local variables
    INT i,k;
    REAL *diaginv = NULL;
    
    // allocate memory   
    diaginv = (REAL *)fasp_mem_calloc(size, sizeof(REAL));
    
    // get all the diagonal sub-blocks   
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i+1]; ++k) {
            if (JA[k] == i)
                memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(REAL));
        }
    }
    
    // compute the inverses of all the diagonal sub-blocks   
    if (nb > 1) {
        for (i = 0; i < ROW; ++i) {
            fasp_blas_smat_inv(diaginv+i*nb2, nb);
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {  
            // zero-diagonal should be tested previously
            diaginv[i] = 1.0 / diaginv[i];
        }
    }
    
    fasp_smoother_dbsr_jacobi1(A, b, u, diaginv);
    fasp_mem_free(diaginv);    
}


/**
 * \fn void fasp_smoother_dbsr_jacobi1 (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv)
 *
 * \brief Jacobi relaxation
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_jacobi1 (dBSRmat *A, 
                                 dvector *b, 
                                 dvector *u, 
                                 REAL *diaginv)    
{    
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // auxiliary array
    REAL *b_tmp = NULL;
    
    // local variables
    INT i,j,k;
    INT pb;
    
    // b_tmp = b_val
    b_tmp = (REAL *)fasp_mem_calloc(size, sizeof(REAL));
    memcpy(b_tmp, b_val, size*sizeof(REAL));
    
    // It's not necessary to assign the smoothing order since the result doesn't depend on it
    if (nb == 1) {
        for (i = 0; i < ROW; ++i) {
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    b_tmp[i] -= val[k]*u_val[j];
            }
        }
    
        for (i = 0; i < ROW; ++i) {
            u_val[i] = b_tmp[i]*diaginv[i]; 
        }      

        fasp_mem_free(b_tmp);   
    }
    else if (nb > 1) {
        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp+pb, nb);
            }
        }
    
        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp+pb, u_val+pb, nb);
        }     

        fasp_mem_free(b_tmp);   
    }
    else {
        printf("### ERROR: nb is illegal!\n");
        exit(RUN_FAIL);
    }
    
}

////////////////////////////////////  Gauss-Seidel  ////////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_gs (dBSRmat *A, dvector *b, dvector *u, INT order, INT *mark)
 *
 * \brief Gauss-Seidel relaxation
 *
 * \param A      Pointer to coefficient matrix
 * \param b      Pointer to right hand side vector
 * \param u      Initial guess (in) and new approximation after one iteration 
 * \param order  Flag to indicate the order for smoothing
 *               If mark = NULL       
 *                    ASCEND       12: in ascending order
 *                    DESCEND      21: in descending order
 *               If mark != NULL:  in the user-defined order
 * \param mark   Pointer to NULL or to the user-defined ordering
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_gs (dBSRmat *A, 
                            dvector *b, 
                            dvector *u, 
                            INT order, 
                            INT *mark )
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb2;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    const REAL   *val = A->val;
    
    // local variables
    INT i,k;
    
    // allocate memory
    REAL *diaginv = (REAL *)fasp_mem_calloc(size, sizeof(REAL));
    
    // get all the diagonal sub-blocks   
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i+1]; ++k) {
            if (JA[k] == i)
                memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(REAL));
        }
    }
    
    // compute the inverses of all the diagonal sub-blocks   
    if (nb > 1) {
        for (i = 0; i < ROW; ++i) {
            fasp_blas_smat_inv(diaginv+i*nb2, nb);
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {  
            // zero-diagonal should be tested previously
            diaginv[i] = 1.0 / diaginv[i];
        }
    }
    
    fasp_smoother_dbsr_gs1(A, b, u, order, mark, diaginv);
    fasp_mem_free(diaginv);
}

/**
 * \fn void fasp_smoother_dbsr_gs1 (dBSRmat *A, dvector *b, dvector *u, INT order, INT *mark, REAL *diaginv)
 *
 * \brief Gauss-Seidel relaxation
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param order    Flag to indicate the order for smoothing
 *                 If mark = NULL       
 *                    ASCEND       12: in ascending order
 *                    DESCEND      21: in descending order
 *                 If mark != NULL:  in the user-defined order
 * \param mark     Pointer to NULL or to the user-defined ordering
 * \param diaginv  Inverses for all the diagonal blocks of A
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_gs1 (dBSRmat *A, 
                             dvector *b, 
                             dvector *u, 
                             INT order, 
                             INT *mark, 
                             REAL *diaginv )    
{    
    if (!mark) {
        if (order == ASCEND) // smooth ascendingly
            {
                fasp_smoother_dbsr_gs_ascend(A, b, u, diaginv);
            }
        else if (order == DESCEND) // smooth descendingly
            {
                fasp_smoother_dbsr_gs_descend(A, b, u, diaginv);
            }
    }
    // smooth according to the order 'mark' defined by user
    else {
        fasp_smoother_dbsr_gs_order1(A, b, u, diaginv, mark);
    }
}

/**
 * \fn void fasp_smoother_dbsr_gs_ascend (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv)
 *
 * \brief Gauss-Seidel relaxation in the ascending order
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A 
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_gs_ascend (dBSRmat *A, 
                                   dvector *b, 
                                   dvector *u, 
                                   REAL *diaginv )
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // local variables
    INT   i,j,k;
    INT   pb;
    REAL  rhs = 0.0;
    
    if (nb == 1) {
        for (i = 0; i < ROW; ++i) {
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {  
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = rhs*diaginv[i];
        }  
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)fasp_mem_calloc(nb, sizeof(REAL));

        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
        }

        fasp_mem_free(b_tmp);
    }
    else {
        printf("### ERROR: nb is illegal!\n");
        exit(RUN_FAIL);
    }

}

/**
 * \fn void fasp_smoother_dbsr_gs_descend (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv)
 *
 * \brief Gauss-Seidel relaxation in the descending order
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A 
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_gs_descend (dBSRmat *A, 
                                    dvector *b, 
                                    dvector *u, 
                                    REAL *diaginv )
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // local variables
    INT i,j,k;
    INT pb;
    REAL rhs = 0.0;
    
    if (nb == 1) {
        for (i = ROW-1; i >= 0; i--) {
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {  
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = rhs*diaginv[i];
        }  
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)fasp_mem_calloc(nb, sizeof(REAL));

        for (i = ROW-1; i >= 0; i--) {
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
        }

        fasp_mem_free(b_tmp);
    }
    else {
        printf("### ERROR: nb is illegal!\n");
        exit(RUN_FAIL);
    }
    
}

/**
 * \fn void fasp_smoother_dbsr_gs_order1 (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv, INT *mark)
 *
 * \brief Gauss-Seidel relaxation in the user-defined order
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A 
 * \param mark     Pointer to the user-defined ordering  
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_gs_order1 (dBSRmat *A, 
                                   dvector *b, 
                                   dvector *u, 
                                   REAL *diaginv, 
                                   INT *mark )
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // local variables
    INT i,j,k;
    INT I,pb;
    REAL rhs = 0.0;
    
    if (nb == 1) {
        for (I = 0; I < ROW; ++I) {
            i = mark[I];
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = rhs*diaginv[i];
        }  
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)fasp_mem_calloc(nb, sizeof(REAL));

        for (I = 0; I < ROW; ++I) {
            i = mark[I];
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            fasp_blas_smat_mxv(diaginv+nb2*i, b_tmp, u_val+pb, nb);
        }

        fasp_mem_free(b_tmp);
    }
    else {
        fasp_chkerr(ERROR_NUM_BLOCKS,"fasp_smoother_dbsr_gs_order1");
    }
    
}

/**
 * \fn void fasp_smoother_dbsr_gs_order2 (dBSRmat *A, dvector *b, dvector *u, INT *mark, REAL *work)
 *
 * \brief Gauss-Seidel relaxation in the user-defined order
 *
 * \param A      Pointer to coefficient matrix
 * \param b      Pointer to right hand side vector
 * \param u      Initial guess (in) and new approximation after one iteration 
 * \param mark   Pointer to the user-defined ordering  
 * \param work   Work temp array
 *
 * \author Zhiyang Zhou
 * \date   2010/11/08 
 *
 * \note The only difference between the functions 'fasp_smoother_dbsr_gs_order2' 
 *       and 'fasp_smoother_dbsr_gs_order1' lies in that we don't have to multiply 
 *       by the inverses of the diagonal blocks in each ROW since matrix A has  
 *       been such scaled that all the diagonal blocks become identity matrices.
 */
void fasp_smoother_dbsr_gs_order2 (dBSRmat *A, 
                                   dvector *b, 
                                   dvector *u, 
                                   INT *mark, 
                                   REAL *work)
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // auxiliary array
    REAL *b_tmp = work;
    
    // local variables
    INT i,j,k,I,pb;
    REAL rhs = 0.0;
    
    if (nb == 1) {
        for (I = 0; I < ROW; ++I) {
            i = mark[I];
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {  
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = rhs;
        }  
    }
    else if (nb > 1) {
        for (I = 0; I < ROW; ++I) {
            i = mark[I];
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            memcpy(u_val+pb, b_tmp, nb*sizeof(REAL));
        }
    }
    else {
        fasp_chkerr(ERROR_NUM_BLOCKS,"fasp_smoother_dbsr_gs_order2");
    }
}

////////////////////////////////////////  SOR  //////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_sor (dBSRmat *A, dvector *b, dvector *u, INT order, INT *mark, REAL weight)
 *
 * \brief SOR relaxation
 *
 * \param A      Pointer to coefficient matrix
 * \param b      Pointer to right hand side vector
 * \param u      Initial guess (in) and new approximation after one iteration 
 * \param order  Flag to indicate the order for smoothing
 *               If mark = NULL       
 *                    ASCEND       12: in ascending order
 *                    DESCEND      21: in descending order
 *               If mark != NULL:  in the user-defined order
 * \param mark   Pointer to NULL or to the user-defined ordering
 * \param weight Over-relaxation weight 
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_sor (dBSRmat *A, 
                             dvector *b, 
                             dvector *u, 
                             INT order, 
                             INT *mark, 
                             REAL weight)
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    size = ROW*nb2;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    
    // local variables
    REAL *diaginv = NULL;
    INT i,k;
    
    // allocate memory
    diaginv = (REAL *)fasp_mem_calloc(size, sizeof(REAL));
    
    // get all the diagonal sub-blocks   
    for (i = 0; i < ROW; ++i) {
        for (k = IA[i]; k < IA[i+1]; ++k) {
            if (JA[k] == i)
                memcpy(diaginv+i*nb2, val+k*nb2, nb2*sizeof(REAL));
        }
    }
    
    // compute the inverses of all the diagonal sub-blocks   
    if (nb > 1) {
        for (i = 0; i < ROW; ++i) {
            fasp_blas_smat_inv(diaginv+i*nb2, nb);
        }
    }
    else {
        for (i = 0; i < ROW; ++i) {  
            // zero-diagonal should be tested previously
            diaginv[i] = 1.0 / diaginv[i];
        }
    }
    
    fasp_smoother_dbsr_sor1(A, b, u, order, mark, diaginv, weight);
    fasp_mem_free(diaginv);
}

/**
 * \fn void fasp_smoother_dbsr_sor1 (dBSRmat *A, dvector *b, dvector *u, INT order, 
 *                                   INT *mark, REAL *diaginv, REAL weight)
 *
 * \brief SOR relaxation
 *
 * \param A       Pointer to coefficient matrix
 * \param b       Pointer to right hand side vector
 * \param u       Initial guess (in) and new approximation after one iteration 
 * \param order   Flag to indicate the order for smoothing
 *                If mark = NULL       
 *                    ASCEND       12: in ascending order
 *                    DESCEND      21: in descending order
 *                If mark != NULL:  in the user-defined order
 * \param mark    Pointer to NULL or to the user-defined ordering
 * \param diaginv Inverses for all the diagonal blocks of A
 * \param weight  Over-relaxation weight 
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_sor1 (dBSRmat *A, 
                              dvector *b, 
                              dvector *u, 
                              INT order, 
                              INT *mark, 
                              REAL *diaginv, 
                              REAL weight)    
{    
    if (!mark) {
        if (order == ASCEND)       // smooth ascendingly
            {
                fasp_smoother_dbsr_sor_ascend(A, b, u, diaginv, weight);
            }
        else if (order == DESCEND) // smooth descendingly
            {
                fasp_smoother_dbsr_sor_descend(A, b, u, diaginv, weight);
            }
    }
    // smooth according to the order 'mark' defined by user
    else {
        fasp_smoother_dbsr_sor_order(A, b, u, diaginv, mark, weight);
    }
}

/**
 * \fn void fasp_smoother_dbsr_sor_ascend (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv, REAL weight)
 *
 * \brief SOR relaxation in the ascending order
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A 
 * \param weight   Over-relaxation weight   
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_sor_ascend (dBSRmat *A, 
                                    dvector *b, 
                                    dvector *u, 
                                    REAL *diaginv, 
                                    REAL weight )
{
    // members of A 
    INT     ROW = A->ROW;
    INT     nb  = A->nb;
    INT    *IA  = A->IA;
    INT    *JA  = A->JA;   
    REAL *val = A->val;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // local variables
    INT nb2 = nb*nb;
    INT i,j,k;
    INT pb;
    REAL rhs = 0.0;
    REAL one_minus_weight = 1.0 - weight;
    
    if (nb == 1) {
        for (i = 0; i < ROW; ++i) {
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {  
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = one_minus_weight*u_val[i] + weight*(rhs*diaginv[i]);         
        }  
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)fasp_mem_calloc(nb, sizeof(REAL));

        for (i = 0; i < ROW; ++i) {
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            fasp_blas_smat_aAxpby(weight, diaginv+nb2*i, b_tmp, one_minus_weight, u_val+pb, nb);         
        }

        fasp_mem_free(b_tmp);
    }
    else {
        fasp_chkerr(ERROR_NUM_BLOCKS,"fasp_smoother_dbsr_sor_ascend");
    }
    
}

/**
 * \fn void fasp_smoother_dbsr_sor_descend (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv, REAL weight)
 *
 * \brief SOR relaxation in the descending order
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A 
 * \param weight   Over-relaxation weight   
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_sor_descend (dBSRmat *A, 
                                     dvector *b, 
                                     dvector *u, 
                                     REAL *diaginv, 
                                     REAL weight)
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    const REAL    one_minus_weight = 1.0 - weight;

    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // local variables
    INT i,j,k;
    INT pb;
    REAL rhs = 0.0;
    
    if (nb == 1) {
        for (i = ROW-1; i >= 0; i--) {
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {  
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = one_minus_weight*u_val[i] + weight*(rhs*diaginv[i]);         
        }  
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)fasp_mem_calloc(nb, sizeof(REAL));

        for (i = ROW-1; i >= 0; i--) {
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            fasp_blas_smat_aAxpby(weight, diaginv+nb2*i, b_tmp, one_minus_weight, u_val+pb, nb);         
        }

        fasp_mem_free(b_tmp);
    }
    else {
        fasp_chkerr(ERROR_NUM_BLOCKS,"fasp_smoother_dbsr_sor_descend");
    }
    
}

/**
 * \fn void fasp_smoother_dbsr_sor_order (dBSRmat *A, dvector *b, dvector *u, REAL *diaginv, 
 *                                        INT *mark, REAL weight)
 *
 * \brief SOR relaxation in the user-defined order
 *
 * \param A        Pointer to coefficient matrix
 * \param b        Pointer to right hand side vector
 * \param u        Initial guess (in) and new approximation after one iteration 
 * \param diaginv  Inverses for all the diagonal blocks of A 
 * \param mark     Pointer to the user-defined ordering   
 * \param weight   Over-relaxation weight   
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_sor_order (dBSRmat *A, 
                                   dvector *b, 
                                   dvector *u, 
                                   REAL *diaginv, 
                                   INT *mark, 
                                   REAL weight )
{
    // members of A 
    const INT     ROW = A->ROW;
    const INT     nb  = A->nb;
    const INT     nb2 = nb*nb;
    const INT    *IA  = A->IA;
    const INT    *JA  = A->JA;   
    REAL         *val = A->val;
    const REAL    one_minus_weight = 1.0 - weight;
    
    // values of dvector b and u
    REAL *b_val = b->val;
    REAL *u_val = u->val;
    
    // local variables
    INT i,j,k;
    INT I,pb;
    REAL rhs = 0.0;
    
    if (nb == 1) {
        for (I = 0; I < ROW; ++I) {
            i = mark[I];
            rhs = b_val[i];
            for (k = IA[i]; k < IA[i+1]; ++k) {  
                j = JA[k];
                if (j != i)
                    rhs -= val[k]*u_val[j];
            }
            u_val[i] = one_minus_weight*u_val[i] + weight*(rhs*diaginv[i]);         
        }  
    }
    else if (nb > 1) {
        REAL *b_tmp = (REAL *)fasp_mem_calloc(nb, sizeof(REAL));

        for (I = 0; I < ROW; ++I) {
            i = mark[I];
            pb = i*nb;
            memcpy(b_tmp, b_val+pb, nb*sizeof(REAL));
            for (k = IA[i]; k < IA[i+1]; ++k) { 
                j = JA[k];
                if (j != i)
                    fasp_blas_smat_ymAx(val+k*nb2, u_val+j*nb, b_tmp, nb);
            }
            fasp_blas_smat_aAxpby(weight, diaginv+nb2*i, b_tmp, one_minus_weight, u_val+pb, nb);         
        }
        
        fasp_mem_free(b_tmp);
    }
    else {
        fasp_chkerr(ERROR_NUM_BLOCKS,"fasp_smoother_dbsr_sor_order");
    }
    
}

////////////////////////////////////////  ILU  //////////////////////////////////////////

/**
 * \fn void fasp_smoother_dbsr_ilu (dBSRmat *A, dvector *b, dvector *x, void *data)
 *
 * \brief ILU method as the smoother in solving Au=b with multigrid method
 *
 * \param A     Pointer to stiffness matrix
 * \param b     Pointer to right hand side
 * \param x     Pointer to current solution
 * \param data  Pointer to user defined data
 *
 * \author Zhiyang Zhou
 * \date   2010/10/25 
 */
void fasp_smoother_dbsr_ilu (dBSRmat *A, 
                             dvector *b, 
                             dvector *x, 
                             void *data)
{
    ILU_data   *iludata=(ILU_data *)data;
    const INT   nb=iludata->nb,m=A->ROW*nb, memneed=6*m;
    
    REAL *xval = x->val, *bval = b->val;    
    REAL *zz = iludata->work + 3*m;
    REAL *zr = zz + m;
    REAL *z  = zr + m;
    
    if (iludata->nwork<memneed) goto MEMERR; 

    /** form residual zr = b - A x */
    fasp_array_cp(m,bval,zr); fasp_blas_dbsr_aAxpy(-1.0,A,xval,zr);
    
    /** solve LU z=zr */
    fasp_precond_dbsr_ilu(zr,z,iludata);
    
    /** x=x+z */
    fasp_blas_array_axpy(m,1,z,xval);    
    
    return;
    
 MEMERR:
    printf("### ERROR: Need %d memory, only %d available!!!\n", memneed, iludata->nwork);
    exit(ERROR_ALLOC_MEM);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
