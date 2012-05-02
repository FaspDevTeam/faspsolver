/*! \file smat.c
 *  \brief Simple operations for small full matrices in row-major format
 *
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_iden_free (idenmat *A)
 *
 * \brief Free idenmat sparse matrix data memeory space
 *
 * \param A   Pointer to the idenmat matrix
 *
 * \author Chensong Zhang
 * \date   2010/04/03 
 */
void fasp_iden_free (idenmat *A)
{    
    unsigned INT i;
    
    if (A==NULL) return;
    
    for (i=0;i<A->row;++i) fasp_mem_free(A->val[i]);
    fasp_mem_free(A->val); A->val = NULL; A->row = 0;
}

/**
 * \fn void fasp_smat_identity_nc2 (REAL *a)
 *
 * \brief Set a 2*2 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 2*2 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2011/11/18
 */
void fasp_smat_identity_nc2 (REAL *a)
{
    memset(a, 0X0, 4*sizeof(REAL));
    
    a[0] = 1.0; a[3] = 1.0; 
}    

/**
 * \fn void fasp_smat_identity_nc3 (REAL *a)
 *
 * \brief Set a 3*3 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 3*3 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity_nc3 (REAL *a)
{
    memset(a, 0X0, 9*sizeof(REAL));
    
    a[0] = 1.0; a[4] = 1.0; a[8] = 1.0;
}    

/**
 * \fn void fasp_smat_identity_nc5 (REAL *a)
 *
 * \brief Set a 5*5 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 5*5 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity_nc5 (REAL *a)
{
    memset(a, 0X0, 25*sizeof(REAL));
    
    a[0] = 1.0; a[6] = 1.0; a[12] = 1.0; a[18] = 1.0; a[24] = 1.0;
}    

/**
 * \fn void fasp_smat_identity_nc7 (REAL *a)
 *
 * \brief Set a 7*7 full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a 7*7 full matrix 
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity_nc7 (REAL *a)
{
    memset(a, 0X0, 49*sizeof(REAL));
    
    a[0] = 1.0; a[8] = 1.0; a[16] = 1.0; a[24] = 1.0;
    a[32] = 1.0; a[40] = 1.0; a[48] = 1.0;
}    

/**
 * \fn void fasp_smat_identity (REAL *a, INT n, INT n2)
 *
 * \brief Set a n*n full matrix to be a identity   
 *
 * \param a      Pointer to the REAL vector which stands for a n*n full matrix 
 * \param n      Size of full matrix
 * \param n2     Length of the REAL vector which stores the n*n full matrix
 *
 * \author Xiaozhe Hu
 * \date   2010/12/25
 */
void fasp_smat_identity (REAL *a, 
                         INT n, 
                         INT n2)
{
    memset(a, 0X0, n2*sizeof(REAL));
    
    switch (n) {

    case 2: {
        a[0] = 1.0;
        a[3] = 1.0;
    }
        break;
    
    case 3: {
        a[0] = 1.0;
        a[4] = 1.0;
        a[8] = 1.0;
    }
        break;
    
    case 4: {
        a[0] = 1.0;
        a[5] = 1.0;
        a[10] = 1.0;
        a[15] = 1.0;
    }
        break;
    
    case 5: {
        a[0] = 1.0;
        a[6] = 1.0;
        a[12] = 1.0;
        a[18] = 1.0;
        a[24] = 1.0;
    }
        break;
    
    case 6: {
        a[0] = 1.0;
        a[7] = 1.0;
        a[14] = 1.0;
        a[21] = 1.0;
        a[28] = 1.0;
        a[35] = 1.0;
    }
        break;
    
    case 7: {
        a[0] = 1.0;
        a[8] = 1.0;
        a[16] = 1.0;
        a[24] = 1.0;
        a[32] = 1.0;
        a[40] = 1.0;
        a[48] = 1.0;
    }
        break;
    
    default: {
        INT l;
        for (l = 0; l < n; l ++) a[l*n+l] = 1.0;
    }
        break;
    }
    
}    

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
