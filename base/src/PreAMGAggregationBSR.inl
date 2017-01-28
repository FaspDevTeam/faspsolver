/*! \file  PreAMGAggregationBSR.inl
 *
 *  \brief Utilities for aggregation methods for BSR matrices
 *
 *  \note This file contains Level-4 (Pre) functions, which are used in
 *        PreAMGSetupSABSR.c and PreAMGSetupUABSR.c
 * 
 *  \warning This file is also used in FASP4BLKOIL!!!
 *
 *  // TODO: Unused functions! --Chensong
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void form_tentative_p_bsr (ivector *vertices, dBSRmat *tentp,
 *                                       AMG_data_bsr *mgl, INT levelNum,
 *                                       INT num_agg)
 *
 * \brief Form tentative prolongation for BSR format matrix
 *
 * \param vertices    Pointer to the aggregation of vertices
 * \param tentp       Pointer to the prolongation operators
 * \param mgl         Pointer to AMG levels
 * \param levelNum    Level number
 * \param num_agg     Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   11/25/2011
 */
static void form_tentative_p_bsr (ivector *vertices,
                                  dBSRmat *tentp,
                                  AMG_data_bsr *mgl,
                                  INT levelNum,
                                  INT num_agg)
{
    INT i, j;
    
    /* Form tentative prolongation */
    tentp->ROW = vertices->row;
    tentp->COL = num_agg;
    tentp->nb  = mgl->A.nb;
    INT nb2    = tentp->nb * tentp->nb;
    
    tentp->IA  = (INT*)fasp_mem_calloc(tentp->ROW+1, sizeof(INT));
    
    // local variables
    INT *IA = tentp->IA;
    INT *JA;
    REAL *val;
    INT *vval = vertices->val;
    
    const INT row = tentp->ROW;
    
    // first run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            j ++;
        }
    }
    IA[row] = j;
    
    // allocate
    tentp->NNZ = j;
    tentp->JA = (INT*)fasp_mem_calloc(tentp->NNZ, sizeof(INT));
    tentp->val = (REAL*)fasp_mem_calloc(tentp->NNZ*nb2, sizeof(REAL));
    
    JA = tentp->JA;
    val = tentp->val;
    
    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            JA[j] = vval[i];
            fasp_smat_identity (&(val[j*nb2]), tentp->nb, nb2);
            j ++;
        }
    }
}

/**
 * \fn static void form_boolean_p_bsr (ivector *vertices, dBSRmat *tentp, AMG_data_bsr *mgl,
 *                                     INT levelNum, INT num_agg)
 *
 * \brief Form boolean prolongations in dBSRmat (only assume constant vector is in the null space)
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param mgl                Pointer to AMG levels
 * \param levelNum           Level number
 * \param num_agg            Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   05/27/2014
 */
static void form_boolean_p_bsr (ivector *vertices,
                                dBSRmat *tentp,
                                AMG_data_bsr *mgl,
                                INT levelNum,
                                INT num_agg)
{
    INT i, j;
    
    /* Form tentative prolongation */
    tentp->ROW = vertices->row;
    tentp->COL = num_agg;
    tentp->nb  = mgl->A.nb;
    INT nb2    = tentp->nb * tentp->nb;
    
    tentp->IA  = (INT*)fasp_mem_calloc(tentp->ROW+1, sizeof(INT));
    
    // local variables
    INT * IA = tentp->IA;
    INT *JA;
    REAL *val;
    INT *vval = vertices->val;
    
    const INT row = tentp->ROW;
    
    // first run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            j ++;
        }
    }
    IA[row] = j;
    
    // allocate
    tentp->NNZ = j;
    
    tentp->JA = (INT*)fasp_mem_calloc(tentp->NNZ, sizeof(INT));
    
    tentp->val = (REAL*)fasp_mem_calloc(tentp->NNZ*nb2, sizeof(REAL));
    
    
    JA = tentp->JA;
    val = tentp->val;
    
    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            JA[j] = vval[i];
            fasp_smat_identity (&(val[j*nb2]), tentp->nb, nb2);
            j ++;
        }
    }
}

/**
 * \fn static void form_tentative_p_bsr1 (ivector *vertices, dBSRmat *tentp,
 *                                        AMG_data_bsr *mgl, INT levelNum,
 *                                        INT num_agg, const INT dim, REAL **basis)
 *
 * \brief Form tentative prolongation for BSR format matrix (use general basis for null space)
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param mgl                Pointer to AMG levels
 * \param levelNum           Level number
 * \param num_agg            Number of aggregations
 * \param dim                Dimension of the near kernel space
 * \param basis              Pointer to the basis of the near kernel space
 *
 * \author Xiaozhe Hu
 * \date   05/27/2014
 */
static void form_tentative_p_bsr1 (ivector *vertices,
                                   dBSRmat *tentp,
                                   AMG_data_bsr *mgl,
                                   INT levelNum,
                                   INT num_agg,
                                   const INT dim,
                                   REAL **basis)
{
    INT i, j, k;
    
    INT p, q;
    
    const INT nnz_row = dim/mgl->A.nb; // nonzeros per row
    
    /* Form tentative prolongation */
    tentp->ROW = vertices->row;
    tentp->COL = num_agg*nnz_row;
    tentp->nb = mgl->A.nb;
    const INT nb = tentp->nb;
    const INT nb2 = nb * nb;
    
    tentp->IA  = (INT*)fasp_mem_calloc(tentp->ROW+1, sizeof(INT));
    
    // local variables
    INT *IA = tentp->IA;
    INT *JA;
    REAL *val;
    INT *vval = vertices->val;
    
    const INT row = tentp->ROW;
    
    // first run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            j = j + nnz_row;
        }
    }
    IA[row] = j;
    
    // allocate
    tentp->NNZ = j;
    tentp->JA = (INT*)fasp_mem_calloc(tentp->NNZ, sizeof(INT));
    tentp->val = (REAL*)fasp_mem_calloc(tentp->NNZ*nb2, sizeof(REAL));
    
    JA = tentp->JA;
    val = tentp->val;
    
    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            
            for (k=0; k<nnz_row; k++){
                
                JA[j] = vval[i]*nnz_row + k;
                
                for (p=0; p<nb; p++){
                    
                    for (q=0; q<nb; q++){
                        
                        val[j*nb2 + p*nb + q] = basis[k*nb+p][i*nb+q];
                        
                    }
                    
                }
                
                j++;
                
            }
        }
    }
}

/**
 * \fn static void smooth_agg_bsr (dBSRmat *A, dBSRmat *tentp, dBSRmat *P,
 *                                 AMG_param *param, INT levelNum, dCSRmat *N)
 *
 * \brief Smooth the tentative prolongation
 *
 * \param A         Pointer to the coefficient matrices (dBSRmat)
 * \param tentp     Pointer to the tentative prolongation operators (dBSRmat)
 * \param P         Pointer to the prolongation operators (dBSRmat)
 * \param param     Pointer to AMG parameters
 * \param levelNum  Current level number
 * \param N         Pointer to strongly coupled neighbors
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
static void smooth_agg_bsr (dBSRmat *A,
                            dBSRmat *tentp,
                            dBSRmat *P,
                            AMG_param *param,
                            INT levelNum,
                            dCSRmat *N)
{
    const SHORT filter = param->smooth_filter;
    const INT   row = A->ROW, col= A->COL, nnz = A->NNZ;
    const INT   nb = A->nb;
    const INT nb2 = nb*nb;
    const REAL  smooth_factor = param->tentative_smooth;
    
    // local variables
    dBSRmat S;
    dvector diaginv;  // diagonal block inv
    
    INT i,j;
    
    REAL *Id   = (REAL *)fasp_mem_calloc(nb2, sizeof(REAL));
    REAL *temp = (REAL *)fasp_mem_calloc(nb2, sizeof(REAL));
    
    fasp_smat_identity(Id, nb, nb2);
    
    /* Step 1. Form smoother */
    
    /* Without filter: Using A for damped Jacobian smoother */
    if ( filter != ON ) {
        
        // copy structure from A
        S = fasp_dbsr_create(row, col, nnz, nb, 0);
        
        for ( i=0; i<=row; ++i ) S.IA[i] = A->IA[i];
        for ( i=0; i<nnz; ++i ) S.JA[i] = A->JA[i];
        
        diaginv = fasp_dbsr_getdiaginv(A);
        
        // for S
        for (i=0; i<row; ++i) {
            
            for (j=S.IA[i]; j<S.IA[i+1]; ++j) {
                
                if (S.JA[j] == i) {
                    
                    fasp_blas_smat_mul(diaginv.val+(i*nb2), A->val+(j*nb2), temp, nb);
                    fasp_blas_smat_add(Id, temp, nb, 1.0, (-1.0)*smooth_factor, S.val+(j*nb2));
                    
                }
                else {
                    
                    fasp_blas_smat_mul(diaginv.val+(i*nb2), A->val+(j*nb2), S.val+(j*nb2), nb);
                    fasp_blas_smat_axm(S.val+(j*nb2), nb, (-1.0)*smooth_factor);
                    
                }
                
            }
            
        }
    }
    
    fasp_dvec_free(&diaginv);
    fasp_mem_free(Id);
    fasp_mem_free(temp);
    
    /* Step 2. Smooth the tentative prolongation P = S*tenp */
    fasp_blas_dbsr_mxm(&S, tentp, P); // Note: think twice about this.
    
    P->NNZ = P->IA[P->ROW];
    
    fasp_dbsr_free(&S);
}

/**
 * \fn static void smooth_agg_bsr1 (dBSRmat *A, dBSRmat *tentp, dBSRmat *P,
 *                                  AMG_param *param, INT levelNum, dCSRmat *N)
 *
 * \brief Smooth the tentative prolongation (without filter)
 *
 * \param A         Pointer to the coefficient matrices (dBSRmat)
 * \param tentp     Pointer to the tentative prolongation operators (dBSRmat)
 * \param P         Pointer to the prolongation operators (dBSRmat)
 * \param param     Pointer to AMG parameters
 * \param levelNum  Current level number
 * \param N         Pointer to strongly coupled neighbors
 *
 * \author Xiaozhe Hu
 * \date   05/26/2014
 */
static void smooth_agg_bsr1 (dBSRmat *A,
                             dBSRmat *tentp,
                             dBSRmat *P,
                             AMG_param *param,
                             INT levelNum,
                             dCSRmat *N)
{
    const INT   row = A->ROW;
    const INT   nb = A->nb;
    const INT   nb2 = nb*nb;
    const REAL  smooth_factor = param->tentative_smooth;
    
    // local variables
    dBSRmat S;
    
    INT i,j;
    
    REAL *Id   = (REAL *)fasp_mem_calloc(nb2, sizeof(REAL));
    REAL *temp = (REAL *)fasp_mem_calloc(nb2, sizeof(REAL));
    
    fasp_smat_identity(Id, nb, nb2);
    
    /* Step 1. D^{-1}A */
    S = fasp_dbsr_diaginv(A);
    
    /* Step 2. -wD^{-1}A */
    fasp_blas_dbsr_axm (&S, (-1.0)*smooth_factor);
    
    /* Step 3. I - wD^{-1}A */
    for (i=0; i<row; ++i) {
        for (j=S.IA[i]; j<S.IA[i+1]; ++j) {
            if (S.JA[j] == i) {
                
                fasp_blas_smat_add(Id, S.val+(j*nb2), nb, 1.0, 1.0, temp);
                fasp_array_cp(nb2, temp, S.val+(j*nb2));
                
            }
        }
    }
    
    fasp_mem_free(Id);
    fasp_mem_free(temp);
    
    /* Step 2. Smooth the tentative prolongation P = S*tenp */
    fasp_blas_dbsr_mxm(&S, tentp, P); // Note: think twice about this.
    
    P->NNZ = P->IA[P->ROW];
    
    fasp_dbsr_free(&S);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
