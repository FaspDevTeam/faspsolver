/*! \file amg_setup_aggregation_bsr.inl
 *  \brief Utilies for multigrid cycles in BSR format
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void form_tentative_p_bsr (ivector *vertices, dBSRmat *tentp, AMG_data_bsr *mgl,
 *                                       INT levelNum, INT num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighborhoods
 *
 * \param vertices           Pointer to the aggregation of vertices
 * \param tentp              Pointer to the prolongation operators
 * \param mgl                Pointer to AMG levele data
 * \param levelNum           Level number
 * \param num_aggregations   Number of aggregations
 *
 * \author Xiaozhe Hu
 * \date   11/25/2011
 */
static void form_tentative_p_bsr (ivector *vertices,
                                  dBSRmat *tentp,
                                  AMG_data_bsr *mgl,
                                  INT levelNum,
                                  INT num_aggregations)
{
    INT i, j;
    
    /* Form tentative prolongation */
    tentp->ROW = vertices->row;
    tentp->COL = num_aggregations;
    tentp->nb = mgl->A.nb;
    INT nb2 = tentp->nb * tentp->nb;
    
    tentp->IA  = (INT*)fasp_mem_calloc(tentp->ROW+1,sizeof(INT));
    
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
    tentp->val = (REAL*)fasp_mem_calloc(tentp->NNZ*nb2,sizeof(REAL));
    
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
