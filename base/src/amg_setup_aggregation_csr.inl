/*! \file amg_setup_aggregation_csr.inl
 *  \brief Utilies for multigrid cycles for CSR matrices
 */

#include <omp.h>
/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void aggregation (dCSRmat *A, ivector *vertices, AMG_param *param, 
 *                              INT levelNum, dCSRmat *Neigh, INT *num_aggregations)
 *
 * \brief Form aggregation based on strong coupled neighborhoods 
 *
 * \param A                 Pointer to the coefficient matrices
 * \param vertices          Pointer to the aggregation of vertics
 * \param param             Pointer to AMG parameters
 * \param levelNum          Level number
 * \param Neigh             Pointer to strongly coupled neighborhoods
 * \param num_aggregations  Pointer to number of aggregations 
 * 
 * \author Xiaozhe Hu
 * \date   09/29/2009
 *
 * \author Chunsheng Feng, Zheng Li
 * \date   09/03/2012
 */
static void aggregation (dCSRmat *A,
                         ivector *vertices, 
                         AMG_param *param, 
                         INT levelNum, 
                         dCSRmat *Neigh, 
                         INT *num_aggregations)
{
    // member of A
    INT row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
    INT *AIA = A->IA;
    INT *AJA = A->JA;
    REAL *Aval = A->val;
    
    // local variable
    REAL strongly_coupled; 
    if (GE(param->tentative_smooth, SMALLREAL)) {
        strongly_coupled= param->strong_coupled * pow(0.5, levelNum-1);
    }
    else {
        strongly_coupled= param->strong_coupled;
    }
    REAL strongly_coupled2 = pow(strongly_coupled,2);
    
    INT i,j,index, row_start, row_end;
    INT *NIA = NULL;
    INT *NJA = NULL;
    REAL *Nval = NULL;
    
    /*------------------------------------------*/
    /* Form strongly coupled neighborhood */
    /*------------------------------------------*/
    dvector diag; 
    fasp_dcsr_getdiag(0, A, &diag);  // get the diagonal entries
    
    fasp_dcsr_alloc(row,col,nnz, Neigh);
    NIA =  Neigh->IA;
    NJA = Neigh->JA;
    Nval = Neigh->val;
    
    //for(i=row+1; i--;) NIA[i] = AIA[i];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, CHUNKSIZE) if(row>OPENMP_HOLDS)
#endif
    for(i=row; i>=0; i--) NIA[i] = AIA[i];
    
    index = 0;
    for (i=0; i<row; ++i) {
        NIA[i] = index;
        row_start = AIA[i]; row_end = AIA[i+1];
        for (j = row_start; j<row_end; ++j) {
            if ((AJA[j] == i) || (pow(Aval[j],2) >= strongly_coupled2 * fabs(diag.val[i]*diag.val[AJA[j]])) ) {
                NJA[index] = AJA[j];
                Nval[index] = Aval[j];
                index++;
            }
        }
    }
    
    NIA[row] = index;
    Neigh->nnz = index;
    
    Neigh->JA = (INT*)fasp_mem_realloc(Neigh->JA, (Neigh->IA[row])*sizeof(INT));
    Neigh->val = (REAL*)fasp_mem_realloc(Neigh->val, (Neigh->IA[row])*sizeof(REAL));
    
    NIA =  Neigh->IA;
    NJA = Neigh->JA;
    Nval = Neigh->val;
    
    fasp_dvec_free(&diag);
    
    /*------------------------------------------*/
    /* Initialization */
    /*------------------------------------------*/
    fasp_ivec_alloc(row, vertices);

    //for (i=row;i--;) vertices->val[i] = -2;
    //for (i=row-1; i>=0; i--) vertices->val[i] = -2;
    fasp_iarray_set(row, vertices->val, -2);
    
    INT num_left = row;
    INT subset;
    INT max_aggregation = param->max_aggregation;
    INT *num_each_aggregation;
    INT count;
    
    *num_aggregations = 0;
    
    /*------------------------------------------*/
    /* Step 1. */
    /*------------------------------------------*/
    for (i=0; i<row; ++i) {
        if ((AIA[i+1] - AIA[i] - 1) == 0) {
            vertices->val[i] = -1;
            num_left--;
        }
        else {
            subset = 1;
            row_start = NIA[i]; row_end = NIA[i+1];
            for (j=row_start; j<row_end; ++j) {
                if (vertices->val[NJA[j]] >= -1) {
                    subset = 0;
                    break;
                }
            }
            if (subset == 1) {
                count = 0;
                vertices->val[i] = *num_aggregations;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for (j=row_start; j<row_end;++j) {
                    if ((NJA[j]!=i) && (count < max_aggregation)){
                        vertices->val[NJA[j]] = *num_aggregations;
                        num_left--;
                        count ++;
                    }
                }
                (*num_aggregations)++;
            }
        }
    }
        
    /*------------------------------------------*/
    /* Step 2. */
    /*------------------------------------------*/
    INT *temp_C = (INT*)fasp_mem_calloc(row,sizeof(INT));
    
    num_each_aggregation = (INT*)fasp_mem_calloc(*num_aggregations,sizeof(INT));
   
    for (i=row;i--;) {
        temp_C[i] = vertices->val[i];  
        if (vertices->val[i] >= 0) {
            num_each_aggregation[vertices->val[i]] ++;
        }
    }
    
    for(i=0; i<row; ++i) {
        if (vertices->val[i] < -1) {
            row_start = NIA[i]; row_end = NIA[i+1];
            for (j=row_start;j<row_end;++j) {
                if (temp_C[NJA[j]] > -1 && num_each_aggregation[temp_C[NJA[j]]] < max_aggregation ) {
                    vertices->val[i] = temp_C[NJA[j]];
                    num_left--;
                    num_each_aggregation[temp_C[NJA[j]]] ++ ;
                    break;
                }
            }
        }
    }
        
    /*------------------------------------------*/
    /* Step 3. */
    /*------------------------------------------*/
    while (num_left > 0) {
        for (i=0; i<row; ++i) {
            if (vertices->val[i] < -1) {
                count = 0;
                vertices->val[i] = *num_aggregations;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for (j=row_start; j<row_end;++j) {
                    if ((NJA[j]!=i) && (vertices->val[NJA[j]] < -1) && (count<max_aggregation) ) {
                        vertices->val[NJA[j]] = *num_aggregations;
                        num_left--;
                        count++;
                    }
                }
                (*num_aggregations)++;
            }
        }
    }
    
    fasp_mem_free(temp_C);
    fasp_mem_free(num_each_aggregation);

}

/**
 * \fn static void form_tentative_p (ivector *vertices, dCSRmat *tentp, AMG_data *mgl, 
 *                                   INT levelNum, INT num_aggregations)
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
 * \date   09/29/2009
 */
static void form_tentative_p (ivector *vertices, 
                              dCSRmat *tentp, 
                              AMG_data *mgl, 
                              INT levelNum, 
                              INT num_aggregations)
{
    INT i, j;
    REAL **basis = mgl->near_kernel_basis;
    
    /* Form tentative prolongation */
    tentp->row = vertices->row;
    tentp->col = num_aggregations;
    tentp->nnz = vertices->row;
    
    tentp->IA  = (INT *)fasp_mem_calloc(tentp->row+1,sizeof(INT));    
    
    // local variables
    INT * IA = tentp->IA;
    INT * JA; 
    REAL *val; 
    INT  *vval = vertices->val;
    
    const INT row = tentp->row;
    
    // first run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            j ++;
        }
    }
    IA[row] = j;
    
    // allocate
    tentp->nnz = j;
    tentp->JA = (INT *)fasp_mem_calloc(tentp->nnz, sizeof(INT));
    tentp->val = (REAL *)fasp_mem_calloc(tentp->nnz, sizeof(REAL));
    
    JA = tentp->JA;
    val = tentp->val;
    
    // second run
    for (i = 0, j = 0; i < row; i ++) {
        IA[i] = j;
        if (vval[i] > -1) {
            JA[j] = vval[i];
            val[j] = basis[0][i];
            j ++;
        }
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
