/*! \file  aggregation_csr.inl
 *  \brief Utilies for multigrid cycles for CSR matrices
 */

#ifdef _OPENMP
#include <omp.h>
#endif

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
 * Modified by Chunsheng Feng, Zheng Li on 09/03/2012
 */
static void aggregation (dCSRmat *A,
                         ivector *vertices, 
                         AMG_param *param, 
                         INT levelNum, 
                         dCSRmat *Neigh, 
                         INT *num_aggregations)
{
    // member of A
    INT    row = A->row, col = A->col, nnz = A->IA[row]-A->IA[0];
    INT  * AIA = A->IA, * AJA = A->JA;
    REAL * Aval = A->val;
    
    // local variables
    REAL strongly_coupled; 
    if ( GE(param->tentative_smooth, SMALLREAL) ) {
        strongly_coupled = param->strong_coupled * pow(0.5, levelNum-1);
    }
    else {
        strongly_coupled = param->strong_coupled;
    }
    REAL strongly_coupled2 = pow(strongly_coupled,2);
    
    INT i, j, index, row_start, row_end;
    INT  * NIA = NULL;
    INT  * NJA = NULL;
    REAL * Nval = NULL;
    
    /*------------------------------------------*/
    /*    Form strongly coupled neighborhood    */
    /*------------------------------------------*/
    dvector diag; 
    fasp_dcsr_getdiag(0, A, &diag);  // get the diagonal entries
        
    fasp_dcsr_alloc(row, col, nnz, Neigh);
    
    NIA  = Neigh->IA;
    NJA  = Neigh->JA;
    Nval = Neigh->val;
    
#ifdef _OPENMP
#pragma omp parallel for if(row>OPENMP_HOLDS)
#endif
    for ( i = row; i >= 0; i-- ) NIA[i] = AIA[i];
    
    for ( index = i = 0; i < row; ++i ) {
        NIA[i] = index;
        row_start = AIA[i]; row_end = AIA[i+1];
        for ( j = row_start; j < row_end; ++j ) {
            if ( (AJA[j] == i) || (pow(Aval[j],2) >= strongly_coupled2 * fabs(diag.val[i]*diag.val[AJA[j]])) ) {
                NJA[index] = AJA[j];
                Nval[index] = Aval[j];
                index++;
            }
        }
    }    
    NIA[row] = index;

    Neigh->nnz = index;
    Neigh->JA  = (INT*) fasp_mem_realloc(Neigh->JA,  (Neigh->IA[row])*sizeof(INT));
    Neigh->val = (REAL*)fasp_mem_realloc(Neigh->val, (Neigh->IA[row])*sizeof(REAL));
    
    NIA  = Neigh->IA;
    NJA  = Neigh->JA;
    Nval = Neigh->val;
    
    fasp_dvec_free(&diag);
    
    /*------------------------------------------*/
    /*             Initialization               */
    /*------------------------------------------*/
    fasp_ivec_alloc(row, vertices);
    fasp_iarray_set(row, vertices->val, -2);
    
    INT num_left = row;
    INT subset, count;
    INT max_aggregation = param->max_aggregation;
    INT *num_each_aggregation;
    
    *num_aggregations = 0;
    
    /*-------------*/
    /*   Step 1.   */
    /*-------------*/
    for ( i = 0; i < row; ++i ) {
        if ( (AIA[i+1] - AIA[i]) == 1 ) {
            vertices->val[i] = -1;
            num_left--;
        }
        else {
            subset = TRUE;
            row_start = NIA[i]; row_end = NIA[i+1];
            for ( j = row_start; j < row_end; ++j ) {
                if ( vertices->val[NJA[j]] >= -1 ) {
                    subset = FALSE;
                    break;
                }
            }
            if ( subset ) {
                count = 0;
                vertices->val[i] = *num_aggregations;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for ( j = row_start; j < row_end; ++j ) {
                    if ( (NJA[j]!=i) && (count < max_aggregation) ) {
                        vertices->val[NJA[j]] = *num_aggregations;
                        num_left--;
                        count ++;
                    }
                }
                (*num_aggregations)++;
            }
        }
    }
        
    /*-------------*/
    /*   Step 2.   */
    /*-------------*/
    INT *temp_C = (INT*)fasp_mem_calloc(row,sizeof(INT));

    if (*num_aggregations == 0){
        printf("WWARNING -- did not find any aggregate in the first round, the matrix might be diagonal matrix\n");
    }
    
    num_each_aggregation = (INT*)fasp_mem_calloc(*num_aggregations,sizeof(INT));
   
    for ( i = row; i--; ) {
        temp_C[i] = vertices->val[i];  
        if ( vertices->val[i] >= 0 ) num_each_aggregation[vertices->val[i]] ++;
    }
    
    for ( i = 0; i < row; ++i ) {
        if ( vertices->val[i] < -1 ) {
            row_start = NIA[i]; row_end = NIA[i+1];
            for ( j = row_start; j < row_end; ++j ) {
                if (temp_C[NJA[j]] > -1 && num_each_aggregation[temp_C[NJA[j]]] < max_aggregation ) {
                    vertices->val[i] = temp_C[NJA[j]];
                    num_left--;
                    num_each_aggregation[temp_C[NJA[j]]] ++ ;
                    break;
                }
            }
        }
    }
        
    /*-------------*/
    /*   Step 3.   */
    /*-------------*/
    while ( num_left > 0 ) {
        for ( i = 0; i < row; ++i ) {
            if ( vertices->val[i] < -1 ) {
                count = 0;
                vertices->val[i] = *num_aggregations;
                num_left--;
                count++;
                row_start = NIA[i]; row_end = NIA[i+1];
                for ( j = row_start; j < row_end; ++j ) {
                    if ( (NJA[j]!=i) && (vertices->val[NJA[j]] < -1) && (count<max_aggregation) ) {
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
    INT  *IA = tentp->IA;
    INT  *vval = vertices->val;    
    const INT row = tentp->row;
    
    // first run
    for ( i = 0, j = 0; i < row; i++ ) {
        IA[i] = j;
        if (vval[i] > -1) j++;
    }
    IA[row] = j;
    
    // allocate memory for P
    tentp->nnz = j;
    tentp->JA  = (INT *)fasp_mem_calloc(tentp->nnz, sizeof(INT));
    tentp->val = (REAL *)fasp_mem_calloc(tentp->nnz, sizeof(REAL));
    
    INT  *JA = tentp->JA;
    REAL *val = tentp->val;
    
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

/**
 * \fn static void pairwise_aggregation (dCSRmat *A,
 *                                       INT pair,
 *                                       ivector *vertices, 
 *                                       INT *num_aggregations)
 *
 * \brief Form aggregation based on pairwise matching 
 *
 * \param A                 Pointer to the coefficient matrices
 * \param pair              Number of pairs in matching
 * \param vertices          Pointer to the aggregation of vertics
 * \param num_aggregations  Pointer to number of aggregations 
 * 
 * \author Xiaoping Li, Zheng Li, Chensong Zhang
 * \date   04/21/2014
 */
static void pairwise_aggregation (const dCSRmat * A,
                                  INT pair,
                                  ivector *vertices,
                                  INT *num_aggregations)
{
    INT row  = A->row;
    INT *AIA = A->IA;
    INT *AJA = A->JA;
    REAL *Aval = A->val;
    INT i, j, row_start, row_end;
    REAL sum;

    REAL k_tg = 6.75;

    /*---------------------------------------------------------- */
    /* Step 1.select very strong diagnoal dominate row (vertices)*/ 
	/*        and store in G0.                                   */
    /*-----------------------------------------------------------*/

    /* G0 : vertices->val[i] =-5 Remain: vertices->val[i] =-1 */

    fasp_ivec_alloc(row, vertices);
    	
    if ( pair == 1 ) {
        for ( i = 0; i < row; i++ ) {
            sum = 0.0;
            row_start = AIA[i]; 
            row_end = AIA[i+1];

            for ( j = row_start+1; j < row_end; j++) sum += ABS(Aval[j]);
        
            if ( Aval[AIA[i]] >= ((k_tg+1.)/(k_tg-1.))*sum) {
                vertices->val[i] = -5;
            }
            else {
                vertices->val[i] = -1;
            }
        }
	}
	else {
		fasp_iarray_set(row, vertices->val, -1);
	}
    
    /*-------------------------------------------------------*/
    /* Step 2. computate row sum (off-diagonal) for each vertice*/
    /*-------------------------------------------------------*/
    REAL *s;

    s = (REAL *)fasp_mem_calloc(row, sizeof(REAL));

    for ( i = 0; i < row; i++ ) {
        s[i] = 0.0;
        row_start = AIA[i];
		row_end = AIA[i+1];

        if (vertices->val[i] == -5) continue;

        for ( j = row_start + 1; j < row_end; j++ ) s[i] -= Aval[j];
    }

    /*------------------------------------------*/
    /* Step 3. start the pairwise aggregation   */
    /*------------------------------------------*/
    REAL mu;
    REAL min_mu;
    REAL aii;
    REAL ajj;
    REAL aij;
    INT  col,index;

    /* temp var for comput mu */
    REAL temp1, temp2, temp3, temp4;
    
    *num_aggregations = 0;

    for ( i = 0; i < row; i++ ) {
        if ( vertices->val[i] != -1) continue;
        
        aij = 0.0;
        min_mu = 1000.0;
        
        row_start = AIA[i]; row_end = AIA[i+1];
        
        aii = Aval[row_start];
        
        for ( j= row_start + 1; j < row_end; j++ ) {
            col = AJA[j];
            if ( vertices->val[col] != -1 ) continue;
            
            aij = Aval[j];
            ajj = Aval[AIA[col]];
            
            temp1 = aii+s[i]+2*aij; 
            temp2 = ajj+s[col]+2*aij;
            temp2 = 1.0/temp1+1.0/temp2;
 
			temp3 = MAX(aii-s[i], SMALLREAL); // avoid temp3 to be zero		
			temp4 = MAX(ajj-s[col], SMALLREAL);	// avoid temp4 to be zero
            temp4 = -aij+1./(1.0/temp3+1.0/temp4); 
            
            mu    = (-aij+1.0/temp2) / temp4;

            if ( min_mu > mu ) {
                min_mu = mu;
                index   = col;
            }
        }

        *num_aggregations += 1;
        if ( min_mu <= k_tg ) {
            vertices->val[i] = *num_aggregations - 1;
            vertices->val[index] = *num_aggregations - 1;
        }
		else {
            vertices->val[i] = *num_aggregations - 1;
        }
    }

    fasp_mem_free(s);
}


/**
 * \fn static void aggregation_coarsening (AMG_data *mgl,
 *                                         AMG_param *param,
 *                                         INT level,
 *                                         ivector *vertice, 
 *                                         INT *num_aggregations)
 *
 * \brief AMG coarsening based on pairwise matching aggregation 
 *
 * \param mgl               Pointer to AMG levele data
 * \param param             Pointer to AMG parameters
 * \param level             Level number
 * \param vertices          Pointer to the aggregation of vertics
 * \param num_aggregations  Pointer to number of aggregations 
 * 
 * \author Xiaoping Li, Zheng Li, Chensong Zhang
 * \date   04/21/2014
 */
static void aggregation_coarsening (AMG_data *mgl, 
                                    AMG_param *param, 
                                    INT level, 
                                    ivector *vertice, 
                                    INT *num_aggregations)
{
    INT i, j, num_agg, aggindex;
    INT pair_number = param->pair_number;

    dCSRmat tmpA = mgl[level].A;

    for ( i = 0; i < pair_number; ++i ) {
        /*-- generate aggregations by pairwise matching --*/
        pairwise_aggregation(&tmpA, i+1, &vertice[level+i], &num_agg);

        if ( i < pair_number-1 ) {
            /*-- Form Prolongation --*/
            form_tentative_p(&vertice[level+i], &mgl[level].P, &mgl[0], level+1, num_agg);
        
            /*-- Perform aggressive coarsening only up to the specified level --*/
            if ( mgl[level].P.col < 20 ) break; // If coarse < 20, stop!!!
        
            /*-- Form resitriction --*/
            fasp_dcsr_trans(&mgl[level].P, &mgl[level].R);
        
            /*-- Form coarse level stiffness matrix --*/
            fasp_blas_dcsr_rap_agg(&mgl[level].R, &tmpA, &mgl[level].P, &mgl[level+1].A);

		    tmpA = mgl[level+1].A;

		    fasp_dcsr_free(&mgl[level].P);
		    fasp_dcsr_free(&mgl[level].R);
		}
	}

    // global aggregation index 
	if ( pair_number > 1 ) {
		for ( i = 0; i < mgl[level].A.row; ++i ) {
			aggindex = vertice[level].val[i];

			if ( aggindex < 0) continue;

			for ( j = 1; j < pair_number; ++j ) aggindex = vertice[level+j].val[aggindex];

			vertice[level].val[i] = aggindex;
		}
	}

   *num_aggregations = num_agg;
  
   /*-- clean memory --*/ 
   fasp_dcsr_free(&mgl[level+1].A);
   for ( i = 1; i < pair_number; ++i ) fasp_ivec_free(&vertice[level+i]); 
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
