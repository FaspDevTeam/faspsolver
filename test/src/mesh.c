/*! \file mesh.c
 *  \brief mesh input, output, refine etc.
 */

#include "fasp.h"
#include "mesh.h"

#define DIM 2

/**
 * \fn int mesh_init (Mesh *mesh, const char *filename)
 *
 * \brief Initialize the mesh info from input file
 *
 * \param *mesh          Mesh info
 * \param *filename      Input mesh filename
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/05/2009
 */
int mesh_init (Mesh *mesh, const char *filename)
{
    int num_node, dim_node, num_elem, dim_elem;
    int i,j;
    
    FILE *inputFile=fopen(filename, "r");
    if (inputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }

    // get the node' coordinates
    fscanf(inputFile, "%d %d", &num_node, &dim_node);
    mesh->node.row=num_node;
    mesh->node.col=dim_node;
    
    // alloc mem we need,
    mesh->node.val=(double **)fasp_mem_calloc(num_node, sizeof(double*)); 
    mesh->node.val[0]=(double *)fasp_mem_calloc(num_node*dim_node, sizeof(double)); 
    double *tmp_node = mesh->node.val[0];
    
    for (i=0;i<num_node;++i) { 
        // re-point to the val
        mesh->node.val[i]=&tmp_node[dim_node*i];
        for (j=0;j<dim_node;++j) fscanf(inputFile, "%lf", &mesh->node.val[i][j]);
    }

    // get triangular grid
    fscanf(inputFile, "%d %d", &num_elem, &dim_elem);
    mesh->elem.row=num_elem;
    mesh->elem.col=dim_elem;
    
    // alloc mem we need,
    mesh->elem.val=(int **)fasp_mem_calloc(num_elem, sizeof(int*)); 
    mesh->elem.val[0]=(int *)fasp_mem_calloc(num_elem*dim_elem, sizeof(int)); 
    int *tmp_elem = mesh->elem.val[0];
    
    for (i=0;i<num_elem;++i) {
        // re-point to the val
        mesh->elem.val[i]=&tmp_elem[dim_elem*i];
        for (j=0;j<dim_elem;++j) {
            fscanf(inputFile, "%d", &mesh->elem.val[i][j]);
            mesh->elem.val[i][j]--;
        }
    }
    
    // get node boundary flag
    mesh->node_bd.row = num_node;
    mesh->node_bd.val = (int *)fasp_mem_calloc(num_node, sizeof(int));
    double p[DIM];
    for (i=0;i<num_node;++i) {
        for (j=0;j<dim_node;j++)
            p[j] = mesh->node.val[i][j];
        mesh->node_bd.val[i] = bd_flag(p);
    }

    fclose(inputFile);
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_init_pro (Mesh *mesh, const char *filename)
 *
 * \brief Initialize the mesh info from input file with bd flag
 *
 * \param *mesh          Mesh info
 * \param *filename      Input mesh filename
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/05/2009
 */
int mesh_init_pro (Mesh *mesh, const char *filename)
{
    int num_node, dim_node, num_elem, dim_elem;
    int i,j;
    
    FILE *inputFile=fopen(filename, "r");
    if (inputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }

    // get the node' coordinates
    fscanf(inputFile, "%d %d", &num_node, &dim_node);
    mesh->node.row=num_node;
    mesh->node.col=dim_node;
    
    // alloc mem we need,
    mesh->node.val=(double **)fasp_mem_calloc(num_node, sizeof(double*)); 
    mesh->node.val[0]=(double *)fasp_mem_calloc(num_node*dim_node, sizeof(double)); 
    double *tmp_node = mesh->node.val[0];
    
    for (i=0;i<num_node;++i) { 
        // re-point to the val
        mesh->node.val[i]=&tmp_node[dim_node*i];
        for (j=0;j<dim_node;++j) fscanf(inputFile, "%lf", &mesh->node.val[i][j]);
    }

    // get triangular grid
    fscanf(inputFile, "%d %d", &num_elem, &dim_elem);
    mesh->elem.row=num_elem;
    mesh->elem.col=dim_elem;
    
    // alloc mem we need,
    mesh->elem.val=(int **)fasp_mem_calloc(num_elem, sizeof(int*)); 
    mesh->elem.val[0]=(int *)fasp_mem_calloc(num_elem*dim_elem, sizeof(int)); 
    int *tmp_elem = mesh->elem.val[0];
    
    for (i=0;i<num_elem;++i) {
        // re-point to the val
        mesh->elem.val[i]=&tmp_elem[dim_elem*i];
        for (j=0;j<dim_elem;++j) {
            fscanf(inputFile, "%d", &mesh->elem.val[i][j]);
            mesh->elem.val[i][j]--;
        }
    }
    
    // get node boundary flag
    fscanf(inputFile, "%d", &num_elem);
    mesh->node_bd.row = num_node;
    mesh->node_bd.val = (int *)fasp_mem_calloc(num_node, sizeof(int));
    for (i=0;i<num_node;++i) {
        fscanf(inputFile, "%d", &mesh->node_bd.val[i]);
    }

    fclose(inputFile);
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_aux_init (Mesh *mesh, Mesh_aux *mesh_aux, const char *filename)
 *
 * \brief Initialize the auxiliary mesh info from mesh info and input file
 *
 * \param *mesh          Mesh info
 * \param *mesh_aux      Auxiliary mesh info
 * \param *filename      Input mesh filename
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2009
 */
int mesh_aux_init (Mesh *mesh, Mesh_aux *mesh_aux, const char *filename)
{
    int num_edge, dim_edge, num_elem2edge, dim_elem2edge;
    int i,j;
    
    FILE *inputFile=fopen(filename, "r");
    if (inputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }

    // get the edge info
    fscanf(inputFile, "%d %d", &num_edge, &dim_edge);
    mesh_aux->edge.row=num_edge;
    mesh_aux->edge.col=dim_edge;
    
    // alloc mem we need,
    mesh_aux->edge.val=(int **)fasp_mem_calloc(num_edge, sizeof(int *)); 
    mesh_aux->edge.val[0]=(int *)fasp_mem_calloc(num_edge*dim_edge, sizeof(int)); 
    int *tmp_edge = mesh_aux->edge.val[0];
    
    for (i=0;i<num_edge;++i) {
        // re-point to the val
        mesh_aux->edge.val[i]=&tmp_edge[dim_edge*i];
        for (j=0;j<dim_edge;++j) {
            fscanf(inputFile, "%d", &mesh_aux->edge.val[i][j]);
            mesh_aux->edge.val[i][j]--;
        }
    }

    // get elem's edge info
    fscanf(inputFile, "%d %d", &num_elem2edge, &dim_elem2edge);
    mesh_aux->elem2edge.row=num_elem2edge;
    mesh_aux->elem2edge.col=dim_elem2edge;
    
    // alloc mem we need,
    mesh_aux->elem2edge.val=(int **)fasp_mem_calloc(num_elem2edge, sizeof(int*)); 
    mesh_aux->elem2edge.val[0]=(int *)fasp_mem_calloc(num_elem2edge*dim_elem2edge, sizeof(int)); 
    int *tmp_elem = mesh_aux->elem2edge.val[0];
    
    for (i=0;i<num_elem2edge;++i) {
        // re-point to the val
        mesh_aux->elem2edge.val[i]=&tmp_elem[dim_elem2edge*i];
        for (j=0;j<dim_elem2edge;++j) {
            fscanf(inputFile, "%d", &mesh_aux->elem2edge.val[i][j]);
            mesh_aux->elem2edge.val[i][j]--;
        }
    }
    
    // get edge boundary flag
    mesh_aux->edge_bd.row = num_edge;
    mesh_aux->edge_bd.val = (int *)fasp_mem_calloc(num_edge, sizeof(int));
    double p[DIM];
    for (j=0;j<DIM;++j) p[j] = 0.0;

    int n, k;
    for (i=0;i<num_edge;++i) {
        for (j=0;j<dim_edge;++j) {
            n = mesh_aux->edge.val[i][j];
            for (k=0;k<DIM;++k) {
                p[k] += mesh->node.val[n][k]/DIM;
            }
        }
        mesh_aux->edge_bd.val[i] = bd_flag(p);
        for (j=0;j<DIM;++j) p[j] = 0.0;
    }

    fclose(inputFile);
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_aux_init_pro (Mesh *mesh, Mesh_aux *mesh_aux, const char *filename)
 *
 * \brief Initialize the auxiliary mesh info from mesh info and input file with bd flag
 *
 * \param *mesh          Mesh info
 * \param *mesh_aux      Auxiliary mesh info
 * \param *filename      Input mesh filename
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2009
 */
int mesh_aux_init_pro (Mesh *mesh, Mesh_aux *mesh_aux, const char *filename)
{
    int num_edge, dim_edge, num_elem2edge, dim_elem2edge;
    int i,j;
    
    FILE *inputFile=fopen(filename, "r");
    if (inputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }

    // get the edge info
    fscanf(inputFile, "%d %d", &num_edge, &dim_edge);
    mesh_aux->edge.row=num_edge;
    mesh_aux->edge.col=dim_edge;
    
    // alloc mem we need,
    mesh_aux->edge.val=(int **)fasp_mem_calloc(num_edge, sizeof(int *)); 
    mesh_aux->edge.val[0]=(int *)fasp_mem_calloc(num_edge*dim_edge, sizeof(int)); 
    int *tmp_edge = mesh_aux->edge.val[0];
    
    for (i=0;i<num_edge;++i) {
        // re-point to the val
        mesh_aux->edge.val[i]=&tmp_edge[dim_edge*i];
        for (j=0;j<dim_edge;++j) {
            fscanf(inputFile, "%d", &mesh_aux->edge.val[i][j]);
            mesh_aux->edge.val[i][j]--;
        }
    }

    // get elem's edge info
    fscanf(inputFile, "%d %d", &num_elem2edge, &dim_elem2edge);
    mesh_aux->elem2edge.row=num_elem2edge;
    mesh_aux->elem2edge.col=dim_elem2edge;
    
    // alloc mem we need,
    mesh_aux->elem2edge.val=(int **)fasp_mem_calloc(num_elem2edge, sizeof(int*)); 
    mesh_aux->elem2edge.val[0]=(int *)fasp_mem_calloc(num_elem2edge*dim_elem2edge, sizeof(int)); 
    int *tmp_elem = mesh_aux->elem2edge.val[0];
    
    for (i=0;i<num_elem2edge;++i) {
        // re-point to the val
        mesh_aux->elem2edge.val[i]=&tmp_elem[dim_elem2edge*i];
        for (j=0;j<dim_elem2edge;++j) {
            fscanf(inputFile, "%d", &mesh_aux->elem2edge.val[i][j]);
            mesh_aux->elem2edge.val[i][j]--;
        }
    }
    
    // get edge boundary flag
    fscanf(inputFile, "%d", &num_edge);
    mesh_aux->edge_bd.row = num_edge;
    mesh_aux->edge_bd.val = (int *)fasp_mem_calloc(num_edge, sizeof(int));

    for (i=0;i<num_edge;++i) {
        fscanf(inputFile, "%d", &mesh_aux->edge_bd.val[i]);
    }

    fclose(inputFile);
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_aux_build (Mesh *mesh, Mesh_aux *mesh_aux)
 *
 * \brief Generate auxiliary mesh info
 *
 * \param *mesh          Mesh info
 * \param *mesh_aux      Auxiliary mesh info
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/05/2009
 */
int mesh_aux_build(Mesh *mesh, Mesh_aux *mesh_aux)
{
    int num_node = mesh->node.row;
    int dim_node = mesh->node.col;
    int num_elem = mesh->elem.row;
    int dim_elem = mesh->elem.col;
    
    int num_elem2edge = num_elem;
    int dim_elem2edge = dim_elem;
    mesh_aux->elem2edge.row = num_elem2edge;
    mesh_aux->elem2edge.col = dim_elem2edge;
    int num_edge = 3*num_node;// pre-define num_edge, actually num_edge < 3*num_node
    int dim_edge = dim_node;
    mesh_aux->edge.col = dim_edge;
    
    int i, j, n1, n2, n1t, n2t, count = 0, edge_c;
    double p[DIM];
    mesh_aux->edge.val = (int **)fasp_mem_calloc(num_edge, sizeof(int*));
    mesh_aux->edge.val[0] = (int *)fasp_mem_calloc(num_edge*dim_edge, sizeof(int));
    int *tmp_edge = mesh_aux->edge.val[0];
    mesh_aux->elem2edge.val = (int **)fasp_mem_calloc(num_elem2edge, sizeof(int*));
    mesh_aux->elem2edge.val[0] = (int *)fasp_mem_calloc(num_elem2edge*dim_elem2edge, sizeof(int));
    int *tmp_elem = mesh_aux->elem2edge.val[0];
    mesh_aux->edge_bd.val = (int *)fasp_mem_calloc(num_edge, sizeof(int));

    for (i=0;i<num_edge;++i) {
        mesh_aux->edge.val[i] = &tmp_edge[dim_edge*i];
    }

    for (i=0;i<num_elem2edge;++i) {
        mesh_aux->elem2edge.val[i] = &tmp_elem[dim_elem2edge*i];
    }

    int *edge_aux = NULL;
    if (num_node < 1e3) { // Why having such a constraint? --Chensong
                          // no reason, :-), just to avoid the shortage of memory
                          // see edge_aux following --Feiteng
        edge_aux = (int *)fasp_mem_calloc(num_node*num_node, sizeof(int));
        for (i=0;i<num_elem;++i) {
            for (j=0;j<dim_elem;++j) {
                n1t = mesh->elem.val[i][(j+1)%dim_elem];
                n2t = mesh->elem.val[i][(j+2)%dim_elem];
                n1 = MIN(n1t, n2t);
                n2 = MAX(n1t, n2t);
                edge_c = n1*num_node+n2;
                if (edge_aux[edge_c] == 0) {
                    mesh_aux->edge.val[count][0] = n1;
                    mesh_aux->edge.val[count][1] = n2;
                    mesh_aux->elem2edge.val[i][j] = count;
                    count++;
                    edge_aux[edge_c] = count;
                }
                else if (edge_aux[edge_c] > 0) {
                    mesh_aux->elem2edge.val[i][j] = edge_aux[edge_c]-1;
                    edge_aux[edge_c] *= -1;
                }
                else {
                    printf("### ERROR during the mesh auxiliary building ... \n");
                    exit(ERROR_MISC);
                }
            }// end of loop j
        }// end of loop i
        for (i=0;i<count;++i) {
            n1 = mesh_aux->edge.val[i][0];
            n2 = mesh_aux->edge.val[i][1];
            edge_c = edge_aux[n1*num_node+n2] - 1;
            if (edge_c >= 0) {
                p[0] = (mesh->node.val[n1][0] + mesh->node.val[n2][0])/2;
                p[1] = (mesh->node.val[n1][1] + mesh->node.val[n2][1])/2;
                mesh_aux->edge_bd.val[edge_c] = bd_flag(p);
            }
        }// end of loop i
    }// end of if

    else { // if there is too many node
        int adj_max = 10; // set a max number of adjacent node
        int k = 0;
        int *edge_map= (int *)fasp_mem_calloc(num_node*adj_max, sizeof(int));
        int *adj_count = (int *)fasp_mem_calloc(num_node, sizeof(int));
        edge_aux = (int *)fasp_mem_calloc(num_node*adj_max, sizeof(int));
        for (i=0;i<num_elem;++i) {
            for (j=0;j<dim_elem;++j) {
                n1t = mesh->elem.val[i][(j+1)%dim_elem];
                n2t = mesh->elem.val[i][(j+2)%dim_elem];
                n1 = MIN(n1t, n2t);
                n2 = MAX(n1t, n2t);
                for (k=0;k<adj_count[n1];++k) {
                    if (edge_map[n1*adj_max+k] == n2) {
                        mesh_aux->elem2edge.val[i][j] = edge_aux[n1*adj_max+k];
                        edge_aux[n1*adj_max+k] *= -1;
                        break;
                    }
                }// end of loop k
                if (k == adj_count[n1]) {
                    mesh_aux->edge.val[count][0] = n1;
                    mesh_aux->edge.val[count][1] = n2;
                    mesh_aux->elem2edge.val[i][j] = count;
                    edge_aux[n1*adj_max+k] = count;
                    edge_map[n1*adj_max+k] = n2;
                    adj_count[n1] += 1;
                    count++;
                }// end of if
            }// end of loop j
        }// end of loop i

        for (i=0;i<count;++i) {
            n1 = mesh_aux->edge.val[i][0];
            n2 = mesh_aux->edge.val[i][1];
            for (k=0;k<adj_count[n1];++k) {
                if (edge_map[n1*adj_max+k] == n2) {
                    edge_c = edge_aux[n1*adj_max+k];
                    break;
                }
                if (k == adj_count[n1]) {
                    printf("### ERROR during the mesh auxiliary building ... (case 2) \n");
                    exit(ERROR_MISC);
                }
            }
            if (edge_c >= 0) {
                p[0] = (mesh->node.val[n1][0] + mesh->node.val[n2][0])/2;
                p[1] = (mesh->node.val[n1][1] + mesh->node.val[n2][1])/2;
                mesh_aux->edge_bd.val[edge_c] = bd_flag(p);
            }
        }// end of loop i

        fasp_mem_free(adj_count);
        fasp_mem_free(edge_map);
    }

    num_edge = count;
    mesh_aux->edge.row = num_edge;
    mesh_aux->edge_bd.row = num_edge;
    mesh_aux->edge.val[0] = (int *)fasp_mem_realloc(mesh_aux->edge.val[0], sizeof(int)*num_edge*dim_edge);
    mesh_aux->edge.val = (int **)fasp_mem_realloc(mesh_aux->edge.val, sizeof(int*)*num_edge);
    mesh_aux->edge_bd.val = (int *)fasp_mem_realloc(mesh_aux->edge_bd.val, sizeof(int)*num_edge);
    tmp_edge = mesh_aux->edge.val[0];
    for (i=0;i<num_edge;++i) {
        mesh_aux->edge.val[i] = &tmp_edge[i*dim_edge];
    }
    
    fasp_mem_free(edge_aux);

    return FASP_SUCCESS;
}

/**
 * \fn int mesh_write (Mesh *mesh, const char *filename)
 *
 * \brief Write the mesh information
 *
 * \param *mesh          Mesh info
 * \param *filename      Output mesh file
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2012
 */
int mesh_write (Mesh *mesh, const char *filename)
{
    int num_node=mesh->node.row;
    int dim_node=mesh->node.col;
    int num_elem=mesh->elem.row;
    int dim_elem=mesh->elem.col;
    int i,j;
    
    FILE *outputFile=fopen(filename, "w");
    if (outputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }
    
    fprintf(outputFile, "%d %d\n", num_node, dim_node);
    
    for (i=0;i<num_node;++i) {
        for (j=0;j<dim_node;++j) {
            fprintf(outputFile, "%lf ", mesh->node.val[i][j]);
        }
        fprintf(outputFile, "\n");
    }
    
    fprintf(outputFile, "%d %d\n", num_elem, dim_elem);
    
    for (i=0;i<num_elem;++i) {
        for (j=0;j<dim_elem;++j) {
            fprintf(outputFile, "%d ", mesh->elem.val[i][j]+1);
        }
        fprintf(outputFile, "\n");
    }
    
    fclose(outputFile);
    
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_write_pro (Mesh *mesh, const char *filename)
 *
 * \brief Write the mesh information with bd flag
 *
 * \param *mesh          Mesh info
 * \param *filename      Output mesh file
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2012
 */
int mesh_write_pro (Mesh *mesh, const char *filename)
{
    int num_node=mesh->node.row;
    int dim_node=mesh->node.col;
    int num_elem=mesh->elem.row;
    int dim_elem=mesh->elem.col;
    int i,j;
    
    FILE *outputFile=fopen(filename, "w");
    if (outputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }
    
    fprintf(outputFile, "%d %d\n", num_node, dim_node);
    
    for (i=0;i<num_node;++i) {
        for (j=0;j<dim_node;++j) {
            fprintf(outputFile, "%lf ", mesh->node.val[i][j]);
        }
        fprintf(outputFile, "\n");
    }
    
    fprintf(outputFile, "%d %d\n", num_elem, dim_elem);
    
    for (i=0;i<num_elem;++i) {
        for (j=0;j<dim_elem;++j) {
            fprintf(outputFile, "%d ", mesh->elem.val[i][j]+1);
        }
        fprintf(outputFile, "\n");
    }

    fprintf(outputFile, "%d\n", num_node);
    for (i=0;i<num_node;++i) {
        fprintf(outputFile, "%d\n", mesh->node_bd.val[i]);
    }
    
    fclose(outputFile);
    
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_aux_write (Mesh_aux *mesh_aux, const char *filename)
 *
 * \brief Write the auxiliary mesh information
 *
 * \param *mesh_aux      Auxiliary mesh info
 * \param *filename      Output mesh file
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2012
 */
int mesh_aux_write (Mesh_aux *mesh_aux, const char *filename)
{
    int num_edge=mesh_aux->edge.row;
    int dim_edge=mesh_aux->edge.col;
    int num_elem2edge=mesh_aux->elem2edge.row;
    int dim_elem2edge=mesh_aux->elem2edge.col;
    int i,j;
    
    FILE *outputFile=fopen(filename, "w");
    if (outputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }
    
    fprintf(outputFile, "%d %d\n", num_edge, dim_edge);
    
    for (i=0;i<num_edge;++i) {
        for (j=0;j<dim_edge;++j) {
            fprintf(outputFile, "%d ", mesh_aux->edge.val[i][j]+1);
        }
        fprintf(outputFile, "\n");
    }
    
    fprintf(outputFile, "%d %d\n", num_elem2edge, dim_elem2edge);
    
    for (i=0;i<num_elem2edge;++i) {
        for (j=0;j<dim_elem2edge;++j) {
            fprintf(outputFile, "%d ", mesh_aux->elem2edge.val[i][j]+1);
        }
        fprintf(outputFile, "\n");
    }
    
    fclose(outputFile);
    
    return FASP_SUCCESS;
    
}

/**
 * \fn int mesh_aux_write_pro (Mesh_aux *mesh_aux, const char *filename)
 *
 * \brief Write the auxiliary mesh information with bd flag
 *
 * \param *mesh_aux      Auxiliary mesh info
 * \param *filename      Output mesh file
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2012
 */
int mesh_aux_write_pro (Mesh_aux *mesh_aux, const char *filename)
{
    int num_edge=mesh_aux->edge.row;
    int dim_edge=mesh_aux->edge.col;
    int num_elem2edge=mesh_aux->elem2edge.row;
    int dim_elem2edge=mesh_aux->elem2edge.col;
    int i,j;
    
    FILE *outputFile=fopen(filename, "w");
    if (outputFile==NULL) {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(ERROR_OPEN_FILE);
    }
    
    fprintf(outputFile, "%d %d\n", num_edge, dim_edge);
    
    for (i=0;i<num_edge;++i) {
        for (j=0;j<dim_edge;++j) {
            fprintf(outputFile, "%d ", mesh_aux->edge.val[i][j]+1);
        }
        fprintf(outputFile, "\n");
    }
    
    fprintf(outputFile, "%d %d\n", num_elem2edge, dim_elem2edge);
    
    for (i=0;i<num_elem2edge;++i) {
        for (j=0;j<dim_elem2edge;++j) {
            fprintf(outputFile, "%d ", mesh_aux->elem2edge.val[i][j]+1);
        }
        fprintf(outputFile, "\n");
    }

    fprintf(outputFile, "%d\n", num_edge);

    for (i=0;i<num_edge;++i) {
        fprintf(outputFile, "%d\n", mesh_aux->edge_bd.val[i]);
    }
    
    fclose(outputFile);
    
    return FASP_SUCCESS;
    
}

/**
 * \fn int mesh_free (Mesh *mesh)
 *
 * \brief free memory for mesh info
 *
 * \param *mesh          Mesh info
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/05/2009
 */
int mesh_free (Mesh *mesh)
{
    fasp_mem_free(mesh->node.val[0]);
    fasp_mem_free(mesh->node.val);
    fasp_mem_free(mesh->elem.val[0]);
    fasp_mem_free(mesh->elem.val);
    fasp_mem_free(mesh->node_bd.val);

    mesh->node.row = 0;
    mesh->node.col = 0;
    mesh->elem.row = 0;
    mesh->elem.col = 0;
    mesh->node_bd.row = 0;
    mesh->node.val = NULL;
    mesh->elem.val = NULL;
    mesh->node_bd.val  = NULL;
    
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_aux_free (Mesh_aux *mesh_aux)
 *
 * \brief free memory for auxiliary mesh info
 *
 * \param *mesh_aux          auxiliary mesh info
 *
 * \return               1 if succeed 0 if fail
 *
 * \author Feiteng Huang
 * \date   04/05/2009
 */
int mesh_aux_free (Mesh_aux *mesh_aux)
{
    fasp_mem_free(mesh_aux->edge.val[0]);
    fasp_mem_free(mesh_aux->edge.val);
    fasp_mem_free(mesh_aux->elem2edge.val[0]);
    fasp_mem_free(mesh_aux->elem2edge.val);
    fasp_mem_free(mesh_aux->edge_bd.val);

    mesh_aux->edge.row = 0;
    mesh_aux->edge.col = 0;
    mesh_aux->elem2edge.row = 0;
    mesh_aux->elem2edge.col = 0;
    mesh_aux->edge_bd.row = 0;
    mesh_aux->edge.val = NULL;
    mesh_aux->elem2edge.val = NULL;
    mesh_aux->edge_bd.val = NULL;
    
    return FASP_SUCCESS;
}

/**
 * \fn int mesh_refine (Mesh *mesh, Mesh_aux *mesh_aux)
 *
 * \brief refine mesh use mesh info and auxiliary mesh info
 *
 * \param *mesh          Mesh info
 * \param *mesh_aux      Auxiliary mesh info
 *
 * \return               FASP_SUCCESS if succeed
 *
 * \author Feiteng Huang
 * \date   04/06/2009
 */
int mesh_refine(Mesh *mesh, Mesh_aux *mesh_aux)
{
    int num_node = mesh->node.row;
    int dim_node = mesh->node.col;
    int num_elem = mesh->elem.row;
    int dim_elem = mesh->elem.col;
    int num_edge = mesh_aux->edge.row;
    int dim_edge = mesh_aux->edge.col;
    
    int num_newnode = num_node + num_edge;
    int num_newedge = 3*num_elem + 2*num_edge;
    int num_newelem = 4*num_elem;
    
    int i, j, k;
    int n[3], e[6];
    
    mesh->node.row = num_newnode;
    mesh->elem.row = num_newelem;
    mesh->node_bd.row = num_newnode;
    mesh_aux->edge.row = num_newedge;
    mesh_aux->elem2edge.row = num_newelem;
    mesh_aux->edge_bd.row = num_newedge;
    
    // realloc node, elem, edge, elem2edge, node_bd, edge_bd
    mesh->node.val[0] = (double *)fasp_mem_realloc(mesh->node.val[0], sizeof(double)*num_newnode*dim_node);
    mesh->node.val = (double **)fasp_mem_realloc(mesh->node.val, sizeof(double *)*num_newnode);
    double *tmp_node = mesh->node.val[0];
    mesh->elem.val[0] = (int *)fasp_mem_realloc(mesh->elem.val[0], sizeof(int)*num_newelem*dim_elem);
    mesh->elem.val = (int **)fasp_mem_realloc(mesh->elem.val, sizeof(int *)*num_newelem);
    int *tmp_elem = mesh->elem.val[0];
    mesh_aux->edge.val[0] = (int *)fasp_mem_realloc(mesh_aux->edge.val[0], sizeof(int)*num_newedge*dim_edge);
    mesh_aux->edge.val = (int **)fasp_mem_realloc(mesh_aux->edge.val, sizeof(int *)*num_newedge);
    int *tmp_edge = mesh_aux->edge.val[0];
    mesh_aux->elem2edge.val[0] = (int *)fasp_mem_realloc(mesh_aux->elem2edge.val[0], sizeof(int)*num_newelem*dim_elem);
    mesh_aux->elem2edge.val = (int **)fasp_mem_realloc(mesh_aux->elem2edge.val, sizeof(int *)*num_newelem);
    int *tmp_elem2edge = mesh_aux->elem2edge.val[0];
    mesh->node_bd.val = (int *)fasp_mem_realloc(mesh->node_bd.val, sizeof(int)*num_newnode);
    mesh_aux->edge_bd.val = (int *)fasp_mem_realloc(mesh_aux->edge_bd.val, sizeof(int)*num_newedge);

    for (i=0;i<num_newnode;++i) {
        mesh->node.val[i] = &tmp_node[i*dim_node];
    }
    for (i=0;i<num_newedge;++i) {
        mesh_aux->edge.val[i] = &tmp_edge[i*dim_edge];
    }
    for (i=0;i<num_newelem;++i) {
        mesh->elem.val[i] = &tmp_elem[i*dim_elem];
        mesh_aux->elem2edge.val[i] = &tmp_elem2edge[i*dim_elem];
    }
    
    // update mesh & mesh_aux info
    for (i=0;i<num_edge;++i) {

        // update node info
        
        // init node value
        for (k=0;k<dim_node;++k) {
            mesh->node.val[i+num_node][k] = 0;
        }
        for (j=0;j<dim_edge;++j) {
            n[0] = mesh_aux->edge.val[i][j];
            for (k=0;k<dim_node;++k) {
                mesh->node.val[i+num_node][k] += mesh->node.val[n[0]][k]/dim_edge;
            }
        }

        // update node_bd info
        mesh->node_bd.val[i+num_node] = mesh_aux->edge_bd.val[i];
        
        // update edge & edge_bd on original edge
        n[0] = num_node + i;
        n[1] = mesh_aux->edge.val[i][1];

        // update auxiliary mesh info
        mesh_aux->edge.val[i][1] = n[0];
        mesh_aux->edge.val[i+num_edge][0] = n[0];
        mesh_aux->edge.val[i+num_edge][1] = n[1];
        
        mesh_aux->edge_bd.val[i+num_edge] = mesh_aux->edge_bd.val[i];
    }
    
    for (i=0;i<num_elem;++i) {

        for (j=0;j<dim_elem;++j) {
            n[j] = mesh_aux->elem2edge.val[i][j] + num_node;
        }

        for (j=0;j<dim_elem;++j) {
            // update edge info on original elem
            mesh_aux->edge.val[2*num_edge+3*i+j][0] = n[(j+1)%dim_elem];
            mesh_aux->edge.val[2*num_edge+3*i+j][1] = n[(j+2)%dim_elem];
            mesh_aux->edge_bd.val[2*num_edge+3*i+j] = INTERIORI;
            if (mesh->elem.val[i][(j+1)%dim_elem] == mesh_aux->edge.val[mesh_aux->elem2edge.val[i][j]][0]) {
                e[2*j] = mesh_aux->elem2edge.val[i][j];
                e[2*j+1] = mesh_aux->elem2edge.val[i][j] + num_edge;
            }
            else {
                e[2*j+1] = mesh_aux->elem2edge.val[i][j];
                e[2*j] = mesh_aux->elem2edge.val[i][j] + num_edge;
            }
            
            // update elem & elem2edge info
            mesh->elem.val[num_elem+3*i+2][j] = n[j];
            mesh_aux->elem2edge.val[num_elem+3*i+2][j] = 2*num_edge+3*i+j;
        }
        
        // update elem info
        mesh->elem.val[num_elem+3*i][0] = mesh->elem.val[i][1];
        mesh->elem.val[num_elem+3*i][1] = n[0];
        mesh->elem.val[num_elem+3*i][2] = n[2];
        
        mesh->elem.val[num_elem+3*i+1][0] = mesh->elem.val[i][2];
        mesh->elem.val[num_elem+3*i+1][1] = n[1];
        mesh->elem.val[num_elem+3*i+1][2] = n[0];
        
        mesh->elem.val[i][1] = n[2];
        mesh->elem.val[i][2] = n[1];
        
        // update elem2edge info
        mesh_aux->elem2edge.val[num_elem+3*i][0] = 2*num_edge+3*i+1;
        mesh_aux->elem2edge.val[num_elem+3*i][1] = e[5];
        mesh_aux->elem2edge.val[num_elem+3*i][2] = e[0];
        
        mesh_aux->elem2edge.val[num_elem+3*i+1][0] = 2*num_edge+3*i+2;
        mesh_aux->elem2edge.val[num_elem+3*i+1][1] = e[1];
        mesh_aux->elem2edge.val[num_elem+3*i+1][2] = e[2];
        
        mesh_aux->elem2edge.val[i][0] = 2*num_edge+3*i;
        mesh_aux->elem2edge.val[i][1] = e[3];
        mesh_aux->elem2edge.val[i][2] = e[4];
    }
    
    return FASP_SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/

