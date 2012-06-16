/*! \file poisson_fem.c
 *  \brief Setup P1 FEM for the Poisson's equation
 */

#include <math.h>

#include "poisson_fem.h"
#include "basis.inl"

#define DIM 2          // dimenision of space
#define MAX_QUAD 49    // max num of quadrature points allowed

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/** 
 * \fn static void localb (double (*node)[DIM],double *b, int num_qp)
 *
 * \brief Form local right-hand side b from triangle node
 *
 * \param (*node)[DIM]   Vertice of the triangule
 * \param *b             Local right-hand side
 * \param num_qp         Number of quad points 
 *
 * \author Xuehai Huang and Feiteng Huang
 * \date   03/29/2009
 * 
 * Modified by Feiteng Huang 02/23/2012: User specified number of quadrature points
 * Modified by Chensong on 05/03/2012
 */
static void localb (double (*node)[DIM],
                    double *b,
                    int num_qp)
{
    double gauss[MAX_QUAD][DIM+1];
    double a,p[DIM];
    int i,j;

#if DIM==2 // 2D case

    const double s=2.0*areaT(node[0][0],node[1][0],node[2][0],
                             node[0][1],node[1][1],node[2][1]);
    
    fasp_gauss2d(num_qp, 2, gauss); // gauss intergation initial    
    
    for (i=0;i<=DIM;++i) b[i]=0; // initialize b
    
    for (i=0;i<num_qp;++i)
        {
            for (j=0;j<DIM;++j)
                p[j]=node[0][j]*gauss[i][0]
                    +node[1][j]*gauss[i][1]
                    +node[2][j]*(1-gauss[i][0]-gauss[i][1]);
        
            a=f(p);
            b[0]+=a*gauss[i][2]*gauss[i][0];
            b[1]+=a*gauss[i][2]*gauss[i][1];
            b[2]+=a*gauss[i][2]*(1-gauss[i][0]-gauss[i][1]);        
        }
    
    b[0]*=s; b[1]*=s; b[2]*=s;

#else

    printf("### ERROR: DIM=%d is not supported currently!\n", DIM);
    exit(RUN_FAIL);

#endif

    return;
}

/**
 * \fn static void assemble_stiffmat (dCSRmat *A, dvector *b, Mesh *mesh, 
 *                                    Mesh_aux *mesh_aux, FEM_param *pt )
 *
 * \brief Assemble stiffness matrix and right-hand side
 *
 * \param *A           Pointer to stiffness matrix
 * \param *b           Pointer to right hand side
 * \param *mesh        Pointer to mesh info 
 * \param *mesh_aux    Pointer to auxiliary mesh info 
 * \param *pt          Pointer to parameter 
 *
 * \author Xuehai Huang and Feiteng Huang
 * \date   03/29/2009
 *
 * \note This subroutine follows Ludmil' note on CSR
 *   
 * Modified by Feiteng Huang on 04/06/2012: restructure assmebling
 */
static void assemble_stiffmat (dCSRmat *A, 
                               dvector *b, 
                               Mesh *mesh, 
                               Mesh_aux *mesh_aux, 
                               FEM_param *pt)
{
    /* yangkai: use the information of edgesTran to get the CSR structure of A directly, 
     * because, for triangle elements, A_ij is nonzero if and only if node i and j are 
     * neighbours in the mesh.
     *
     * Feiteng: mesh_aux->edge do the same thing. :-)
     */
    
    // assemble option
    int mat = 0, rhs = 0;
    if ( strcmp(pt->option,"ab") == 0 ) { mat = 1; rhs = 1; }
    if ( strcmp(pt->option,"a") == 0 ) mat = 1;
    if ( strcmp(pt->option,"b") == 0 ) rhs = 1;
    
    const int num_node = mesh->node.row;
    const int num_edge = mesh_aux->edge.row;
    const int nnz = num_node + 2*num_edge;
    const double epsilon=1;
    
    double T[3][2],phi1[2],phi2[2];
    double gauss[MAX_QUAD][DIM+1];
    double s;
    
    int i,j,k;
    int k1,n1,n2,i1;
    double btmp[3], tmp_a;
    int tmp, edge_c;
    
    fasp_gauss2d(pt->num_qp_mat, 2, gauss); // Gauss intergation initial
    
    // alloc mem for A & b
    A->row = A->col = num_node;
    A->nnz = nnz;
    A->IA = (int*)fasp_mem_calloc(A->row+1, sizeof(int));
    A->JA = (int*)fasp_mem_calloc(nnz, sizeof(int));
    A->val = (double*)fasp_mem_calloc(nnz, sizeof(double));
    b->row = num_node;
    b->val = (double*)fasp_mem_calloc(num_node, sizeof(double));
    total_alloc_mem += (A->row+1)*sizeof(int);
    total_alloc_mem += nnz*sizeof(int);
    total_alloc_mem += nnz*sizeof(double);
    total_alloc_mem += num_node*sizeof(double);
    int *count = (int*)fasp_mem_calloc(num_node, sizeof(int));
    int *edge2idx_g1 = (int*)fasp_mem_calloc(num_edge, sizeof(int));
    int *edge2idx_g2 = (int*)fasp_mem_calloc(num_edge, sizeof(int));//edge to global idx of A->val
    
    // get IA
    for (i=0;i<num_edge;++i) {
        n1 = mesh_aux->edge.val[i][0];
        n2 = mesh_aux->edge.val[i][1];
        A->IA[n1+1] += 1;
        A->IA[n2+1] += 1;
    }
    for (i=0;i<num_node;++i) {
        A->IA[i+1] += A->IA[i] + 1;
        A->JA[A->IA[i]] = i;
    }
    
    // get JA 
    for (i=0;i<num_edge;++i) {
        n1 = mesh_aux->edge.val[i][0];
        n2 = mesh_aux->edge.val[i][1];
        count[n1]++;
        count[n2]++;
        A->JA[A->IA[n1]+count[n1]] = n2;
        A->JA[A->IA[n2]+count[n2]] = n1;
        edge2idx_g1[i] = A->IA[n1]+count[n1];
        edge2idx_g2[i] = A->IA[n2]+count[n2];
    }
    fasp_mem_free(count);
    
    // Loop element by element and compute the actual entries storing them in A
    if (rhs && mat) {
        for (k=0;k<mesh->elem.row;++k) { 

            for (i=0;i<mesh->elem.col;++i) {
                j=mesh->elem.val[k][i];
                T[i][0]=mesh->node.val[j][0];
                T[i][1]=mesh->node.val[j][1];
            } // end for i

            s=areaT(T[0][0], T[1][0], T[2][0], T[0][1], T[1][1], T[2][1]);
            localb(T,btmp,pt->num_qp_rhs);
            tmp=0;
            
            // mesh->elem.col == mesh_aux->elem2edge.col
            // "if" should be out of "for" --Chensong
        // fixed? --Feiteng
            for (k1=0;k1<mesh->elem.col;++k1) {

                edge_c = mesh_aux->elem2edge.val[k][k1];
                i=mesh->elem.val[k][k1];

                b->val[i]+=btmp[tmp++];

                // for the diag entry
                gradBasisP1(T, s, k1, phi1);
                tmp_a = 2*s*(phi1[0]*phi1[0]+phi1[1]*phi1[1]*epsilon);
                for (i1=0;i1<pt->num_qp_mat;i1++) {
                    A->val[A->IA[i]] += gauss[i1][2]*tmp_a;
                } // end for i1

                // for the off-diag entry
                i = edge2idx_g1[edge_c];
                j = edge2idx_g2[edge_c];
                n1 = (k1+1)%3;//local node index in the elem
                n2 = (k1+2)%3;
                gradBasisP1(T, s, n1, phi1);
                gradBasisP1(T, s, n2, phi2);
                tmp_a = 2*s*(phi1[0]*phi2[0]+phi1[1]*phi2[1]*epsilon);
                for (i1=0;i1<pt->num_qp_mat;i1++) {
                    A->val[i] += gauss[i1][2]*tmp_a;
                    A->val[j] += gauss[i1][2]*tmp_a; 
                } // end for i1
            } // end for k1
        } // end for k
    }
    else if (rhs) {
        for (k=0;k<mesh->elem.row;++k) { 

            for (i=0;i<mesh->elem.col;++i) {
                j=mesh->elem.val[k][i];
                T[i][0]=mesh->node.val[j][0];
                T[i][1]=mesh->node.val[j][1];
            } // end for i

            s=areaT(T[0][0], T[1][0], T[2][0], T[0][1], T[1][1], T[2][1]);
            localb(T,btmp,pt->num_qp_rhs);
            tmp=0;
            
            // mesh->elem.col == mesh_aux->elem2edge.col
            for (k1=0;k1<mesh->elem.col;++k1) {

                edge_c = mesh_aux->elem2edge.val[k][k1];
                i=mesh->elem.val[k][k1];

                b->val[i]+=btmp[tmp++];
            } // end for k1
        } // end for k
    }
    else if (mat) {
        for (k=0;k<mesh->elem.row;++k) { 

            for (i=0;i<mesh->elem.col;++i) {
                j=mesh->elem.val[k][i];
                T[i][0]=mesh->node.val[j][0];
                T[i][1]=mesh->node.val[j][1];
            } // end for i

            s=areaT(T[0][0], T[1][0], T[2][0], T[0][1], T[1][1], T[2][1]);
            
            for (k1=0;k1<mesh->elem.col;++k1) {

                edge_c = mesh_aux->elem2edge.val[k][k1];
                i=mesh->elem.val[k][k1];

                // for the diag entry
                gradBasisP1(T, s, k1, phi1);
                tmp_a = 2*s*(phi1[0]*phi1[0]+phi1[1]*phi1[1]*epsilon);
                for (i1=0;i1<pt->num_qp_mat;i1++) {
                    A->val[A->IA[i]] += gauss[i1][2]*tmp_a;
                } // end for i1

                // for the off-diag entry
                i = edge2idx_g1[edge_c];
                j = edge2idx_g2[edge_c];
                n1 = (k1+1)%3;//local node index in the elem
                n2 = (k1+2)%3;
                gradBasisP1(T, s, n1, phi1);
                gradBasisP1(T, s, n2, phi2);
                tmp_a = 2*s*(phi1[0]*phi2[0]+phi1[1]*phi2[1]*epsilon);
                for (i1=0;i1<pt->num_qp_mat;i1++) {
                    A->val[i] += gauss[i1][2]*tmp_a;
                    A->val[j] += gauss[i1][2]*tmp_a; 
                } // end for i1
            } // end for k1
        } // end for k
    }
    else {
    printf(" ###You are not supposed to see this message ...\n");
    }

    fasp_mem_free(edge2idx_g1);
    fasp_mem_free(edge2idx_g2);
}

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int setup_poisson (dCSRmat *A, dvector *b, Mesh *mesh, Mesh_aux *mesh_aux, 
 *                        FEM_param *pt, dvector *ptr_uh, ivector *dof)
 *
 * \brief Setup P1 FEM for the Poisson's equation
 *
 * \param *A                     Pointer to stiffness matrix
 * \param *b                     Pointer to right hand side
 * \param *mesh                  Pointer to mesh info 
 * \param *mesh_aux              Pointer to auxiliary mesh info 
 * \param *pt                    Pointer to parameter 
 * \param *ptr_uh                 Discrete solution
 * \param *ptr_dof                 DOF index
 *
 * \return                       SUCCESS if succeed
 *
 * \author Xuehai Huang, Kai Yang, Chensong Zhang and Feiteng Huang
 * \date   08/10/2010
 * 
 * Modified by Feiteng Huang on 02/21/2012: restructure the code
 * Modified by Feiteng Huang on 04/01/2012, output node, elem, uh, and dof for l2 error
 * Modified by Feiteng Huang on 04/09/2012, restructure the code
 */
int setup_poisson (dCSRmat *A, 
                   dvector *b, 
                   Mesh *mesh, 
                   Mesh_aux *mesh_aux, 
                   FEM_param *pt, 
                   dvector *ptr_uh, 
                   ivector *dof)
{
    dCSRmat Stiff;
    dvector rhs;
    ivector dirichlet,nondirichlet,index;
    int i,j,k;
    int dirichlet_count = 0;

    // assemble A and b
    assemble_stiffmat(&Stiff, &rhs, mesh, mesh_aux, pt);    
    
    // get information to deal with Dirichlet boundary condition
    for (i=0;i<mesh->node_bd.row;++i) {
        if (mesh->node_bd.val[i] == DIRICHLET) dirichlet_count++;
    }
    
    dirichlet = fasp_ivec_create(dirichlet_count); 
    nondirichlet = fasp_ivec_create(mesh->node.row-dirichlet_count); 
    index = fasp_ivec_create(mesh->node.row);
    
    j = k = 0;
    for (i=0;i<mesh->node_bd.row;++i) {
        if(mesh->node_bd.val[i]==DIRICHLET) { //  Dirichlet boundary node
            dirichlet.val[k]=i;
            index.val[i]=k;
            ++k;
        }
        else { // degree of freedom
            nondirichlet.val[j]=i;
            index.val[i]=j;
            ++j;
        }
    }
    
    // set initial boundary value
    dvector uh = fasp_dvec_create(Stiff.row);
    double p[DIM];
    for (i=0;i<uh.row;++i) {
        if(mesh->node_bd.val[i]==DIRICHLET) { // the node is on the boundary
            for (j=0;j<DIM;++j)
                p[j] = mesh->node.val[i][j];
            uh.val[i]=u(p);
        }
    }
    
    extractNondirichletMatrix(&Stiff, &rhs, A, b, &(mesh->node_bd), 
                              &dirichlet, &nondirichlet, &index, &uh);
        
    // output info for l2 error
    ivec_output( dof, &nondirichlet );
    dvec_output( ptr_uh, &uh );
    
    // clean up memory
    fasp_dcsr_free(&Stiff);
    fasp_dvec_free(&rhs);
    fasp_ivec_free(&dirichlet);
    fasp_ivec_free(&index);
    
    return SUCCESS;
}

/** 
 * \fn double get_l2_error_poisson (ddenmat *node,
 *                                  iCSRmat *elem,
 *                                  dvector *uh,
 *                                  int num_qp)
 *
 * \brief get l2 error of fem.
 *
 * \param *node                 Node info
 * \param *elem                 Elem info
 * \param *ptr_uh               Discrete solution
 * \param num_qp                Number of quad points
 *
 * \author Feiteng Huang
 * \date   03/30/2012
 */
double get_l2_error_poisson (ddenmat *node,
                             idenmat *elem,
                             dvector *uh,
                             int num_qp)
{
    double l2error = 0.0;
    double T[3][2];
    double gauss[MAX_QUAD][DIM+1]; 
    double s, uh_local[3], l2, a, p[DIM], uh_p;
    
    int i,j,k;
    
    fasp_gauss2d(num_qp, 2, gauss); // Gauss intergation initial
    
    for (k=0;k<elem->row;++k) { 

        for (i=0;i<elem->col;++i) {
            j=elem->val[k][i];
            T[i][0]=node->val[j][0];
            T[i][1]=node->val[j][1];
            uh_local[i] = uh->val[j];
        } // end for i

        s=2.0*areaT(T[0][0], T[1][0], T[2][0], T[0][1], T[1][1], T[2][1]);

        for (i=0;i<num_qp;++i) {
            l2 = 1 - gauss[i][0] - gauss[i][1];
            for (j=0;j<DIM;++j)
                p[j]=T[0][j]*gauss[i][0]+T[1][j]*gauss[i][1]+T[2][j]*l2;
            a=u(p);
            
            // Improve readability, revise this line --Chensong
        // fixed? --Feiteng
        uh_p = uh_local[0]*gauss[i][0] + uh_local[1]*gauss[i][1] + uh_local[2]*l2;
            l2error+=s*gauss[i][2]*((a - uh_p)*(a - uh_p));
        } // end for i

    } // end for k
    
    l2error = sqrt(l2error);
    
    return l2error;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
