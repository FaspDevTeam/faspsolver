/*! \file poisson_fem.c
 *  \brief Setup P1 FEM for the Poisson's equation
 */

#include "poisson_fem.h"
#include "basis.inl"
#include "refinement2D.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/** 
 * \fn static void localb (double (*nodes)[2],double *b, int num_qp)
 *
 * \brief Form local right-hand side b from triangle nodes
 *
 * \param (*nodes)[2]  the vertice of the triangule
 * \param *b           local right-hand side
 * \param num_qp       number of quad point 
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 * \note   Modify by Feiteng Huang 02/23/2012: User specified number of quadrature points
 */
static void localb (double (*nodes)[2],
                    double *b,
                    int num_qp)
{
    const double s=2.0*areaT(nodes[0][0],nodes[1][0],nodes[2][0],
                             nodes[0][1],nodes[1][1],nodes[2][1]);
    double x,y,a;
    double gauss[num_qp][3];
    int i;
    
    fasp_init_Gauss(num_qp, 2, gauss); // gauss intergation initial	
	
    for (i=0;i<3;++i) b[i]=0;
	
    for (i=0;i<num_qp;++i) {
        x=nodes[0][0]*gauss[i][0]+nodes[1][0]*gauss[i][1]+nodes[2][0]*(1-gauss[i][0]-gauss[i][1]);
        y=nodes[0][1]*gauss[i][0]+nodes[1][1]*gauss[i][1]+nodes[2][1]*(1-gauss[i][0]-gauss[i][1]);
        a=f(x,y);
		
        b[0]+=a*gauss[i][2]*gauss[i][0];
        b[1]+=a*gauss[i][2]*gauss[i][1];
        b[2]+=a*gauss[i][2]*(1-gauss[i][0]-gauss[i][1]);		
    }
	
    b[0]*=s; b[1]*=s; b[2]*=s;
}

/**
 * \fn static void assemble_stiffmat ( iCSRmat *elements, 
 *                                     iCSRmat *elementsTran, 
 *                                     ddenmat *nodes, 
 *                                     dCSRmat *A, 
 *                                     iCSRmat *edgesTran, 
 *                                     dvector *b, 
 *                                     const char *assemble_option,
 *                                     int num_qp_rhs,
 *                                     int num_qp_mat )
 *
 * \brief Assemble stiffness matrix and right-hand side
 *
 * \param *elements                    pointer to the structure of the triangulation
 * \param *elementsTran                pointer to the transpose of *elements
 * \param *nodes                       pointer to the nodes location of the triangulation
 * \param *A                           pointer to stiffness matrix
 * \param *b                           pointer to right hand side
 * \param *assemble_option             pointer to assemble option
 * \param *num_qp_rhs                  number of quad point for rhs 
 * \param *num_qp_mat                  number of quad point for mat 
 *
 *
 * \note This subroutine follows Ludmil' note on CSR.
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 * \note   Modified by Feiteng Huang on 23/02/2012: restructure assmebling
 */
static void assemble_stiffmat (iCSRmat *elements, 
                               iCSRmat *elementsTran, 
                               ddenmat *nodes, 
                               dCSRmat *A, 
                               iCSRmat *edgesTran, 
                               dvector *b, 
                               const char *assemble_option,
                               int num_qp_rhs,
                               int num_qp_mat)
{
    /* yangkai: use the information of edgesTran to get the CSR structure of A directly, 
     * because, for triangle elements, A_ij is nonzero if and only if nodes i and j are 
     * neighbours in the mesh. 
     */
    int mat = 0, rhs = 0;
    
    if ( strcmp(assemble_option,"a&b") == 0 ) {
        mat = 1;
        rhs = 1;
    }
    
    if ( strcmp(assemble_option,"a") == 0 ) mat = 1;
    
    if ( strcmp(assemble_option,"b") == 0 ) rhs = 1;
    
    const int nodenum = elementsTran->row;
    const int Arow=nodenum;
    const int JAlen= edgesTran->IA[edgesTran->row]-edgesTran->IA[0]+edgesTran->row+1;
    const double epsilon=1;
	
    double T[3][2],phi1[2],phi2[2];
    double gauss[49][3]; // Need to make max number of Gauss a const!!! --Chensong
    double s;
	
    int i,j,k,count;
    int k1,k2,i1,j1;
    double btemp[3];
	
    fasp_init_Gauss(num_qp_mat, 2, gauss); // Gauss intergation initial
	
    A->row = A->col = nodenum;
	
    A->IA=(int*)fasp_mem_calloc(A->row+1, sizeof(int));
    total_alloc_mem += (A->row+1)*sizeof(int);
	
    A->JA=(int*)fasp_mem_calloc(JAlen, sizeof(int));
    total_alloc_mem += JAlen*sizeof(int);
	
    // count is for counting the JA number in A since A-JA and edgesTran-val are different in number 	
    for (count=i=0;i<Arow;++i) {
        A->IA[i]=edgesTran->IA[i]+i;
        A->JA[count]=i; count++;
        for (j=edgesTran->IA[i];j<edgesTran->IA[i+1];++j) {
            A->JA[count]=edgesTran->val[j]; count++;
        }
    }
    A->IA[i]=edgesTran->IA[i]+i;
    A->JA[count]=i; count++;
    
    // Loop element by element and compute the actual entries storing them in A
    A->nnz = A->IA[Arow]-A->IA[0]; // number of nonzeros	
	
    A->val=(double*)fasp_mem_calloc(A->nnz, sizeof(double));
    total_alloc_mem += A->nnz*sizeof(int);
	
    for (k=0;k<elements->row;++k) { 
		
        for (i=elements->IA[k];i<elements->IA[k+1];++i) {
            j=elements->JA[i];
            T[i-elements->IA[k]][0]=nodes->val[j][0];
            T[i-elements->IA[k]][1]=nodes->val[j][1];
        } // end for i
        s=areaT(T[0][0], T[1][0], T[2][0], T[0][1], T[1][1], T[2][1]);
        
        if (rhs) localb(T,btemp,num_qp_rhs);
        
        int tmp=0;
        for (k1=elements->IA[k];k1<elements->IA[k+1];k1++)
        {
            i=elements->JA[k1];
            if (rhs) b->val[i]+=btemp[tmp++];
            if (mat) {
                for (k2=elements->IA[k];k2<elements->IA[k+1];k2++) {
                    j=elements->JA[k2];
                    for (j1=A->IA[i];j1<A->IA[i+1];j1++) {
                        if (A->JA[j1]==j) {
                            for (i1=0;i1<num_qp_mat;i1++) {
                                basisP1(T, s, k1-elements->IA[k], phi1);
                                basisP1(T, s, k2-elements->IA[k], phi2);
                                A->val[j1]+=2*s*gauss[i1][2]*(phi1[0]*phi2[0]+phi1[1]*phi2[1]*epsilon);
                            } // end for i1
                            break;
                        } // end if
                    } // end for j1 
                } // end for k2
            }
        } // end for k1
	} // end for k
	
}

/**
 * \fn static void extractNondirichletMatrix (dCSRmat *A, 
 *                                            dvector *b, 
 *                                            dCSRmat *A11, 
 *                                            dvector *b1, 
 *                                            ivector *isInNode, 
 *                                            ivector *dirichlet, 
 *                                            ivector *nondirichlet, 
 *                                            ivector *index, 
 *                                            dvector *uh)
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
 * \author Xuehai Huang
 * \date   03/29/2009
 */
static void extractNondirichletMatrix (dCSRmat *A, 
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
	
    // add the information in dirichlet and nondirichlet into isInNode
    for (i=0;i<dirichlet->row;++i)
    {
        if (isInNode->val[dirichlet->val[i]]==0) {
            printf("### ERROR: Wrong boundary condition!\n"); exit(RUN_FAIL);
        }
        else isInNode->val[dirichlet->val[i]]=1;
    }
	
    // form A11->IA
    for (i=0;i<A11->row;++i)
    {
        l=nondirichlet->val[i];
        for (k=A->IA[l];k<A->IA[l+1];++k)
        {
            j=A->JA[k];
            if (isInNode->val[j]==0) A11->IA[i+1]++;
        }
    }
	
    for (i=0;i<A11->row;++i) A11->IA[i+1]+=A11->IA[i];
	
    // form A11->JA
    A11->JA=(int*)fasp_mem_calloc(A11->IA[A11->row]+1, sizeof(int));
	
    count=0;
    for (i=0;i<A11->row;++i)
    {
        l=nondirichlet->val[i];
        for (k=A->IA[l];k<A->IA[l+1];++k)
        {
            j=A->JA[k];
            if (isInNode->val[j]==0)
            {
                A11->JA[count]=index->val[j];
                ++count;
            }
        }
    }
	
    // form A11->val
    A11->val=(double*)fasp_mem_calloc(A11->IA[A11->row]+1, sizeof(double));
	
    for (i1=0;i1<A11->row;++i1)
    {
        i=nondirichlet->val[i1];
        for (k=A11->IA[i1];k<A11->IA[i1+1];++k)
        {
            j1=A11->JA[k];
            j=nondirichlet->val[j1];
			
            for (l=A->IA[i];l<A->IA[i+1];l++)
            {
                if (A->JA[l]==j)
                {
                    A11->val[k]=A->val[l];
                    break;
                }
            }
        }
    }
    A11->nnz=A11->IA[A11->row]-A11->IA[0];
	
    // yangkai: In the following part, I have changed the loop order of previous code and reduced time. 	
    for (i=0;i<A->row;++i) {
        for (k=A->IA[i];k<A->IA[i+1];++k) {
            if (isInNode->val[A->JA[k]]==1) b->val[i]-=A->val[k]*uh->val[A->JA[k]]; // dirichlet point
        }
    }
	
    for (i1=0;i1<b1->row;++i1) {
        i=nondirichlet->val[i1];
        b1->val[i1]=b->val[i];
    }
	
}


/**
 * \fn static iCSRmat getCSRfromIntT (idenmat *T, int nNode)
 *
 * \brief Get CSR structure of trianglution T
 *
 * \param *T     pointer to the idenmat
 * \param nNode  the number of nodes
 *
 * \return       trianglution T in CSR format
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 */
static iCSRmat getCSRfromIntT (idenmat *T, 
                               int nNode)
{
    iCSRmat element;
    int i,j;
	
    element.val=NULL;
    element.JA=NULL;
    element.row=T->row;
    element.col=nNode;
    element.nnz=T->col*T->row;
	
    element.IA=(int*)fasp_mem_calloc(element.row+1, sizeof(int));
	
    for (i=0;i<element.row;++i) element.IA[i+1]=element.IA[i]+T->col;
	
    element.JA=(int*)fasp_mem_realloc(element.JA, (element.IA[element.row])*sizeof(int));
	
    for (i=0;i<T->row;++i) {
        for (j=0;j<T->col;++j) {
            element.JA[N2C(element.IA[i])+j]=T->val[i][j];
        }
    }
    
    return element;
}

/**
 * \fn static int getMeshInfo (ddenmat *nodes, iCSRmat *elements, idenmat *edges, 
 *                             iCSRmat *elementsTran, iCSRmat *edgesTran, ivector *isInNode, 
 *                             ivector *nodeCEdge, const char *filename)
 *
 * \brief generate the mesh information and store it into point, element, edge respectively
 *
 * \param *nodes         the first column stores the x coordinate of points
 *                       the second column stores the y coordinate of points
 * \param *elements      using CSR format to store 3 nodes corresponding to the element
 * \param *edges         the first two columns store the two vertice
 *                       the third column stores the affiliated element
 * \param *elementsTran  pointer to the transpose of *elements
 * \param *edgesTran     the relation between nodes and edges. JA = edge index, A = another vertex
 * \param *isInNode      if the node is interior node, it will be 0
 *                       if the node is on the boundary, it will be -1
 * \param *nodeCEdge     record the index of coarse edge which the node belong to
 *                       if the node is located in the coarset grid, it will be set -1
 * \param *filename      input mesh file
 *
 * \return               1 if succeed 0 if fail
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 * \note   Modify by Feiteng Huang 02/23/2012: alloc memory out of loop.
 */
static int getMeshInfo (ddenmat *nodes, 
                        iCSRmat *elements, 
                        idenmat *edges, 
                        iCSRmat *elementsTran, 
                        iCSRmat *edgesTran, 
                        ivector *isInNode, 
                        ivector *nodeCEdge, 
                        const char *filename)
{
    int NumNode, Ncoor, NumElem, Nnode;
    int i,j;
    
    FILE *inputFile=fopen(filename, "r");
    if (inputFile==NULL)
    {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(-2);
    }
	
    // get the nodes' coordinates
    fscanf(inputFile, "%d %d", &NumNode, &Ncoor);
    nodes->row=NumNode;
    nodes->col=Ncoor;
    // alloc mem we need,
    nodes->val=(double **)fasp_mem_calloc(nodes->row, sizeof(double*)); 
    nodes->val[0]=(double *)fasp_mem_calloc(nodes->col*nodes->row, sizeof(double)); 
    double *tmp = nodes->val[0];
    
    for (i=0;i<nodes->row;++i)// re point to the val
    {
        nodes->val[i]=&tmp[nodes->col*i];		
        for (j=0;j<nodes->col;++j) fscanf(inputFile, "%lf", &nodes->val[i][j]);
    }
	
    // get triangular grid
    fscanf(inputFile, "%d %d", &NumElem, &Nnode);
	
    elements->val=NULL;
    elements->JA=NULL;
    elements->row=NumElem;
    elements->col=NumNode;
    elements->nnz=Nnode*NumElem;
	
    elements->IA=(int*)fasp_mem_calloc(elements->row+1, sizeof(int));
	
    for (i=0;i<elements->row;++i) elements->IA[i+1]=elements->IA[i]+Nnode;
	
    elements->JA=(int*)fasp_mem_realloc(elements->JA, (elements->IA[elements->row])*sizeof(int));
	
    for (i=0;i<NumElem;++i) {
        for (j=0;j<Nnode;++j) {
            fscanf(inputFile, "%d", &(elements->JA[N2C(elements->IA[i])+j]));
            (elements->JA[N2C(elements->IA[i])+j])--; 
            // the data from matlab start with 1 while from c start with 0
        }
    }	
	
    fclose(inputFile);
    
    fasp_icsr_trans(elements,elementsTran);
	
    getEdgeInfo(elements, elementsTran, edges, edgesTran, NumNode);
	
    // generate isInNode
    isInNode->row=nodes->row;
    isInNode->val=(int*)fasp_mem_calloc(isInNode->row, sizeof(int));
	
    for (i=0;i<edges->row;++i)
    {
        if (edges->val[i][3]==-1) // case the edge is on boundary
        {
            isInNode->val[edges->val[i][0]]=-1;
            isInNode->val[edges->val[i][1]]=-1;
        }
    }
	
    // get nodeCEdge
    nodeCEdge->row=nodes->row;
    nodeCEdge->val=(int*)fasp_mem_calloc(nodeCEdge->row, sizeof(int));
	
    for (i=0;i<nodeCEdge->row;++i) nodeCEdge->val[i]=-1;
	
    return SUCCESS;
}

/**
 * \fn static void write_mesh_info (ddenmat *nodes, iCSRmat *elements, const char *filename)
 *
 * \brief Write the mesh information
 *
 * \param *nodes         the first column stores the x coordinate of points, 
 *                       the second column stores the y coordinate of points
 * \param *elements      using CSR format to store 3 nodes corresponding to the element
 * \param *filename      filename of output file
 *
 * \author Feiteng Huang
 * \date   02/22/2012
 */
static void write_mesh_info (ddenmat *nodes, 
                             iCSRmat *elements, 
                             const char *filename)
{
    int NumNode=nodes->row;
    int Ncoor=nodes->col;
    int NumElem=elements->row;
    int Nnode=elements->IA[1]-elements->IA[0];
    int i,j;
    
    FILE *outputFile=fopen(filename, "w");
    if (outputFile==NULL)
    {
        printf("### ERROR: Opening file %s fails!\n", filename);
        exit(-2);
    }
	
    fprintf(outputFile, "%d %d\n", NumNode, Ncoor);
    
    for (i=0;i<NumNode;++i)
    {
        for (j=0;j<Ncoor;++j) {
            fprintf(outputFile, "%lf ", nodes->val[i][j]);
        }
        fprintf(outputFile, "\n");
    }
	
    fprintf(outputFile, "%d %d\n", NumElem, Nnode);
    
    for (i=0;i<NumElem;++i) {
        for (j=0;j<Nnode;++j) {
            fprintf(outputFile, "%d ", elements->JA[N2C(elements->IA[i])+j]+1);
        }
        fprintf(outputFile, "\n");
    }	
	
    fclose(outputFile);
}

/**
 * \fn int setup_poisson (dCSRmat *ptr_A, 
 *                        dvector *ptr_b, 
 *                        int levelNum, 
 *                        const char *meshIn,
 *                        const char *meshOut, 
 *                        int mo, 
 *                        const char *assemble_option, 
 *                        int num_qp_rhs, 
 *                        int num_qp_mat)
 *
 * \brief Setup P1 FEM for the Poisson's equation
 *
 * \param *ptr_A                 pointer to stiffness matrix
 * \param *ptr_b                 pointer to right hand side
 * \param levelNum               total level number of grid
 * \param *meshIn                input mesh file 
 * \param *meshOut               output mesh file 
 * \param *mo                    output mesh info option 
 * \param *assemble_option       assemble option
 * \param num_qp_rhs             number of quad point for rhs 
 * \param num_qp_mat             number of quad point for mat 
 *
 * \return                       SUCCESS if succeed
 *
 * \author Xuehai Huang, Kai Yang, Chensong Zhang
 * \date   08/10/2010
 * \note   Modified by Feiteng Huang on 02/21/2012: restructure the code
 */
int setup_poisson (dCSRmat *ptr_A, 
                   dvector *ptr_b, 
                   int levelNum, 
                   const char *meshIn, 
                   const char *meshOut, 
                   int mo, 
                   const char *assemble_option,
                   int num_qp_rhs,
                   int num_qp_mat)
{
    clock_t all_start=clock();
    
    ddenmat nodes; 
    // the first column = the x coordinate of points, the second column = the y coordinate of points
    
    iCSRmat elements[MAX_REFINE_LVL]; 
    // triangulation: store 3 points corresponding to the element in each row
    
    idenmat edges[MAX_REFINE_LVL]; 
    // the first two columns store the two vertice, the third column stores the affiliated element
    
    iCSRmat elementsTran[MAX_REFINE_LVL]; // the transpose of elements. 
    
    iCSRmat edgesTran[MAX_REFINE_LVL]; 
    // the tranpose of elements, used to get restriction operator. 
    // The relation between nodes and edges. JA stores edge index, A stores another vertex
    
    ivector nodeCEdge; 
    // record the index of coarse edge which the node belong to
    // if the node is located in the coarset grid, it will be set -1
    
    ivector isInNode; 
    // 0: if the node is interior node; -1: if the node is on the boundary.
    
    ivector dirichlet,nondirichlet,index;
    int i,j,k,l;
	
    if (levelNum>=MAX_REFINE_LVL) {
        printf("### WARNING: Refinement level %d excceds max refinement level = %d\n", 
               levelNum, MAX_REFINE_LVL);
        levelNum = MAX_REFINE_LVL-1;
    }
	
    // get mesh info from a disk file
    clock_t meshRead_start=clock();	
    
    if ( getMeshInfo(&nodes, &elements[0], &edges[0], &elementsTran[0], 
                     &edgesTran[0], &isInNode, &nodeCEdge, meshIn) != SUCCESS ) 
    {
        printf("### ERROR: Read mesh information failed!\n");
        exit(RUN_FAIL);
    }
    
    clock_t meshRead_end=clock();	
    double meshRead_time = (double)(meshRead_end - meshRead_start)/(double)(CLOCKS_PER_SEC);
    printf("Reading mesh costs ...%8.4f seconds.\n", meshRead_time); 
	
    // refinements
    clock_t refine_start=clock();	
	
    for (i=0;i<levelNum-1;++i) {	
        
        fasp_ivec_free(&isInNode);
        refine(&nodes, &elements[i], &edges[i], &elementsTran[i], &elements[i+1], 
               &edges[i+1], &elementsTran[i+1], &edgesTran[i+1], &isInNode, &nodeCEdge);
        
        if (mo)
        {
            char filename[120];
            sprintf(filename, "%s%d.dat", meshOut, i+1);
            write_mesh_info(&nodes, &elements[i+1], filename);
        }
        
    }
	
    clock_t refine_end=clock();
	
    double refine_time = (double)(refine_end - refine_start)/(double)(CLOCKS_PER_SEC);
    printf("Refinement costs .....%8.4f seconds.\n", refine_time); 
	
    // get information to deal with Dirichlet boundary condition
    int dirichlet_count[levelNum];
    dirichlet_count[0]=0;
    for (j=0;j<elementsTran[0].row;++j)
    {
        if (isInNode.val[j]==-1) dirichlet_count[0]++;
    }
    for (i=1;i<levelNum;++i)
    {
        dirichlet_count[i]=dirichlet_count[i-1];
        for (j=elementsTran[i-1].row;j<elementsTran[i].row;++j)
        {
            if(isInNode.val[j]==-1) dirichlet_count[i]++;
        }
    }
	
    dirichlet = fasp_ivec_create(dirichlet_count[levelNum-1]); 
    nondirichlet = fasp_ivec_create(elementsTran[levelNum-1].row-dirichlet_count[levelNum-1]); 
    index = fasp_ivec_create(elementsTran[levelNum-1].row);
	
    j = k = 0;
    for (i=0;i<elementsTran[levelNum-1].row;++i)
    {
        if(isInNode.val[i]==-1) //  Dirichlet boundary node
        {
            dirichlet.val[k]=i;
            index.val[i]=k;
            ++k;
        }
        else // free variable
        {
            nondirichlet.val[j]=i;
            index.val[i]=j;
            ++j;
        }
    }
	
    dirichlet.row=k; nondirichlet.row=j;
	
    // assemble stiffness matrix A
    dCSRmat A;
    dvector b = fasp_dvec_create(elementsTran[levelNum-1].row);
	
    clock_t assemblestif_start=clock();
	
    assemble_stiffmat(&elements[levelNum-1], &elementsTran[levelNum-1], &nodes, &A, 
                      &edgesTran[levelNum-1], &b, assemble_option, num_qp_rhs, num_qp_mat);
	
    clock_t assemblestif_end=clock();
	
    double assemblestif_time = (double)(assemblestif_end - assemblestif_start)/(double)(CLOCKS_PER_SEC);
	
    printf("Assemble %3s costs ...%8.4f seconds.\n", assemble_option, assemblestif_time); 
    
    // set initial boundary value
    dvector uh = fasp_dvec_create(A.row);
	
    for (i=0;i<uh.row;++i)
    {
        if(isInNode.val[i]==-1) // the node is on the boundary
        { 
            uh.val[i]=u(nodes.val[i][0], nodes.val[i][1]);
        }
    }
	
    extractNondirichletMatrix(&A, &b, ptr_A, ptr_b, &isInNode, 
                              &dirichlet, &nondirichlet, &index, &uh);
	
    // clean up memory
    fasp_mem_free(nodes.val[0]);
    fasp_mem_free(nodes.val);
	
    for (l=0;l<levelNum;l++)
    {
        fasp_icsr_free(&elements[l]);
        fasp_icsr_free(&elementsTran[l]);
        fasp_icsr_free(&edgesTran[l]);
        fasp_mem_free(edges[l].val[0]);
        fasp_mem_free(edges[l].val);
    }
	
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&uh);
    fasp_ivec_free(&dirichlet);
    fasp_ivec_free(&nondirichlet);
    fasp_ivec_free(&index);
    fasp_ivec_free(&isInNode);
    fasp_ivec_free(&nodeCEdge);
	
    clock_t all_end=clock();
	
    double all_time = (double)(all_end - all_start)/(double)(CLOCKS_PER_SEC);
    printf("FEM setup costs ......%8.4f seconds.\n", all_time);
    
    return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
