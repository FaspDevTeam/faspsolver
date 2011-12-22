/*! \file assemble_info.c
 *  \brief Subroutines for assembling purpose
 */

/**
 * \fn iCSRmat getCSRfromIntT(idenmat *T, int nNode)
 * \brief Get CSR structure of trianglution T.
 *
 * \param *T     pointer to the idenmat
 * \param nNode  the number of nodes
 * \return       trianglution T in CSR format
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
static iCSRmat getCSRfromIntT(idenmat *T, int nNode)
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
	fasp_mem_check(element.JA,"getCSRfromIntT: cannot allocate memory!\n",ERROR_ALLOC_MEM);
	
	for (i=0;i<T->row;++i) {
		for (j=0;j<T->col;++j) {
			element.JA[N2C(element.IA[i])+j]=T->val[i][j];
		}
	}	
	return element;
}

/**
 * \fn static void assemble_stiffmat(iCSRmat *elements, iCSRmat *elementsTran, ddenmat *nodes, dCSRmat *A)
 * \brief assemble stiffness matrix 
 *
 * \param *elements      pointer to the structure of the triangulation
 * \param *elementsTran  pointer to the transpose of *elements
 * \param *nodes         pointer to the nodes location of the triangulation
 * \param *A             pointer to stiffness matrix
 *
 *	 \note this subroutine follows Ludmil' note about CSR.
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
static void assemble_stiffmat(iCSRmat *elements, iCSRmat *elementsTran, ddenmat *nodes, dCSRmat *A, iCSRmat *edgesTran)
{
	/* yangkai: use the information of edgesTran to get the CSR structure of A directly, 
	 * because, for triangle elements, A_ij is nonzero if and only if nodes i and j are 
	 * neighbours in the mesh. 
	 */
	const int nodenum = elementsTran->row;
	const int Arow=nodenum;
	const int JAlen= edgesTran->IA[edgesTran->row]-edgesTran->IA[0]+edgesTran->row+1;
	const int num_qp=1; // the number of numerical intergation points
	const double epsilon=1;
	
	double T[3][2],phi1[2],phi2[2];
	double gauss[49][3]; // Need to make max number of Gauss a const!!! --Chensong
	double s;
	
	int i,j,k,count;
	int k1,k2,i1,j1;
	
	fasp_init_Gauss(num_qp, 2, gauss); // gauss intergation initial
	
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
		s=area(T[0][0], T[1][0], T[2][0], T[0][1], T[1][1], T[2][1]);
		
		for (k1=elements->IA[k];k1<elements->IA[k+1];k1++)
		{
			i=elements->JA[k1];
			for (k2=elements->IA[k];k2<elements->IA[k+1];k2++) {
				j=elements->JA[k2];
				for (j1=A->IA[i];j1<A->IA[i+1];j1++) {
					if (A->JA[j1]==j) {
						for (i1=0;i1<num_qp;i1++) {
							basis(T, s, k1-elements->IA[k], phi1);
							basis(T, s, k2-elements->IA[k], phi2);
							A->val[j1]+=2*s*gauss[i1][2]*(phi1[0]*phi2[0]+phi1[1]*phi2[1]*epsilon);
						} // end for i1
						break;
					} // end if
				} // end for j1 
			} // end for k2
		} // end for k1
		
	} // end for k
	
}

/**
 * \fn void getEdgeInfo(iCSRmat *elements, iCSRmat *elementsTran, idenmat *edges, iCSRmat *edgesTran)
 * \brief get edges information from elements
 *
 * \param *elements      pointer to the structure of the triangulation
 * \param *elementsTran  pointer to the transpose of *elements
 * \param *edges         the first two columns store the two vertice corresponding to the edge; 
 *				                the 3rd and 4th columns store the triangles which the edge belongs to;
 *				                if the edge is a boundary, the 4th column will stores -1;
 *				                the first column is in ascend order.
 * \param *edgesTran     the relation between nodes and edges. JA stores edge index, A stores another vertex
 *
 * \note  ALgorithm: node-->element-->edge
 *
 * TODO: modify the algorithm for finding edges, alloc should not be used in the loop; 
 * first count the edges and then allocate. --Chensong
 *
 * \author Xuehai Huang, Kai Yang
 * \date 03/29/2009
 */
static void getEdgeInfo(iCSRmat *elements, iCSRmat *elementsTran, idenmat *edges, iCSRmat *edgesTran)
{
	int i,j,k,l;
	int element, len;	
	int iSTART=0, point1=0, point2=0;
	
	edges->row=0; edges->col=4; edges->val=NULL;
    
	edgesTran->row=elementsTran->row;
	
	edgesTran->IA=(int*)fasp_mem_calloc(edgesTran->row+1, sizeof(int));
	
	for (i=0;i<elementsTran->row;++i)
	{
		for (j=elementsTran->IA[i];j<elementsTran->IA[i+1];++j)
		{
			element=elementsTran->JA[j]; // achieve the element
			for (k=elements->IA[element];k<elements->IA[element+1];++k)
			{
				if (elements->JA[k]==i)
				{
					len=elements->IA[element+1]-elements->IA[element];
					point1=elements->JA[elements->IA[element]+(k-elements->IA[element]-1+len)%len];
					point2=elements->JA[elements->IA[element]+(k-elements->IA[element]+1)%len];
					break;
				}
			} // k
			
			if (i<point1) // the case i>point1 already exists in edges
			{
				for (l=iSTART;l<edges->row;l++)
				{					
					if ((edges->val[l][0]==i) && (edges->val[l][1]==point1))
					{
						break;
					}
				}
				if (l==edges->row) // the edge(i, point1) is new 
				{
					edges->val=(int**)fasp_mem_realloc(edges->val, sizeof(int *)*(edges->row+1));
					
					edges->val[edges->row]=(int*)fasp_mem_calloc(edges->col, sizeof(int));
					
					edges->row++;
					edges->val[edges->row-1][0]=i;
					edges->val[edges->row-1][1]=point1;
					edges->val[edges->row-1][2]=element;
					edges->val[edges->row-1][3]=-1;
					edgesTran->IA[i+1]++;
					edgesTran->IA[point1+1]++;
				}
				else // the edge(i, point1) already exits
				{
					edges->val[l][3]=element;
				}				
			}
			
			if (i<point2) // the case i>point2 already exists in edges
			{
				for (l=iSTART;l<edges->row;l++)
				{					
					if ((edges->val[l][0]==i) && (edges->val[l][1]==point2))
					{
						break;
					}
				}
				if (l==edges->row) // the edge(i, point2) is new
				{
					edges->val=(int**)fasp_mem_realloc(edges->val, sizeof(int *)*(edges->row+1));
					
					edges->val[edges->row]=(int*)fasp_mem_calloc(edges->col, sizeof(int));
										
					edges->row++;
					edges->val[edges->row-1][0]=i;
					edges->val[edges->row-1][1]=point2;
					edges->val[edges->row-1][2]=element;
					edges->val[edges->row-1][3]=-1;
					edgesTran->IA[i+1]++;
					edgesTran->IA[point2+1]++;
				}
				else // the edge(i, point2) already exits
				{
					edges->val[l][3]=element;
				}				
			}
		} // j
		
		iSTART=edges->row;
	} // i
	
	// generate egdesTran
	edgesTran->col=edges->row;
	for (i=0;i<edgesTran->row;++i)
	{
		edgesTran->IA[i+1]=edgesTran->IA[i+1]+edgesTran->IA[i];
	}
	
	edgesTran->JA=(int*)fasp_mem_calloc(edgesTran->IA[edgesTran->row], sizeof(int));
	
	edgesTran->val=(int*)fasp_mem_calloc(edgesTran->IA[edgesTran->row], sizeof(int));
	
	int *count; // for counting in each row
	count=(int*)fasp_mem_calloc(edgesTran->row, sizeof(int));
	
	for (i=0;i<edges->row;++i)
	{
		point1=edges->val[i][0];
		point2=edges->val[i][1];
		edgesTran->JA[edgesTran->IA[point1]+count[point1]]=i;
		edgesTran->val[edgesTran->IA[point1]+count[point1]]=point2;
		count[point1]++;
		edgesTran->JA[edgesTran->IA[point2]+count[point2]]=i;
		edgesTran->val[edgesTran->IA[point2]+count[point2]]=point1;
		count[point2]++;
	}
	
	fasp_mem_free(count);
}

/**
 * \fn int getCoarseInfo(ddenmat *nodes, iCSRmat *elements, idenmat *edges, iCSRmat *elementsTran, iCSRmat *edgesTran, ivector *isInNode, ivector *nodeCEdge)
 * \brief generate the coarse grid information and store it into point, element, edge respectively
 *
 * \param *nodes         the first column stores the x coordinate of points, the second column stores the y coordinate of points
 * \param *elements      using CSR format to store 3 nodes corresponding to the element
 * \param *edges         the first two columns store the two vertice, the third column stores the affiliated element
 * \param *elementsTran  pointer to the transpose of *elements
 * \param *edgesTran     the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *isInNode      if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *nodeCEdge     record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \return               1 if succeed 0 if fail
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
static int getCoarseInfo(ddenmat *nodes, iCSRmat *elements, idenmat *edges, iCSRmat *elementsTran, iCSRmat *edgesTran, ivector *isInNode, ivector *nodeCEdge)
{
	int NumNode, Ncoor, NumElem, Nnode;
	int i,j;

	// Should change this: get data from an input file!!!
	char *filename="./data/testmesh.dat";
	
	FILE *inputFile=fopen(filename, "r");
	if (inputFile==NULL)
	{
		printf("Opening file %s fails!\n", filename);
		exit(-2);
	}
	
	// get the nodes' coordinates
	fscanf(inputFile, "%d %d", &NumNode, &Ncoor);
	nodes->row=NumNode;
	nodes->col=Ncoor;
	nodes->val=(double **)fasp_mem_calloc(nodes->row, sizeof(double*)); 
		
	for (i=0;i<nodes->row;++i)
	{
		nodes->val[i]=(double*)fasp_mem_calloc(nodes->col, sizeof(double));		
		for (j=0;j<nodes->col;++j) fscanf(inputFile, "%lf", &nodes->val[i][j]);
	}
	
	// get triangular grid
	idenmat T;
	fscanf(inputFile, "%d %d", &NumElem, &Nnode);
	
	T.row=NumElem;
	T.col=Nnode;	
	T.val=(int **)fasp_mem_calloc(T.row, sizeof(int*)); 
	
	for (i=0;i<T.row;++i)
	{
		T.val[i]=(int*)fasp_mem_calloc(T.col, sizeof(int)); 
		
		for (j=0;j<T.col;++j)
		{
			fscanf(inputFile, "%d", &T.val[i][j]);
			T.val[i][j]--; // the data from matlab start with 1 while from c start with 0
		}
	}	
	
	fclose(inputFile);
	
	(*elements)=getCSRfromIntT(&T,nodes->row);
	
	fasp_icsr_trans(elements,elementsTran);
	
	// get edge information
	getEdgeInfo(elements, elementsTran, edges, edgesTran);
	
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
	
	fasp_iden_free(&T);
	
	return 1;
}

/**
 * \fn void refine(ddenmat *nodes, iCSRmat *Celements, idenmat *Cedges, iCSRmat *CelementsTran, iCSRmat *Felements, idenmat *Fedges, iCSRmat *FelementsTran, iCSRmat *FedgesTran, ivector *isInNode, ivector *nodeCEdge)
 * \brief generate fine grid using regular section
 *
 * \param *nodes          pointer to the nodes location of the triangulation
 * \param *Celements      pointer to the structure of the triangulation on coasre grid
 * \param *Cedges         pointer to the edge information of the triangulation on coasre grid
 * \param *CelementsTran  pointer to the transpose of *Celements
 * \param *Felements      pointer to the structure of the triangulation on fine grid
 * \param *Fedges         pointer to the edge information of the triangulation on fine grid
 * \param *FelementsTran  pointer to the transpose of *Felements
 * \param *FedgesTran     the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *isInNode       if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *nodeCEdge      record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 *
 * \author Xuehai Huang, Kai Yang
 * \date 03/29/2009
 */
static void refine(ddenmat *nodes, iCSRmat *Celements, idenmat *Cedges, iCSRmat *CelementsTran, iCSRmat *Felements, \
						idenmat *Fedges, iCSRmat *FelementsTran, iCSRmat *FedgesTran, ivector *isInNode, ivector *nodeCEdge)
{
	const int NumCNodes=nodes->row;
	const int Crow=Celements->row;
	const int NUM_NODE=3;

	int i,j,k,l,rowstart;
	int point1, point2, location;
	
	// generate fine gird's nodes information: store the 3 middle point of each tirangle
	int *midElements=(int*)fasp_mem_calloc(Crow*NUM_NODE, sizeof(int));
	
	nodes->row=NumCNodes+Cedges->row;
	nodes->val=(double**)fasp_mem_realloc(nodes->val, sizeof(double *)*(nodes->row));
	
	for (i=NumCNodes;i<nodes->row;++i) {
		nodes->val[i]=(double*)fasp_mem_calloc(nodes->col, sizeof(double));
	}
	
	nodeCEdge->row=nodes->row;
	nodeCEdge->val=(int*)fasp_mem_realloc(nodeCEdge->val, sizeof(int)*(nodeCEdge->row));
	
	for (i=0;i<Cedges->row;++i)
	{
		point1=Cedges->val[i][0];
		point2=Cedges->val[i][1];
		nodes->val[NumCNodes+i][0]=(nodes->val[point1][0]+nodes->val[point2][0])/2;
		nodes->val[NumCNodes+i][1]=(nodes->val[point1][1]+nodes->val[point2][1])/2;
		nodeCEdge->val[NumCNodes+i]=i;
		
		// get the relation between triangle and edge
		l=Cedges->val[i][2];
		for (j=0;j<(Celements->IA[l+1]-Celements->IA[l]);++j)
		{
			if (Celements->JA[Celements->IA[l]+j]==point1)
			{
				break;
			}
		}
		for (k=0;k<(Celements->IA[l+1]-Celements->IA[l]);++k)
		{
			if (Celements->JA[Celements->IA[l]+k]==point2)
			{
				break;
			}
		}
		location=3-j-k;
		midElements[l*NUM_NODE+location]=i;
		if (Cedges->val[i][3]>-1) // case this edge is not a boundary edge
		{
			l=Cedges->val[i][3];
			for (j=0;j<(Celements->IA[l+1]-Celements->IA[l]);++j)
			{
				if (Celements->JA[Celements->IA[l]+j]==point1)
				{
					break;
				}
			}
			for (k=0;k<(Celements->IA[l+1]-Celements->IA[l]);++k)
			{
				if (Celements->JA[Celements->IA[l]+k]==point2)
				{
					break;
				}
			}
			location=3-j-k;
			midElements[l*NUM_NODE+location]=i;
		}
	}
	
	// generate fine grid trianglues information
	Felements->row=Crow*4;
	Felements->col=nodes->row;	
	Felements->val=NULL;
	Felements->JA=NULL;
	Felements->IA=(int*)fasp_mem_calloc(Felements->row+1, sizeof(int));
	
	for (i=0;i<Felements->row;++i)
	{
		Felements->IA[i+1]=Felements->IA[i]+3;  // triangluation in 2-dimension case
	}
	Felements->nnz = Felements->IA[Felements->row];
	Felements->JA=(int*)fasp_mem_realloc(Felements->JA, (Felements->nnz)*sizeof(int));
	
	for (i=0;i<Crow;++i)
	{
		rowstart=Celements->IA[i];
		
		// Give an option for choosing between bisection and regular refinement!!!	
#if 0 // bisection
		 Felements->JA[i*12] = midElements[i*NUM_NODE+2]+NumCNodes;
		 Felements->JA[i*12+1] = midElements[i*NUM_NODE+0]+NumCNodes;
		 Felements->JA[i*12+2] = Celements->JA[rowstart];
		 Felements->JA[i*12+3] = midElements[i*NUM_NODE+2]+NumCNodes;
		 Felements->JA[i*12+4] = Celements->JA[rowstart+1];
		 Felements->JA[i*12+5] = midElements[i*NUM_NODE]+NumCNodes;
		 Felements->JA[i*12+6] = midElements[i*NUM_NODE+1]+NumCNodes;
		 Felements->JA[i*12+7] = midElements[i*NUM_NODE]+NumCNodes;
		 Felements->JA[i*12+8] = Celements->JA[rowstart+2];
		 Felements->JA[i*12+9] = midElements[i*NUM_NODE+1]+NumCNodes;
		 Felements->JA[i*12+10] = Celements->JA[rowstart];
		 Felements->JA[i*12+11] = midElements[i*NUM_NODE]+NumCNodes;
#else // regular section start
		Felements->JA[i*12] = Celements->JA[rowstart];
		Felements->JA[i*12+1] = midElements[i*NUM_NODE+2]+NumCNodes;
		Felements->JA[i*12+2] = midElements[i*NUM_NODE+1]+NumCNodes;
		Felements->JA[i*12+3] = Celements->JA[rowstart+1];
		Felements->JA[i*12+4] = midElements[i*NUM_NODE]+NumCNodes;
		Felements->JA[i*12+5] = midElements[i*NUM_NODE+2]+NumCNodes;
		Felements->JA[i*12+6] = Celements->JA[rowstart+2];
		Felements->JA[i*12+7] = midElements[i*NUM_NODE+1]+NumCNodes;
		Felements->JA[i*12+8] = midElements[i*NUM_NODE]+NumCNodes;
		Felements->JA[i*12+9] = midElements[i*NUM_NODE]+NumCNodes;
		Felements->JA[i*12+10] = midElements[i*NUM_NODE+1]+NumCNodes;
		Felements->JA[i*12+11] = midElements[i*NUM_NODE+2]+NumCNodes;  
#endif	
  }
	
	fasp_icsr_trans(Felements,FelementsTran);

	getEdgeInfo(Felements, FelementsTran, Fedges, FedgesTran);
	
	isInNode->row=nodes->row;
	isInNode->val=(int*)fasp_mem_calloc(isInNode->row, sizeof(int));
	
	for (i=0;i<Fedges->row;++i)
	{
		if (Fedges->val[i][3]==-1) // in case the edge is on boundary
		{
			isInNode->val[Fedges->val[i][0]]=-1;
			isInNode->val[Fedges->val[i][1]]=-1;
		}
	}
			
	fasp_mem_free(midElements);		
}

/**
 * \fn void extractNondirichletMatrix(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index, dvector *uh)
 * \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition  
 *
 * \param *A            pointer to the stiffness matrix with dirichelt boundary condition(without removed)
 * \param *b            pointer to the right hand side with dirichelt boundary condition(without removed)
 * \param *A11          pointer to the stiffness matrix without dirichelt boundary condition(removed)
 * \param *b1           pointer to the right hand side without dirichelt boundary condition(removed)
 * \param *isInNode     if the node is interior node, it will be 0; if the node is on the boundary, it will be -1---- mod: if it is dirichlet boundary it will be 1
 * \param *dirichlet    pointer to the indicator of the dirichlet boundary
 * \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
 * \param *index        pointer to the transpose of *dirichlet and *nondirichlet
 * \param *uh           pointer to the dirichlet boundary value
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
static void extractNondirichletMatrix(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index, dvector *uh)
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
			printf("Fatal error\n"); exit(RUN_FAIL);
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
