/*! \file refinement2D.inl
 *  \brief Uniform refinement for 2D triangular meshes
 */

/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void getEdgeInfo (iCSRmat *elements, iCSRmat *elementsTran, idenmat *edges, 
 *                              iCSRmat *edgesTran, int nNode)
 *
 * \brief Construct edge information from elements
 *
 * \param *elements      pointer to the structure of the triangulation
 * \param *elementsTran  pointer to the transpose of *elements
 * \param *edges         the first two columns store the two vertice corresponding to the edge; 
 *				         the 3rd and 4th columns store the triangles which the edge belongs to;
 *				         if the edge is a boundary, the 4th column will stores -1;
 *				         the first column is in ascend order.
 * \param *edgesTran     the relation between nodes and edges: JA = edge index, A = another vertex
 *
 * \note ALgorithm: node-->element-->edge
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 * \note   Modified by Kai Yang 03/10/2011: optimize code a little
 * \note   Modify by Feiteng Huang 02/23/2012: pre-alloc memory out of loop
 */
static void getEdgeInfo (iCSRmat *elements, 
                         iCSRmat *elementsTran, 
                         idenmat *edges, 
                         iCSRmat *edgesTran, 
                         int nNode)
{
    int i,j,k,l;
    int element, len;	
    int iSTART=0, point1=0, point2=0;
	
    edges->row=0; edges->col=4; edges->val=NULL;
    
    edgesTran->row=elementsTran->row;
	
    //pre-alloc memory, num of edge < 3*num of node
    edgesTran->IA=(int*)fasp_mem_calloc(edgesTran->row+1, sizeof(int));
    edges->val=(int**)fasp_mem_calloc(3*nNode, sizeof(int *));
    edges->val[0]=(int*)fasp_mem_calloc(4*3*nNode, sizeof(int));
    int *tmp=edges->val[0];
	
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
					
                    edges->val[edges->row]=&tmp[edges->row*4];
					
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
                    edges->val[edges->row]=&tmp[edges->row*4];
                    
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
		
    // realloc edges, we don't really have 3*numOfNode edges
    edges->val=(int**)fasp_mem_realloc(edges->val, sizeof(int *)*(edges->row));
    edges->val[0]=(int*)fasp_mem_realloc(edges->val[0], sizeof(int *)*(edges->row)*(edges->col));
		tmp = edges->val[0];
		for (i=0;i<edges->row;++i)
		{
			edges->val[i]=&tmp[i*edges->col];// re point to the val
		}
	
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
 * \fn static void refine (ddenmat *nodes, iCSRmat *Celements, idenmat *Cedges, 
 *                         iCSRmat *CelementsTran, iCSRmat *Felements, idenmat *Fedges, 
 *                         iCSRmat *FelementsTran, iCSRmat *FedgesTran, ivector *isInNode, 
 *                         ivector *nodeCEdge)
 *
 * \brief Generate fine grid using bisection or regular refinement
 *
 * \param *nodes          pointer to the nodes location of the triangulation
 * \param *Celements      pointer to the structure of the triangulation on coasre grid
 * \param *Cedges         pointer to the edge information of the triangulation on coasre grid
 * \param *CelementsTran  pointer to the transpose of *Celements
 * \param *Felements      pointer to the structure of the triangulation on fine grid
 * \param *Fedges         pointer to the edge information of the triangulation on fine grid
 * \param *FelementsTran  pointer to the transpose of *Felements
 * \param *FedgesTran     the relation between nodes and edges. JA = edge index, A = another vertex
 * \param *isInNode       if the node is interior node, it will be 0
 *                        if the node is on the boundary, it will be -1
 * \param *nodeCEdge      record the index of coarse edge which the node belong to
 *                        if the node is located in the coarset grid, it will be set -1
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 * \note   Modified by Kai Yang 03/10/2011: optimize code a little
 * \note   Modified by Feiteng Huang 02/23/2012: alloc memory out of loop
 */
static void refine (ddenmat *nodes, 
                    iCSRmat *Celements, 
                    idenmat *Cedges, 
                    iCSRmat *CelementsTran, 
                    iCSRmat *Felements,
                    idenmat *Fedges, 
                    iCSRmat *FelementsTran, 
                    iCSRmat *FedgesTran, 
                    ivector *isInNode, 
                    ivector *nodeCEdge)
{
    const int NumCNodes=nodes->row;
    const int Crow=Celements->row;
    const int NUM_NODE=3;
    
    int i,j,k,l,rowstart;
    int point1, point2, location;
	
    // generate fine gird's nodes information: store the 3 middle point of each tirangle
    int *midElements=(int*)fasp_mem_calloc(Crow*NUM_NODE, sizeof(int));

    nodes->row=NumCNodes+Cedges->row;
		// realloc the mem we need
    nodes->val=(double**)fasp_mem_realloc(nodes->val, sizeof(double *)*(nodes->row));
		nodes->val[0]=(double*)fasp_mem_realloc(nodes->val[0], sizeof(double)*(nodes->row)*(nodes->col));
		double *tmp=nodes->val[0];

    for (i=0;i<nodes->row;++i) {
        nodes->val[i]=&tmp[i*nodes->col];// re point to the val
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
#if 0 // bisection refinement
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
#else // regular refinement
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

    getEdgeInfo(Felements, FelementsTran, Fedges, FedgesTran, nodes->row);
	
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
