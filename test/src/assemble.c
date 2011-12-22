/*! \file assemble.c
 *  \brief Assembling for stiffness matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

#include "basisP1.inl"
#include "assemble_util.inl"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn int assemble(dCSRmat *ptr_A, dvector *ptr_b, int levelNum)
 * \brief assemble stiffness matrix *ptr_A and righ hand side *ptr_b
 *
 * \param *ptr_A    pointer to stiffness matrix
 * \param *ptr_b    pointer to right hand side
 * \param levelNum  total level number of grid
 * \return          1 if succeed 0 if fail
  *
 * \author Xuehai Huang, Kai Yang, Chensong Zhang
 * \date 08/10/2010
 */
int assemble(dCSRmat *ptr_A, dvector *ptr_b, int levelNum)
{
	ddenmat nodes; // the first column stores the x coordinate of points, the second column stores the y coordinate of points
	iCSRmat elements[MAX_REFINE_LVL]; // triangulation: store 3 points corresponding to the element in each row
	idenmat edges[MAX_REFINE_LVL]; // the first two columns store the two vertice, the third column stores the affiliated element
	iCSRmat elementsTran[MAX_REFINE_LVL]; // the transpose of elements. 
	iCSRmat edgesTran[MAX_REFINE_LVL]; // the tranpose of elements, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
	ivector nodeCEdge; // record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
	ivector isInNode; // 0: if the node is interior node; -1: if the node is on the boundary.
	ivector dirichlet,nondirichlet,index;
	int i,j,k,l;
	
	if (levelNum>=MAX_REFINE_LVL) {
		levelNum = MAX_REFINE_LVL-1;
		printf("Warning: Refinement level %d excceds max refinement level = %d\n", levelNum, MAX_REFINE_LVL);
	}
	
	// get coarse grid from a disk file
	int IsExist=getCoarseInfo(&nodes, &elements[0], &edges[0], &elementsTran[0], &edgesTran[0], &isInNode, &nodeCEdge);
	if (IsExist==0) {
		printf("Error: constructing coarse grid fails!\n");
		return RUN_FAIL;
	}
	
	// refinements
	clock_t refine_start=clock();	
	
	for (i=0;i<levelNum-1;++i) {	
		fasp_ivec_free(&isInNode);
		refine(&nodes, &elements[i], &edges[i], &elementsTran[i], &elements[i+1], &edges[i+1], &elementsTran[i+1], &edgesTran[i+1], &isInNode, &nodeCEdge);
	}
	
	clock_t refine_end=clock();
	
	double refine_time = (double)(refine_end - refine_start)/(double)(CLOCKS_PER_SEC);
	printf("Refinement costs %f seconds.\n", refine_time); 
	
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
	
	clock_t assemblestif_start=clock();
	
	assemble_stiffmat(&elements[levelNum-1], &elementsTran[levelNum-1], &nodes, &A, &edgesTran[levelNum-1]);
	
	clock_t assemblestif_end=clock();
	
	double assemblestif_time = (double)(assemblestif_end - assemblestif_start)/(double)(CLOCKS_PER_SEC);
	
	printf("Assemble A costs %f seconds.\n", assemblestif_time); 
	
	// set initial boundary value
	dvector uh = fasp_dvec_create(A.row);
	
	for (i=0;i<uh.row;++i)
	{
		if(isInNode.val[i]==-1) // the node is on the boundary
		{ 
			uh.val[i]=u(nodes.val[i][0], nodes.val[i][1]);
		}
	}
	
	clock_t assemblerhs_start=clock();
	
	// assemble right-hand side
	dvector b = fasp_dvec_create(A.row);
	
	double nodetemp[3][2],btemp[3];
	for (l=0;l<elements[levelNum-1].row;l++)
	{
		i=elements[levelNum-1].JA[elements[levelNum-1].IA[l]];
		j=elements[levelNum-1].JA[elements[levelNum-1].IA[l]+1];
		k=elements[levelNum-1].JA[elements[levelNum-1].IA[l]+2];
		nodetemp[0][0]=nodes.val[i][0];
		nodetemp[0][1]=nodes.val[i][1];
		nodetemp[1][0]=nodes.val[j][0];
		nodetemp[1][1]=nodes.val[j][1];
		nodetemp[2][0]=nodes.val[k][0];
		nodetemp[2][1]=nodes.val[k][1];
		
		localb(nodetemp,btemp);
		
		b.val[i]+=btemp[0];
		b.val[j]+=btemp[1];
		b.val[k]+=btemp[2];
	}
	
	clock_t assemblerhs_end=clock();
	
	double assemblerhs_time = (double)(assemblerhs_end - assemblerhs_start)/(double)(CLOCKS_PER_SEC);
	
	printf("Assemble b costs %f seconds.\n", assemblerhs_time); 
	
	extractNondirichletMatrix(&A, &b, ptr_A, ptr_b, &isInNode, &dirichlet, &nondirichlet, &index, &uh);
	
	// clean up memory
	for (i=0;i<nodes.row;++i) fasp_mem_free(nodes.val[i]);
	fasp_mem_free(nodes.val);
	
	for (l=0;l<levelNum;l++)
	{
		fasp_icsr_free(&elements[l]);
		fasp_icsr_free(&elementsTran[l]);
		fasp_icsr_free(&edgesTran[l]);
		fasp_iden_free(&edges[l]);
	}
	
	fasp_dcsr_free(&A);
	fasp_dvec_free(&b);
	fasp_dvec_free(&uh);
	fasp_ivec_free(&dirichlet);
	fasp_ivec_free(&nondirichlet);
	fasp_ivec_free(&index);
	fasp_ivec_free(&isInNode);
	fasp_ivec_free(&nodeCEdge);
	
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
