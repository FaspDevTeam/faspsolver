/*! \file basis.c
 *  \brief Basis functions and problem information
 */
 
#include <stdio.h>
#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn double f(double x, double y)
 * \brief load f of 2D Laplace equation: u_xx+u_yy=-f
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value
 *
 * \note Right hand side when the true solution u is -cos(x)*sin(y)
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
double f(double x, double y)
{
	return -2*cos(x)*sin(y);
}

/**
 * \fn double u(double x, double y)
 * \brief true solution u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value 
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
double u(double x, double y)
{
	return -cos(x)*sin(y);
}

/**
 * \fn double u_x(double x, double y)
 * \brief the x-directional partial derivative of true solution u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    x-directional partial derivative of true solution u 
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
double u_x(double x, double y)
{
	return sin(x)*sin(y);
}

/**
 * \fn double u_y(double x, double y)
 * \brief the y-directional partial derivative of true solution u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    y-directional partial derivative of true solution u 
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
double u_y(double x, double y)
{
	return -cos(x)*cos(y);
}

/** 
 * \fn void basis(double nodes[3][2], double s, int index, double phi[2])
 * \brief basis function of Lagrange element, i.e. area coordiante
 *
 * \param nodes[3][2]   the vertice of the triangule
 * \param s             the area of the triangule
 * \param index         the indicator of the basis function
 * \param phi[2]        basis function
 * \return void
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
void basis(double nodes[3][2], double s, int index, double phi[2])
{
	const int node1 = (index+1)%3, node2 = (index+2)%3;
	phi[0]=(nodes[node1][1]-nodes[node2][1])/(2.0*s);
	phi[1]=(nodes[node2][0]-nodes[node1][0])/(2.0*s);
}

/**
 * \fn double area(double x1,double x2,double x3,double y1,double y2,double y3)
 * \brief get area for triangle p1(x1,y1),p2(x2,y2),p3(x3,y3)
 *
 * \param x1   the x-axis value of the point p1
 * \param x2   the x-axis value of the point p2
 * \param x3   the x-axis value of the point p3
 * \param y1   the y-axis value of the point p1
 * \param y2   the y-axis value of the point p2
 * \param y3   the y-axis value of the point p3
 * \return     area of the trianle
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
double area(double x1,double x2,double x3,double y1,double y2,double y3)
{
	return ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))/2;
}

/** 
 * \fn void localb(double (*nodes)[2],double *b)
 * \brief get local right-hand side b from triangle nodes
 *
 * \param (*nodes)[2]  the vertice of the triangule
 * \param *b           local right-hand side
 *
 * \author Xuehai Huang
 * \date 03/29/2009
 */
void localb(double (*nodes)[2],double *b)
{
	const double s=2.0*area(nodes[0][0],nodes[1][0],nodes[2][0],nodes[0][1],nodes[1][1],nodes[2][1]);
	const int num_qp=16; // the number of numerical intergation points
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

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
