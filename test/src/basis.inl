/*! \file basis.inl
 *  \brief Finite Element Basis Functions
 */
 
#include <stdio.h>
#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/** 
 * \fn void basisP1(double nodes[3][2], double s, int index, double phi[2])
 * \brief Basis function of Lagrange element
 *
 * \param nodes[3][2]   the vertice of the triangule
 * \param s             the area of the triangule
 * \param index         the indicator of the basis function
 * \param phi[2]        basis function
 * \return void
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 */
void basisP1(double nodes[3][2], double s, int index, double phi[2])
{
  const int node1 = (index+1)%3, node2 = (index+2)%3;
  phi[0]=(nodes[node1][1]-nodes[node2][1])/(2.0*s);
  phi[1]=(nodes[node2][0]-nodes[node1][0])/(2.0*s);
}

/**
 * \fn double areaT(double x1,double x2,double x3,double y1,double y2,double y3)
 * \brief Get area for triangle p1(x1,y1),p2(x2,y2),p3(x3,y3)
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
 * \date   03/29/2009
 */
double areaT(double x1,double x2,double x3,double y1,double y2,double y3)
{
  return ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))/2;
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
