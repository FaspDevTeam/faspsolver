/*! \file testfct_poisson.c
 *
 *  \brief Test functions for the Poisson's equation
 */
 
#include <math.h>

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn double f (double x, double y)
 *
 * \brief Load f of 2D Poisson's equation: - u_xx - u_yy = f
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value
 *
 * \note Right hand side when the true solution u is -cos(x)*sin(y)
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 */
double f (double x, double y)
{
	return -2*cos(x)*sin(y);
}

/**
 * \fn double u (double x, double y)
 *
 * \brief Exact solution u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    function value 
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 */
double u (double x, double y)
{
	return -cos(x)*sin(y);
}

/**
 * \fn double u_x (double x, double y)
 *
 * \brief X-directional partial derivative of u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    x-directional partial derivative of true solution u 
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 */
double u_x (double x, double y)
{
	return sin(x)*sin(y);
}

/**
 * \fn double u_y(double x, double y)
 *
 * \brief Y-directional partial derivative of u
 *
 * \param x   the x-axis value of the point
 * \param y   the y-axis value of the point
 * \return    y-directional partial derivative of true solution u 
 *
 * \author Xuehai Huang
 * \date   03/29/2009
 */
double u_y (double x, double y)
{
	return -cos(x)*cos(y);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
