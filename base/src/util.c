/*! \file util.c
 *  \brief Useful utilities.
 */

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn unsigned long fasp_aux_change_endian4 (unsigned long x)
 *
 * \brief Swap order for different endian systems
 *
 * \param x   An unsigned long integer
 *
 * \return    Unsigend long ineger after swapping
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
unsigned long fasp_aux_change_endian4(unsigned long x)
{
    unsigned char *ptr = (unsigned char *)&x;
    return (ptr[0] << 24) | (ptr[1] << 16) | (ptr[2] << 8) | ptr[3];
}

/**
 * \fn REAL fasp_aux_change_endian8 (REAL x)
 *
 * \brief Swap order for different endian systems
 *
 * \param x   A unsigned long integer
 *
 * \return    Unsigend long ineger after swapping
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
REAL fasp_aux_change_endian8(REAL x)
{   
    REAL dbl;   
    unsigned char *bytes, *buffer;   
    
    buffer=(unsigned char *)&dbl;
    bytes=(unsigned char *)&x;
    
    buffer[0]=bytes[7];   
    buffer[1]=bytes[6];   
    buffer[2]=bytes[5];   
    buffer[3]=bytes[4];   
    buffer[4]=bytes[3];   
    buffer[5]=bytes[2];   
    buffer[6]=bytes[1];   
    buffer[7]=bytes[0];   
    return dbl;   
}

/**
 * \fn REAL fasp_aux_bbyteTolREAL (unsigned char bytes[])   
 *
 * \brief Swap order for different endian systems
 *
 * \param bytes  A unsigned char
 *
 * \return       Unsigend long ineger after swapping
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
REAL fasp_aux_bbyteTolREAL(unsigned char bytes[])
{   
    REAL dbl;   
    unsigned char *buffer;   
    buffer=(unsigned char *)&dbl;   
    buffer[0]=bytes[7];   
    buffer[1]=bytes[6];   
    buffer[2]=bytes[5];   
    buffer[3]=bytes[4];   
    buffer[4]=bytes[3];   
    buffer[5]=bytes[2];   
    buffer[6]=bytes[1];   
    buffer[7]=bytes[0];   
    return dbl;   
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
