/*! \file convert.c
 *
 *  \brief Some utilities for format conversion.
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
unsigned long fasp_aux_change_endian4 (unsigned long x)
{
    unsigned char *ptr = (unsigned char *)&x;
    return (ptr[0] << 24) | (ptr[1] << 16) | (ptr[2] << 8) | ptr[3];
}

/**
 * \fn double fasp_aux_change_endian8 (double x)
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
double fasp_aux_change_endian8 (double x)
{
    double dbl;
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
 * \fn double fasp_aux_bbyteToldouble (unsigned char bytes[])
 *
 * \brief Swap order of double-precision float for different endian systems
 *
 * \param bytes  A unsigned char
 *
 * \return       Unsigend long ineger after swapping
 *
 * \author Chensong Zhang
 * \date   11/16/2009
 */
double fasp_aux_bbyteToldouble (unsigned char bytes[])
{
    double dbl;
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

/**
 * \fn INT endian_convert_int (const INT inum, INT ilength, INT endianflag)
 *
 * \brief Swap order of an INT number
 *
 * \param index       An INT value
 * \param ilength     Length of INT: 2 for short, 4 for int, 8 for long
 * \param endianflag  If endianflag = 1, it returns inum itself
 *                    If endianflag = 2, it returns the swapped inum
 *
 * \return Value of inum or swapped inum
 *
 * \author Ziteng Wang
 * \date   2012-12-24
 */
INT endian_convert_int (const INT inum,
                        const INT ilength,
                        const INT endianflag)
{
    INT iretVal,i;
    char *intToConvert = ( char* ) & inum;
    char *returnInt = ( char* ) & iretVal;
    
    if (endianflag==1) return inum;
    else {
        for (i = 0; i < ilength; i++) {
            returnInt[i] = intToConvert[ilength-i-1];
        }
        return iretVal;
    }
}

/**
 * \fn REAL endian_convert_real (const REAL rnum, INT ilength, INT endianflag)
 *
 * \brief Swap order of a REAL number
 *
 * \param index       An REAL value
 * \param ilength     Length of INT: 2 for short, 4 for int, 8 for long
 * \param endianflag  If endianflag = 1, it returns rnum itself
 *                    If endianflag = 2, it returns the swapped rnum
 *
 * \return Value of rnum or swapped rnum
 *
 * \author Ziteng Wang
 * \date   2012-12-24
 */
REAL endian_convert_real (const REAL rnum,
                          INT vlength,
                          INT endianflag)
{
    REAL dretVal;
    char *realToConvert = (char *) & rnum;
    char *returnReal    = (char *) & dretVal;
	INT  i;
    
    if (endianflag==1) return rnum;
    else {
        for (i = 0; i < vlength; i++) {
            returnReal[i] = realToConvert[vlength-i-1];
        }
        return dretVal;
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
