/*! \file ilu.c
 *
 *  \brief Incomplete LU decomposition: ILUk, ILUt, ILUtp
 *
 *  \note  This is a translation from SPARSEKIT Fortran version
 *
 *  Translated by Chunsheng Feng, 09/03/2016
 */

#include <math.h>
#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

/**
 * \fn void void fasp_qsplit (REAL *a, INT *ind, INT n, INT ncut)
 *
 * \brief Get  a quick-sort split of a real array
 *
 * \param a  a real array. on output a(1:n) is permuted such that
 *           its elements satisfy: abs(a(i)) .ge. abs(a(ncut)) for
 *           i .lt. ncut and abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut.
 * \param ind  is an integer array which permuted in the same way as a(*).
 * \param n  size of array a.
 * \param ncut  integer.
 *
 * \author Chunsheng Feng
 * \date   09/06/2016
 */
void fasp_qsplit (REAL *a,
                  INT *ind,
                  INT n,
                  INT ncut)
{
    /*-----------------------------------------------------------------------
     does a quick-sort split of a real array.
     on input a(1:n). is a real array
     on output a(1:n) is permuted such that its elements satisfy:
     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
     ind(1:n) is an integer array which permuted in the same way as a(*).
     -----------------------------------------------------------------------*/
    REAL tmp, abskey;
    INT itmp, first, last, mid, j;
    
    /* Parameter adjustments */
    --ind;
    --a;
    
    first = 1;
    last = n;
    if ((ncut  <  first) || (ncut  >  last)) return;
    
    //    outer loop -- while mid .ne. ncut do
F161:
    mid = first;
    abskey = ABS(a[mid]);
    for (j = first + 1; j <= last; ++j ) {
        if (ABS(a[j])  >  abskey) {
            ++mid;
            //     interchange
            tmp = a[mid];
            itmp = ind[mid];
            a[mid] = a[j];
            ind[mid] = ind[j];
            a[j] = tmp;
            ind[j] = itmp;
        }
    }
    
    //     interchange
    tmp = a[mid];
    a[mid] = a[first];
    a[first] = tmp;
    //
    itmp = ind[mid];
    ind[mid] = ind[first];
    ind[first] = itmp;
    
    //    test for while loop
    if (mid  ==  ncut) {
        ++ind;
        ++a;
        return;
    }
    
    if (mid  >  ncut)  {
        last = mid - 1;
    } else   {
        first = mid + 1;
    }
    
    goto F161;
    /*----------------end-of-qsplit------------------------------------------*/
}

/**
 * \fn void fasp_iluk (INT n, REAL *a,INT *ja, INT *ia, INT lfil,
 *                     REAL *alu, INT *jlu, INT iwk, INT *ierr, INT *nzlu)
 *
 * \brief Get ILU factorization with level of fill-in k (ilu(k)) for a CSR matrix A
 *
 * \param n   row number of A
 * \param a   nonzero entries of A
 * \param ja  integer array of column for A
 * \param ia  integer array of row pointers for A
 * \param lfil  integer. The fill-in parameter. Each row of L and each row
 *              of U will have a maximum of lfil elements (excluding the diagonal
 *              element). lfil must be .ge. 0.
 * \param droptol  real*8. Sets the threshold for dropping small terms in the
 *                 factorization. See below for details on dropping strategy.
 * \param alu,jlu  matrix stored in Modified Sparse Row (MSR) format containing
 *                 the L and U factors together. The diagonal (stored in
 *                 alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
 *                 contains the i-th row of L (excluding the diagonal entry=1)
 *                 followed by the i-th row of U.
 * \param jlu  integer array of length n containing the pointers to
 *             the beginning of each row of U in the matrix alu,jlu.
 * \param iwk  integer. The minimum length of arrays alu, jlu, and levs.
 * \param ierr  integer pointer. Return error message with the following meaning.
 *                0  --> successful return.
 *               >0  --> zero pivot encountered at step number ierr.
 *               -1  --> Error. input matrix may be wrong.
 *                       (The elimination process has generated a
 *                       row in L or U whose length is .gt.  n.)
 *               -2  --> The matrix L overflows the array al.
 *               -3  --> The matrix U overflows the array alu.
 *               -4  --> Illegal value for lfil.
 *               -5  --> zero row encountered.
 * \param nzlu  integer pointer. Return number of nonzero entries for alu and jlu
 *
 * \note:  All the diagonal elements of the input matrix must be nonzero.
 *
 * \author Chunsheng Feng
 * \date   09/06/2016
 */
void fasp_iluk (INT n,
                REAL *a,
                INT *ja,
                INT *ia,
                INT lfil,
                REAL *alu,
                INT *jlu,
                INT iwk,
                INT *ierr,
                INT *nzlu)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: %s (ILUk) ...... [Start]\n", __FUNCTION__);
#endif
    /*----------------------------------------------------------------------
     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k))
     ----------------------------------------------------------------------
     
     on entry:
     ==========
     n       = integer. The row dimension of the matrix A. The matrix
     
     a,ja,ia = matrix stored in Compressed Sparse Row format.
     
     lfil    = integer. The fill-in parameter. Each element whose
     leve-of-fill exceeds lfil during the ILU process is dropped.
     lfil must be .ge. 0
     
     iwk     = integer. The minimum length of arrays alu, jlu, and levs.
     
     On return:
     ===========
     
     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
     the L and U factors together. The diagonal (stored in
     alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
     contains the i-th row of L (excluding the diagonal entry=1)
     followed by the i-th row of U.
     
     jlu     = integer array of length n containing the pointers to
     the beginning of each row of U in the matrix alu,jlu.
     
     levs    = integer (work) array of size iwk -- which contains the
     levels of each element in alu, jlu.
     
     ierr    = integer. Error message with the following meaning.
     ierr  = 0    --> successful return.
     ierr .gt. 0  --> zero pivot encountered at step number ierr.
     ierr  = -1   --> Error. input matrix may be wrong.
     (The elimination process has generated a
     row in L or U whose length is .gt.  n.)
     ierr  = -2   --> The matrix L overflows the array al.
     ierr  = -3   --> The matrix U overflows the array alu.
     ierr  = -4   --> Illegal value for lfil.
     ierr  = -5   --> zero row encountered in A or U.
     
     work arrays:
     =============
     jw      = integer work array of length 3*n.
     w       = real work array of length n
     
     ----------------------------------------------------------------------
     w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = U-part]
     jw(n+1:2n)  stores the nonzero indicator.
     
     Notes:
     ------
     All the diagonal elements of the input matrix must be nonzero.
     
     ---------------------------------------------------------------------- */
    
    //     locals
    INT ju0, k, j1, j2, j, ii, i, lenl, lenu, jj, jrow, jpos, n2, jlev, NE;
    REAL t, s, fact;
    SHORT cinindex=0;
    REAL *w;
    INT *ju, *jw, *levs;
    
    if (lfil  <  0) goto F998;
    
    w = (REAL *)fasp_mem_calloc(n, sizeof(REAL));
    ju = (INT *)fasp_mem_calloc(n, sizeof(INT));
    jw = (INT *)fasp_mem_calloc(3*n, sizeof(INT));
    levs = (INT *)fasp_mem_calloc(iwk, sizeof(INT));
    
    --jw;
    --w;
    --ju;
    --jlu;
    --alu;
    --ia;
    --ja;
    --a;
    --levs;
    
    /*-----------------------------------------------------------------------
     shift index for C routines
     -----------------------------------------------------------------------*/
    if (ia[1]  ==  0) cinindex=1 ;
    if (cinindex)
    {
        NE = n + 1; //modify by chunsheng 2012, Sep, 1;
        for (i=1; i<=NE; ++i)  ++ia[i];
        NE = ia[n+1] - 1;
        for (i=1; i<=NE; ++i)  ++ja[i];
    }
    
    /*-----------------------------------------------------------------------
     initialize ju0 (points to next element to be added to alu,jlu)
     and pointer array.
     -----------------------------------------------------------------------  */
    n2 = n + n;
    ju0 = n + 2;
    jlu[1] = ju0;
    
    // initialize nonzero indicator array + levs array --
    for(j = 1; j<=2*n; ++j) jw[j] = 0;
    
    /*-----------------------------------------------------------------------
     beginning of main loop.
     ----------------------------------------------------------------------- */
    for(ii = 1; ii <= n; ++ii) 	{  //500
        j1 = ia[ii];
        j2 = ia[ii + 1] - 1;
        
        //	 unpack L-part and U-part of row of A in arrays w
        lenu = 1;
        lenl = 0;
        jw[ii] = ii;
        w[ii] = 0.0;
        jw[n + ii] = ii;
        
        //
        for(j = j1; j <= j2; ++j) 	{ //170
            k = ja[j];
            t = a[j];
            if (t  ==  0.0) continue;  //goto g170;
            if (k  <  ii) 	{
                ++lenl;
                jw[lenl] = k;
                w[lenl] = t;
                jw[n2 + lenl] = 0;
                jw[n + k] = lenl;
            } else if (k  ==  ii) {
                w[ii] = t;
                jw[n2 + ii] = 0;
            } else 	{
                ++lenu;
                jpos = ii + lenu - 1;
                jw[jpos] = k;
                w[jpos] = t;
                jw[n2 + jpos] = 0;
                jw[n + k] = jpos;
            }
            
        }  //170
        
        jj = 0;
        //	 eliminate previous rows
        
    F150:
        ++jj;
        if (jj  >  lenl) goto F160;
        
        /*-----------------------------------------------------------------------
         in order to do the elimination in the correct order we must select
         the smallest column index among jw(k), k=jj+1, ..., lenl.
         -----------------------------------------------------------------------*/
        
        jrow = jw[jj];
        k = jj;
        
        //	 determine smallest column index
        for(j = jj + 1; j <= lenl; ++j) 	{ //151
            if (jw[j]  <  jrow) {
                jrow = jw[j];
                k = j;
            }
        } //151
        
        if (k  !=  jj) {
            //     exchange in jw
            j = jw[jj];
            jw[jj] = jw[k];
            jw[k] = j;
            //     exchange in jw(n+  (pointers/ nonzero indicator).
            jw[n + jrow] = jj;
            jw[n + j] = k;
            //     exchange in jw(n2+  (levels)
            j = jw[n2 + jj];
            jw[n2 + jj] = jw[n2 + k];
            jw[n2 + k] = j;
            //     exchange in w
            s = w[jj];
            w[jj] = w[k];
            w[k] = s;
        }
        
        //	 zero out element in row by resetting jw(n+jrow) to zero.
        jw[n + jrow] = 0;
        
        //	 get the multiplier for row to be eliminated (jrow) + its level
        fact = w[jj]*alu[jrow];
        jlev = jw[n2 + jj];
        if (jlev  >  lfil) goto F150;
        
        //	 combine current row and row jrow
        for(k = ju[jrow]; k <= jlu[jrow + 1] - 1; ++k ) { // 203
            s = fact*alu[k];
            j = jlu[k];
            jpos = jw[n + j];
            if (j  >=  ii) {
                //	 dealing with upper part.
                if (jpos  ==  0) {
                    //	 this is a fill-in element
                    ++lenu;
                    if (lenu  >  n) goto F995;
                    i = ii + lenu - 1;
                    jw[i] = j;
                    jw[n + j] = i;
                    w[i] = -s;
                    jw[n2 + i] = jlev + levs[k] + 1;
                } else 	{
                    //	 this is not a fill-in element
                    w[jpos] = w[jpos] - s;
                    jw[n2 + jpos] = MIN(jw[n2 + jpos], jlev + levs[k] + 1);
                }
            } else 	{
                //	 dealing with lower part.
                if (jpos  ==  0) 	{
                    //	 this is a fill-in element
                    ++lenl;
                    if (lenl  >  n) goto F995;
                    jw[lenl] = j;
                    jw[n + j] = lenl;
                    w[lenl] = -s;
                    jw[n2 + lenl] = jlev + levs[k] + 1;
                } else {
                    //	 this is not a fill-in element
                    w[jpos] = w[jpos] - s;
                    jw[n2 + jpos] = MIN(jw[n2 + jpos], jlev + levs[k] + 1);
                }
            }
            
        } //203
        w[jj] = fact;
        jw[jj] = jrow;
        goto F150;
        
    F160:
        //  reset double-pointer to zero (U-part)
        for(k = 1; k <= lenu; ++k)  jw[n + jw[ii + k - 1]] = 0;
        
        //	 update l-matrix
        for(k = 1; k <= lenl; ++k ) {   //204
            if (ju0  >  iwk) goto F996;
            if (jw[n2 + k]  <=  lfil)  {
                alu[ju0] = w[k];
                jlu[ju0] = jw[k];
                ++ju0;
            }
        } //204
        
        //	 save pointer to beginning of row ii of U
        ju[ii] = ju0;
        
        //	 update u-matrix
        for(k = ii + 1; k <= ii + lenu - 1; ++k ) {  //302
            if (ju0  >  iwk) goto F997;
            
            if (jw[n2 + k]  <=  lfil) {
                jlu[ju0] = jw[k];
                alu[ju0] = w[k];
                levs[ju0] = jw[n2 + k];
                ++ju0;
            }
            
        } //302
        
        if (w[ii]  ==  0.0) goto F999;
        //
        alu[ii] = 1.0/w[ii];
        
        //	 update pointer to beginning of next row of U.
        jlu[ii + 1] = ju0;
        /*-----------------------------------------------------------------------
         end main loop
         -----------------------------------------------------------------------*/
    } //500
    
    *nzlu = ju[n] - 1;
    
    if (cinindex)  {
        for ( i = 1; i <= *nzlu; ++i ) 	--jlu[i];
    }
    
    *ierr = 0;
    
F100:
    ++jw;
    ++w;
    ++ju;
    ++jlu;
    ++alu;
    ++ia;
    ++ja;
    ++a;
    ++levs;
    
    fasp_mem_free(w);
    fasp_mem_free(ju);
    fasp_mem_free(jw);
    fasp_mem_free(levs);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return;
    
    //	 incomprehensible error. Matrix must be wrong.
F995:
    printf("### ERROR: incomprehensible error. Matrix must be wrong.\n");
    *ierr = -1;
    goto F100;
    
    // insufficient storage in L.
F996:
    printf("### ERROR: insufficient storage in L.\n");
    *ierr = -2;
    goto F100;
    
    // insufficient storage in U.
F997:
    printf("### ERROR: insufficient storage in U.\n");
    *ierr = -3;
    goto F100;
    
    //  illegal lfil entered.
F998:
    printf("### ERROR: illegal lfil entered\n");
    *ierr = -4;
    return;
    
    // zero row encountered in A or U.
F999:
    printf("### ERROR: zero row encountered in A or U.\n");
    *ierr = -5;
    goto F100;
    /*----------------end-of-iluk--------------------------------------------
     ---------------------------------------------------------------------- */
}

/**
 * \fn void fasp_ilut (INT n, REAL *a, INT *ja, INT *ia, INT lfil, REAL droptol,
 *                     REAL *alu, INT *jlu, INT iwk, INT *ierr, INT *nz)
 *
 * \brief Get incomplete LU factorization with dual truncations of a CSR matrix A
 *
 * \param n   row number of A
 * \param a   nonzero entries of A
 * \param ja  integer array of column for A
 * \param ia  integer array of row pointers for A
 * \param lfil  integer. The fill-in parameter. Each row of L and each row
 *              of U will have a maximum of lfil elements (excluding the diagonal
 *              element). lfil must be .ge. 0.
 * \param droptol  real*8. Sets the threshold for dropping small terms in the
 *                 factorization. See below for details on dropping strategy.
 * \param alu,jlu  matrix stored in Modified Sparse Row (MSR) format containing
 *                 the L and U factors together. The diagonal (stored in
 *                 alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
 *                 contains the i-th row of L (excluding the diagonal entry=1)
 *                 followed by the i-th row of U.
 * \param iwk  integer. The lengths of arrays alu and jlu. If the arrays
 *             are not big enough to store the ILU factorizations, ilut
 *             will stop with an error message.
 * \param ierr  integer pointer. Return error message with the following meaning.
 *                0  --> successful return.
 *               >0  --> zero pivot encountered at step number ierr.
 *               -1  --> Error. input matrix may be wrong.
 *                       (The elimination process has generated a
 *                       row in L or U whose length is .gt.  n.)
 *               -2  --> The matrix L overflows the array al.
 *               -3  --> The matrix U overflows the array alu.
 *               -4  --> Illegal value for lfil.
 *               -5  --> zero row encountered.
 * \param nz  integer pointer. Return number of nonzero entries for alu and jlu
 *
 * \note  All the diagonal elements of the input matrix must be nonzero.
 *
 * \author Chunsheng Feng
 * \date   09/06/2016
 */
void fasp_ilut (INT n,
                REAL *a,
                INT  *ja,
                INT  *ia,
                INT lfil,
                REAL droptol,
                REAL *alu,
                INT *jlu,
                INT iwk,
                INT *ierr,
                INT *nz)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: %s (ILUt) ...... [Start]\n", __FUNCTION__);
#endif
    /*--------------------------------------------------------------------*
     *** ILUT preconditioner ***                                      *
     incomplete LU factorization with dual truncation mechanism       *
     ----------------------------------------------------------------------*
     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
     ----------------------------------------------------------------------*
     PARAMETERS
     -----------
     
     on entry:
     ==========
     n       = integer. The row dimension of the matrix A. The matrix
     
     a,ja,ia = matrix stored in Compressed Sparse Row format.
     
     lfil    = integer. The fill-in parameter. Each row of L and each row
     of U will have a maximum of lfil elements (excluding the diagonal
     element). lfil must be .ge. 0.
     
     droptol = real*8. Sets the threshold for dropping small terms in the
     factorization. See below for details on dropping strategy.
     
     iwk     = integer. The lengths of arrays alu and jlu. If the arrays
     are not big enough to store the ILU factorizations, ilut
     will stop with an error message.
     
     On return:
     ===========
     
     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
     the L and U factors together. The diagonal (stored in
     alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
     contains the i-th row of L (excluding the diagonal entry=1)
     followed by the i-th row of U.
     
     ju      = integer array of length n containing the pointers to
     the beginning of each row of U in the matrix alu,jlu.
     
     ierr    = integer. Error message with the following meaning.
     ierr  = 0    --> successful return.
     ierr .gt. 0  --> zero pivot encountered at step number ierr.
     ierr  = -1   --> Error. input matrix may be wrong.
     (The elimination process has generated a
     row in L or U whose length is .gt.  n.)
     ierr  = -2   --> The matrix L overflows the array al.
     ierr  = -3   --> The matrix U overflows the array alu.
     ierr  = -4   --> Illegal value for lfil.
     ierr  = -5   --> zero row encountered.
     
     work arrays:
     =============
     jw      = integer work array of length 2*n.
     w       = real work array of length n+1.
     
     ----------------------------------------------------------------------
     w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
     jw(n+1:2n)  stores nonzero indicators
     
     Notes:
     ------
     The diagonal elements of the input matrix must be nonzero (at least
     'structurally').
     
     ----------------------------------------------------------------------*
     ---- Dual drop strategy works as follows.                             *
     *
     1) Theresholding in L and U as set by droptol. Any element whose *
     magnitude is less than some tolerance (relative to the abs       *
     value of diagonal element in u) is dropped.                      *
     *
     2) Keeping only the largest lfil elements in the i-th row of L   *
     and the largest lfil elements in the i-th row of U (excluding    *
     diagonal elements).                                              *
     *
     Flexibility: one can use droptol=0 to get a strategy based on    *
     keeping the largest elements in each row of L and U. Taking      *
     droptol .ne. 0 but lfil=n will give the usual threshold strategy *
     (however, fill-in is then unpredictible).                        *
     ---------------------------------------------------------------------- */
    
    //     locals
    INT ju0, k, j1, j2, j, ii, i, lenl, lenu, jj, jrow, jpos, NE, len;
    REAL t, s, fact, tmp;
    SHORT cinindex=0;
    REAL *w, *tnorm;
    INT  *ju, *jw;
    
    if (lfil  <  0) goto F998;
    
    ju = (INT *)fasp_mem_calloc(n, sizeof(INT));
    jw = (INT *)fasp_mem_calloc(2*n, sizeof(INT));
    w = (REAL *)fasp_mem_calloc(n+1, sizeof(REAL));
    tnorm = (REAL *)fasp_mem_calloc(n, sizeof(REAL));
    
    --jw;
    --ju;
    --w;
    --tnorm;
    --jlu;
    --alu;
    --ia;
    --ja;
    --a;
    
    if (ia[1]  ==  0) cinindex=1 ;
    
    if (cinindex)
    {
        NE = n + 1; //modify by chunsheng 2012, Sep, 1;
        for (i=1; i<=NE; ++i)  ++ia[i];
        NE = ia[n+1] - 1;
        for (i=1; i<=NE; ++i)  ++ja[i];
    }
    
    /*-----------------------------------------------------------------------
     initialize ju0 (points to next element to be added to alu,jlu)
     and pointer array.
     -----------------------------------------------------------------------*/
    ju0 = n + 2;
    jlu[1] = ju0;
    
    // initialize nonzero indicator array.
    for (j = 1; j<=n; ++j)  jw[n + j] = 0;
    
    
    /*-----------------------------------------------------------------------
     beginning of main loop.
     -----------------------------------------------------------------------*/
    for (ii = 1; ii <= n; ++ii ) {
        j1 = ia[ii];
        j2 = ia[ii + 1] - 1;
        tmp = 0.0;
        for ( k = j1; k<= j2; ++k) 	tmp = tmp + ABS(a[k]);
        tmp = tmp/(REAL)(j2 - j1 + 1);
        tnorm[ii] = tmp*droptol;;
    }
    
    for (ii = 1; ii<=n; ++ii) {
        j1 = ia[ii];
        j2 = ia[ii + 1] - 1;
        
        //   unpack L-part and U-part of row of A in arrays w
        lenu = 1;
        lenl = 0;
        jw[ii] = ii;
        w[ii] = 0.0;
        jw[n + ii] = ii;
        
        for(j = j1; j<=j2; ++j) {
            k = ja[j];
            t = a[j];
            if (k  <  ii) {
                ++lenl;
                jw[lenl] = k;
                w[lenl] = t;
                jw[n + k] = lenl;
            } else if (k == ii) {
                w[ii] = t;
            } else 	{
                ++lenu ;
                jpos = ii + lenu - 1;
                jw[jpos] = k;
                w[jpos] = t;
                jw[n + k] = jpos;
            }
        }
        jj = 0;
        len = 0;
        
        //     eliminate previous rows
    F150:
        ++jj;
        if (jj  >  lenl) goto F160;
        
        /*-----------------------------------------------------------------------
         in order to do the elimination in the correct order we must select
         the smallest column index among jw(k), k=jj+1, ..., lenl.
         -----------------------------------------------------------------------*/
        jrow = jw[jj];
        k = jj;
        
        /*
         determine smallest column index
         */
        for(j = jj + 1; j<=lenl; ++j)	{  //151
            if (jw[j]  <  jrow) {
                jrow = jw[j];
                k = j;
            }
        }    //151
        
        if (k  !=  jj) {
            // exchange in jw
            j = jw[jj];
            jw[jj] = jw[k];
            jw[k] = j;
            // exchange in jr
            jw[n + jrow] = jj;
            jw[n + j] = k;
            // exchange in w
            s = w[jj];
            w[jj] = w[k];
            w[k] = s;
        }
        
        //     zero out element in row by setting jw(n+jrow) to zero.
        jw[n + jrow] = 0;
        
        //    get the multiplier for row to be eliminated (jrow).
        fact = w[jj]*alu[jrow];
        
        if (ABS(fact)  <=  droptol) goto F150;
        
        //     combine current row and row jrow
        for ( k = ju[jrow]; k <= jlu[jrow + 1] - 1; ++k) {   //203
            s = fact*alu[k];
            j = jlu[k];
            jpos = jw[n + j];
            if (j  >=  ii)  {
                //     dealing with upper part.
                if (jpos  ==  0)
                {
                    //     this is a fill-in element
                    ++lenu;
                    if (lenu  >  n) goto F995;
                    i = ii + lenu - 1;
                    jw[i] = j;
                    jw[n + j] = i;
                    w[i] = -s;
                } else 	{
                    //    this is not a fill-in element
                    w[jpos] = w[jpos] - s;
                }
            } else {
                //     dealing  with lower part.
                if (jpos  ==  0) {
                    //     this is a fill-in element
                    ++lenl;
                    if (lenl  >  n) goto F995;
                    jw[lenl] = j;
                    jw[n + j] = lenl;
                    w[lenl] = -s;
                } else 	{
                    //    this is not a fill-in element
                    w[jpos] = w[jpos] - s;
                }
            }
        }  //203
        
        /*
         store this pivot element -- (from left to right -- no danger of
         overlap with the working elements in L (pivots).
         */
        ++len;
        w[len] = fact;
        jw[len] = jrow;
        goto F150;
        
    F160:
        //     reset double-pointer to zero (U-part)
        for (k = 1; k <= lenu; ++k ) jw[n + jw[ii + k - 1]] = 0;  //308
        
        //     update L-matrix
        lenl = len;
        len = MIN(lenl, lfil);
        
        //     sort by quick-split
        fasp_qsplit(&w[1], &jw[1], lenl, len);
        
        //     store L-part
        for (k = 1; k <= len; ++k ) 	{   //204
            if (ju0  >  iwk) goto F996;
            alu[ju0] = w[k];
            jlu[ju0] = jw[k];
            ++ju0;
        }
        
        //     save pointer to beginning of row ii of U
        ju[ii] = ju0;
        
        //     update U-matrix -- first apply dropping strategy
        len = 0;
        for (k = 1; k <= lenu - 1; ++k) {
            //		if ( ABS(w[ii + k])  >  droptol*tnorm )
            if ( ABS(w[ii + k])  >  tnorm[ii] ) {
                ++len;
                w[ii + len] = w[ii + k];
                jw[ii + len] = jw[ii + k];
            }
        }
        
        lenu = len + 1;
        len = MIN(lenu, lfil);
        
        fasp_qsplit(&w[ii + 1], &jw[ii + 1], lenu - 1, len);
        
        //     copy
        t = ABS(w[ii]);
        if (len + ju0  >  iwk) goto F997;
        for (k = ii + 1; k<=ii + len - 1; ++k)  {  //302
            jlu[ju0] = jw[k];
            alu[ju0] = w[k];
            t = t + ABS(w[k]);
            ++ju0;
        }
        
        //     store inverse of diagonal element of u
        //     if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
        if (w[ii]  ==  0.0) w[ii] = tnorm[ii];
        
        alu[ii] = 1.0/w[ii];
        
        //     update pointer to beginning of next row of U.
        jlu[ii + 1] = ju0;
        /*-----------------------------------------------------------------------
         end main loop
         ----------------------------------------------------------------------- */
    }
    
    *nz = ju[n] - 1;
    
    if (cinindex) {
        for(i = 1; i <= *nz; ++i)  --jlu[i];
    }
    
    *ierr = 0;
    
F100:
    ++jw;
    ++ju;
    ++w;
    ++tnorm;
    ++jlu;
    ++alu;
    ++ia;
    ++ja;
    ++a;
    
    fasp_mem_free(ju);
    fasp_mem_free(jw);
    fasp_mem_free(w);
    fasp_mem_free(tnorm);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return;
    
    //     incomprehensible error. Matrix must be wrong.
F995:
    printf("### ERROR: input matrix may be wrong \n");
    *ierr = -1;
    goto F100;
    
    //     insufficient storage in L.
F996:
    printf("### ERROR: insufficient storage in L\n");
    *ierr = -2;
    goto F100;
    
    //     insufficient storage in U.
F997:
    printf("### ERROR: insufficient storage in U\n");
    *ierr = -3;
    goto F100;
    
    //     illegal lfil entered.
F998:
    *ierr = -4;
    printf("### ERROR: illegal lfil entered\n");
    return;
    /*----------------end-of-ilut--------------------------------------------
     -----------------------------------------------------------------------*/
}


/**
 * \fn void fasp_ilutp (INT n, REAL *a, INT *ja, INT *ia, INT lfil, REAL droptol,
 *                      REAL permtol, INT mbloc, REAL *alu, INT *jlu, INT iwk,
 *                      INT *ierr, INT *nz)
 *
 * \brief Get incomplete LU factorization with pivoting dual truncations of
 *        a CSR matrix A
 *
 * \param n   row number of A
 * \param a   nonzero entries of A
 * \param ja  integer array of column for A
 * \param ia  integer array of row pointers for A
 * \param lfil  integer. The fill-in parameter. Each row of L and each row
 *              of U will have a maximum of lfil elements (excluding the diagonal
 *              element). lfil must be .ge. 0.
 * \param droptol  real*8. Sets the threshold for dropping small terms in the
 *                 factorization. See below for details on dropping strategy.
 * \param permtol  tolerance ratio used to  determne whether or not to permute
 *                 two columns.  At step i columns i and j are permuted when
 *                 abs(a(i,j))*permtol .gt. abs(a(i,i))
 *                 [0 --> never permute; good values 0.1 to 0.01]
 * \param mbloc  integer.If desired, permuting can be done only within the diagonal
 *               blocks of size mbloc. Useful for PDE problems with several
 *               degrees of freedom.. If feature not wanted take mbloc=n.
 * \param alu,jlu  matrix stored in Modified Sparse Row (MSR) format containing
 *                 the L and U factors together. The diagonal (stored in
 *                 alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
 *                 contains the i-th row of L (excluding the diagonal entry=1)
 *                 followed by the i-th row of U.
 * \param iwk  integer. The lengths of arrays alu and jlu. If the arrays
 *             are not big enough to store the ILU factorizations, ilut
 *             will stop with an error message.
 * \param ierr  integer pointer. Return error message with the following meaning.
 *                0  --> successful return.
 *               >0  --> zero pivot encountered at step number ierr.
 *               -1  --> Error. input matrix may be wrong.
 *                       (The elimination process has generated a
 *                       row in L or U whose length is .gt.  n.)
 *               -2  --> The matrix L overflows the array al.
 *               -3  --> The matrix U overflows the array alu.
 *               -4  --> Illegal value for lfil.
 *               -5  --> zero row encountered.
 * \param nz  integer pointer. Return number of nonzero entries for alu and jlu
 *
 * \note:  All the diagonal elements of the input matrix must be nonzero.
 *
 * \author Chunsheng Feng
 * \date   09/06/2016
 */
void fasp_ilutp (INT n,
                 REAL *a,
                 INT  *ja,
                 INT  *ia,
                 INT  lfil,
                 REAL droptol,
                 REAL permtol,
                 INT  mbloc,
                 REAL *alu,
                 INT  *jlu,
                 INT  iwk,
                 INT  *ierr,
                 INT  *nz)
{
#if DEBUG_MODE > 0
    printf("### DEBUG: %s (ILUtp) ...... [Start]\n", __FUNCTION__);
#endif
    /*----------------------------------------------------------------------*
     *** ILUTP preconditioner -- ILUT with pivoting  ***              *
     incomplete LU factorization with dual truncation mechanism       *
     ----------------------------------------------------------------------*
     author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996. *
     ----------------------------------------------------------------------*
     on entry:
     ==========
     n       = integer. The dimension of the matrix A.
     
     a,ja,ia = matrix stored in Compressed Sparse Row format.
     ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR
     DETAILS.
     
     lfil    = integer. The fill-in parameter. Each row of L and each row
     of U will have a maximum of lfil elements (excluding the
     diagonal element). lfil must be .ge. 0.
     ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
     EARLIER VERSIONS.
     
     droptol = real*8. Sets the threshold for dropping small terms in the
     factorization. See below for details on dropping strategy.
     
     lfil    = integer. The fill-in parameter. Each row of L and
     each row of U will have a maximum of lfil elements.
     
     permtol = tolerance ratio used to  determne whether or not to permute
     two columns.  At step i columns i and j are permuted when
     
     abs(a(i,j))*permtol .gt. abs(a(i,i))
     
     [0 --> never permute; good values 0.1 to 0.01]
     
     mbloc   = if desired, permuting can be done only within the diagonal
     blocks of size mbloc. Useful for PDE problems with several
     degrees of freedom.. If feature not wanted take mbloc=n.
     
     iwk     = integer. The lengths of arrays alu and jlu. If the arrays
     are not big enough to store the ILU factorizations, ilut
     will stop with an error message.
     
     On return:
     ===========
     
     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
     the L and U factors together. The diagonal (stored in
     alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
     contains the i-th row of L (excluding the diagonal entry=1)
     followed by the i-th row of U.
     
     ju      = integer array of length n containing the pointers to
     the beginning of each row of U in the matrix alu,jlu.
     
     iperm   = contains the permutation arrays.
     iperm(1:n) = old numbers of unknowns
     iperm(n+1:2*n) = reverse permutation = new unknowns.
     
     ierr    = integer. Error message with the following meaning.
     ierr  = 0    --> successful return.
     ierr .gt. 0  --> zero pivot encountered at step number ierr.
     ierr  = -1   --> Error. input matrix may be wrong.
     (The elimination process has generated a
     row in L or U whose length is .gt.  n.)
     ierr  = -2   --> The matrix L overflows the array al.
     ierr  = -3   --> The matrix U overflows the array alu.
     ierr  = -4   --> Illegal value for lfil.
     ierr  = -5   --> zero row encountered.
     
     work arrays:
     =============
     jw      = integer work array of length 2*n.
     w       = real work array of length n
     
     IMPORTANR NOTE:
     --------------
     TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE,
     THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
     changed]. SIMILARLY FOR THE U MATRIX.
     To permute the matrix back to its original state use the loop:
     
     do k=ia(1), ia(n+1)-1
     ja(k) = iperm(ja(k))
     enddo
     
     -----------------------------------------------------------------------*/
    
    //     local variables
    INT k, i, j, jrow, ju0, ii, j1, j2, jpos, len, imax, lenu, lenl, jj, icut,NE;
    REAL s, tmp, tnorm, xmax, xmax0, fact, t;
    SHORT cinindex=0;
    REAL  *w;
    INT  *ju, *jw, *iperm;
    
    if (lfil  <  0) goto F998;
    
    
    ju = (INT *)fasp_mem_calloc(n, sizeof(INT));
    jw = (INT *)fasp_mem_calloc(2*n, sizeof(INT));
    iperm = (INT *)fasp_mem_calloc(2*n, sizeof(INT));
    w = (REAL *)fasp_mem_calloc(n+1, sizeof(REAL));
    
    --ju;
    --jw;
    --iperm;
    --w;
    --jlu;
    --alu;
    --ia;
    --ja;
    --a;
    
    /*-----------------------------------------------------------------------
     shift index for C routines
     -----------------------------------------------------------------------*/
    if (ia[1]  ==  0) cinindex=1 ;
    
    if (cinindex)
    {
        NE = n + 1; //modify by chunsheng 2012, Sep, 1;
        for (i=1; i<=NE; ++i)  ++ia[i];
        NE = ia[n+1] - 1;
        for (i=1; i<=NE; ++i)  ++ja[i];
    }
    
    /*-----------------------------------------------------------------------
     initialize ju0 (points to next element to be added to alu,jlu)
     and pointer array.
     -----------------------------------------------------------------------*/
    ju0 = n + 2;
    jlu[1] = ju0;
    
    
    //	 integer double pointer array.
    for ( j = 1; j <= n; ++j ) { //1
        jw[n + j] = 0;
        iperm[j] = j;
        iperm[n + j] = j;
    } //1
    
    /*-----------------------------------------------------------------------
     beginning of main loop.
     -----------------------------------------------------------------------*/
    for (ii = 1; ii <= n; ++ii ) 	{ //500
        j1 = ia[ii];
        j2 = ia[ii + 1] - 1;
        
        tnorm = 0.0;
        for (k = j1; k <= j2; ++k ) 	tnorm = tnorm + ABS( a[k] ); //501
        if (tnorm  ==  0.0) goto F999;
        tnorm = tnorm/(REAL)(j2 - j1 + 1);
        
        //	 unpack L-part and U-part of row of A in arrays  w  --
        lenu = 1;
        lenl = 0;
        jw[ii] = ii;
        w[ii] = 0.0;
        jw[n + ii] = ii;
        //
        for (j = j1; j <= j2; ++j )  { // 170
            k = iperm[n + ja[j]];
            t = a[j];
            if (k  <  ii) {
                ++lenl;
                jw[lenl] = k;
                w[lenl] = t;
                jw[n + k] = lenl;
            } else if (k  ==  ii) {
                w[ii] = t;
            } else {
                ++lenu;
                jpos = ii + lenu - 1;
                jw[jpos] = k;
                w[jpos] = t;
                jw[n + k] = jpos;
            }
        }  //170
        
        jj = 0;
        len = 0;
        
        
        //	 eliminate previous rows
    F150:
        ++jj;
        if (jj  >  lenl) goto F160;
        
        /*-----------------------------------------------------------------------
         in order to do the elimination in the correct order we must select
         the smallest column index among jw(k), k=jj+1, ..., lenl.
         -----------------------------------------------------------------------*/
        jrow = jw[jj];
        k = jj;
        
        //	 determine smallest column index
        for (j = jj + 1; j <= lenl; ++j) {  //151
            if (jw[j]  <  jrow) {
                jrow = jw[j];
                k = j;
            }
        }
        
        if (k  !=  jj) 	{
            //     exchange in jw
            j = jw[jj];
            jw[jj] = jw[k];
            jw[k] = j;
            //     exchange in jr
            jw[n + jrow] = jj;
            jw[n + j] = k;
            //     exchange in w
            s = w[jj];
            w[jj] = w[k];
            w[k] = s;
        }
        
        //	 zero out element in row by resetting jw(n+jrow) to zero.
        jw[n + jrow] = 0;
        
        //	 get the multiplier for row to be eliminated: jrow
        fact = w[jj]*alu[jrow];
        
        //	 drop term if small
        if (ABS(fact)  <=  droptol) goto F150;
        
        //	 combine current row and row jrow
        
        for ( k = ju[jrow]; k <= jlu[jrow + 1] - 1; ++k ) {  //203
            s = fact*alu[k];
            // new column number
            j = iperm[n + jlu[k]];
            jpos = jw[n + j];
            if (j  >=  ii) {
                //	 dealing with upper part.
                if (jpos  ==  0) {
                    //	 this is a fill-in element
                    ++lenu;
                    i = ii + lenu - 1;
                    if (lenu  >  n) goto F995;
                    jw[i] = j;
                    jw[n + j] = i;
                    w[i] = -s;
                } else {
                    //     no fill-in element --
                    w[jpos] = w[jpos] - s;
                }
                
            } else {
                //	 dealing with lower part.
                if (jpos  ==  0) {
                    //	 this is a fill-in element
                    ++lenl;
                    if (lenl  >  n) goto F995;
                    jw[lenl] = j;
                    jw[n + j] = lenl;
                    w[lenl] = -s;
                } else {
                    //	 this is not a fill-in element
                    w[jpos] = w[jpos] - s;
                }
            }
        }  //203
        
        /*
         store this pivot element -- (from left to right -- no danger of
         overlap with the working elements in L (pivots).
         */
        
        ++len;
        w[len] = fact;
        jw[len] = jrow;
        goto F150;
        
    F160:
        //	 reset double-pointer to zero (U-part)
        for ( k = 1; k <= lenu; ++k ) jw[n + jw[ii + k - 1]] = 0;  //308
        
        //	 update L-matrix
        lenl = len;
        len = MIN(lenl, lfil);
        
        //	 sort by quick-split
        fasp_qsplit(&w[1], &jw[1], lenl, len);
        
        //	 store L-part -- in original coordinates ..
        for ( k = 1; k <= len; ++k ) {  // 204
            if (ju0  >  iwk) goto F996;
            alu[ju0] = w[k];
            jlu[ju0] = iperm[jw[k]];
            ++ju0;
        }  //204
        
        //	 save pointer to beginning of row ii of U
        ju[ii] = ju0;
        
        //	 update U-matrix -- first apply dropping strategy
        len = 0;
        for(k = 1; k <= lenu - 1; ++k ) {
            if ( ABS(w[ii + k])  >  droptol*tnorm) {
                ++len;
                w[ii + len] = w[ii + k];
                jw[ii + len] = jw[ii + k];
            }
        }
        
        lenu = len + 1;
        len = MIN(lenu, lfil);
        fasp_qsplit(&w[ii + 1], &jw[ii + 1], lenu-1, len);
        
        //	 determine next pivot --
        imax = ii;
        xmax = ABS(w[imax]);
        xmax0 = xmax;
        icut = ii - 1 + mbloc - (ii - 1)%mbloc;
        
        for ( k = ii + 1; k <= ii + len - 1; ++k ) {
            t = ABS(w[k]);
            if ((t  >  xmax) && (t*permtol  >  xmax0) && (jw[k]  <=  icut)) {
                imax = k;
                xmax = t;
            }
        }
        
        //	 exchange w's
        tmp = w[ii];
        w[ii] = w[imax];
        w[imax] = tmp;
        
        //	 update iperm and reverse iperm
        j = jw[imax];
        i = iperm[ii];
        iperm[ii] = iperm[j];
        iperm[j] = i;
        
        //	 reverse iperm
        iperm[n + iperm[ii]] = ii;
        iperm[n + iperm[j]] = j;
        
        //-----------------------------------------------------------------------
        if (len + ju0  >  iwk) goto F997;
        
        
        //	 copy U-part in original coordinates
        for ( k = ii + 1; k <= ii + len - 1; ++k ) { //302
            jlu[ju0] = iperm[jw[k]];
            alu[ju0] = w[k];
            ++ju0;
        }
        
        //	 store inverse of diagonal element of u
        if (w[ii]  ==  0.0) w[ii] = (1.0e-4 + droptol)*tnorm;
        alu[ii] = 1.0/w[ii];
        
        //	 update pointer to beginning of next row of U.
        jlu[ii + 1] = ju0;
        
        /*-----------------------------------------------------------------------
         end main loop
         -----------------------------------------------------------------------*/
    }  //500
    
    //	 permute all column indices of LU ...
    for ( k = jlu[1]; k <= jlu[n + 1] - 1; ++k ) 	jlu[k] = iperm[n + jlu[k]];
    
    //	 ...and of A
    for ( k = ia[1]; k <= ia[n + 1] - 1; ++k )	ja[k] = iperm[n + ja[k]];
    
    *nz = ju[n]- 1;
    
    if (cinindex)  {
        for (i = 1; i <= *nz; ++i ) --jlu[i];
    }
    
    *ierr = 0;
    
F100:
    ++jw;
    ++ju;
    ++iperm;
    ++w;
    ++jlu;
    ++alu;
    ++ia;
    ++ja;
    ++a;
    
    fasp_mem_free(ju);
    fasp_mem_free(jw);
    fasp_mem_free(iperm);
    fasp_mem_free(w);
    
#if DEBUG_MODE > 0
    printf("### DEBUG: %s ...... [Finish]\n", __FUNCTION__);
#endif

    return;
    
    //	 incomprehensible error. Matrix must be wrong.
F995:
    printf("### ERROR: input matrix may be wrong \n");
    *ierr = -1;
    goto F100;
    
    //	 insufficient storage in L.
F996:
    printf("### ERROR: insufficient storage in L\n");
    *ierr = -2;
    goto F100;
    
    //	 insufficient storage in U.
F997:
    printf("### ERROR: insufficient storage in U\n");
    *ierr = -3;
    goto F100;
    
    //	 illegal lfil entered.
F998:
    printf("### ERROR: illegal lfil entered\n");
    *ierr = -4;
    return;
    
    
    //	 zero row encountered
F999:
    printf("### ERROR: zero row encountered\n");
    *ierr = -5;
    goto F100;
    //----------------end-of-ilutp-------------------------------------------
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
