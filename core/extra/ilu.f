c     
c     ilu.f
c     
c-----------------------------------------------------------------------
c     Created by Shiquan Zhang on 12/27/09.
c     Modified by Chensong Zhang on 05/25/10.
c-----------------------------------------------------------------------

!> \file ilu.f
!> \brief ILU routines for preconditioning adapted from SparseKit
!>
!> Methods: ILUt, ILUtp, ILUk, ILUs

c-----------------------------------------------------------------------
      subroutine iluk(n,a,ja,ia,lfil,alu,jlu,iwk,ierr,nzlu)
c-----------------------------------------------------------------------
      implicit none 
      integer n,nzlu,lfil,iwk,ierr
      real*8 a(*),alu(*),w(n)
      integer ja(*),ia(n+1),jlu(*),ju(n),levs(iwk),jw(3*n)
c----------------------------------------------------------------------
c     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k)) 
c----------------------------------------------------------------------
c     
c     on entry:
c==========
c     n       = integer. The row dimension of the matrix A. The matrix 
c     
c     a,ja,ia = matrix stored in Compressed Sparse Row format.              
c     
c     lfil    = integer. The fill-in parameter. Each element whose
c     leve-of-fill exceeds lfil during the ILU process is dropped.
c     lfil must be .ge. 0 
c     
c     iwk     = integer. The minimum length of arrays alu, jlu, and levs.
c     
c     On return:
c===========
c     
c     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c     the L and U factors together. The diagonal (stored in
c     alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c     contains the i-th row of L (excluding the diagonal entry=1)
c     followed by the i-th row of U.
c     
c     ju      = integer array of length n containing the pointers to
c     the beginning of each row of U in the matrix alu,jlu.
c     
c     levs    = integer (work) array of size iwk -- which contains the 
c     levels of each element in alu, jlu.
c     
c     ierr    = integer. Error message with the following meaning.
c     ierr  = 0    --> successful return.
c     ierr .gt. 0  --> zero pivot encountered at step number ierr.
c     ierr  = -1   --> Error. input matrix may be wrong.
c     (The elimination process has generated a
c     row in L or U whose length is .gt.  n.)
c     ierr  = -2   --> The matrix L overflows the array al.
c     ierr  = -3   --> The matrix U overflows the array alu.
c     ierr  = -4   --> Illegal value for lfil.
c     ierr  = -5   --> zero row encountered in A or U.
c     
c     work arrays:
c=============
c     jw      = integer work array of length 3*n.
c     w       = real work array of length n 
c     
c----------------------------------------------------------------------
c     w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = U-part] 
c     jw(n+1:2n)  stores the nonzero indicator. 
c     
c     Notes:
c     ------
c     All the diagonal elements of the input matrix must be nonzero.
c     
c----------------------------------------------------------------------* 
c     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,n2,
     *     jlev, min 
      real*8 t, s, fact
      logical cinindex

      if (lfil .lt. 0) goto 998

c-----------------------------------------------------------------------
c     shift index for C routines
c-----------------------------------------------------------------------
      cinindex = ( ia(1) .eq. 0 )

      if (cinindex) then 
         do i=1,n+1
            ia(i)=ia(i)+1
         end do
         do i = 1,ia(n+1)-1
            ja(i)=ja(i)+1
         end do
      end if

c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      n2 = n+n 
      ju0 = n+2
      jlu(1) = ju0
c     
c     initialize nonzero indicator array + levs array -- 
c     
      do 1 j=1,2*n 
         jw(j) = 0
 1    continue
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
c     
c     unpack L-part and U-part of row of A in arrays w 
c     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c     
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (t .eq. 0.0) goto 170 
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n2+lenl) = 0 
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
               jw(n2+ii) = 0 
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n2+jpos) = 0 
               jw(n+k) = jpos
            endif
 170     continue
c     
         jj = 0
c     
c     eliminate previous rows
c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c     
c     determine smallest column index
c     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c     
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jw(n+  (pointers/ nonzero indicator).
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in jw(n2+  (levels) 
            j = jw(n2+jj) 
            jw(n2+jj)  = jw(n2+k) 
            jw(n2+k) = j
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c     
c     zero out element in row by resetting jw(n+jrow) to zero.
c     
         jw(n+jrow) = 0
c     
c     get the multiplier for row to be eliminated (jrow) + its level
c     
         fact = w(jj)*alu(jrow)
         jlev = jw(n2+jj) 
         if (jlev .gt. lfil) goto 150
c     
c     combine current row and row jrow
c     
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
c     
c     dealing with upper part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
                  jw(n2+i) = jlev+levs(k)+1 
               else
c     
c     this is not a fill-in element 
c     
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            else
c     
c     dealing with lower part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
                  jw(n2+lenl) = jlev+levs(k)+1 
               else
c     
c     this is not a fill-in element 
c     
                  w(jpos) = w(jpos) - s
                  jw(n2+jpos) = min(jw(n2+jpos),jlev+levs(k)+1)
               endif
            endif
 203     continue
         w(jj) = fact
         jw(jj)  = jrow
         goto 150 
 160     continue 
c     
c     reset double-pointer to zero (U-part) 
c     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
c     
c     update l-matrix
c     
         do 204 k=1, lenl 
            if (ju0 .gt. iwk) goto 996
            if (jw(n2+k) .le. lfil) then
               alu(ju0) =  w(k)
               jlu(ju0) =  jw(k)
               ju0 = ju0+1
            endif
 204     continue
c     
c     save pointer to beginning of row ii of U
c     
         ju(ii) = ju0
c     
c     update u-matrix
c     
         do 302 k=ii+1,ii+lenu-1 
            if (ju0 .gt. iwk) goto 997
            if (jw(n2+k) .le. lfil) then
               jlu(ju0) = jw(k)
               alu(ju0) = w(k)
               levs(ju0) = jw(n2+k) 
               ju0 = ju0+1
            endif
 302     continue

         if (w(ii) .eq. 0.0) goto 999 
c     
         alu(ii) = 1.0d0/ w(ii) 
c     
c     update pointer to beginning of next row of U.
c     
         jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue

      nzlu = ju(n)-1

      if (cinindex) then
         do i=1,nzlu
            jlu(i)=jlu(i)-1
         end do
      end if 

      ierr = 0
      return
c     
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      return
c     
c     zero row encountered in A or U. 
c     
 999  ierr = -5
      return
c----------------end-of-iluk--------------------------------------------
c-----------------------------------------------------------------------
      end
c     
c-----------------------------------------------------------------------
      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,iwk,ierr,nz)
c-----------------------------------------------------------------------
      implicit none 
      integer n 
      real*8 a(1),alu(1),w(n+1),droptol
      integer ja(1),ia(1),jlu(1),ju(n),jw(2*n),lfil,iwk,ierr,nz
c----------------------------------------------------------------------*
c     *** ILUT preconditioner ***                                      *
c     incomplete LU factorization with dual truncation mechanism       *
c----------------------------------------------------------------------*
c     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
c----------------------------------------------------------------------*
c     PARAMETERS                                                           
c-----------
c     
c     on entry:
c==========
c     n       = integer. The row dimension of the matrix A. The matrix 
c     
c     a,ja,ia = matrix stored in Compressed Sparse Row format.              
c     
c     lfil    = integer. The fill-in parameter. Each row of L and each row
c     of U will have a maximum of lfil elements (excluding the diagonal
c     element). lfil must be .ge. 0.
c     
c     droptol = real*8. Sets the threshold for dropping small terms in the
c     factorization. See below for details on dropping strategy.
c     
c     iwk     = integer. The lengths of arrays alu and jlu. If the arrays
c     are not big enough to store the ILU factorizations, ilut
c     will stop with an error message. 
c     
c     On return:
c===========
c     
c     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c     the L and U factors together. The diagonal (stored in
c     alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c     contains the i-th row of L (excluding the diagonal entry=1)
c     followed by the i-th row of U.
c     
c     ju      = integer array of length n containing the pointers to
c     the beginning of each row of U in the matrix alu,jlu.
c     
c     ierr    = integer. Error message with the following meaning.
c     ierr  = 0    --> successful return.
c     ierr .gt. 0  --> zero pivot encountered at step number ierr.
c     ierr  = -1   --> Error. input matrix may be wrong.
c     (The elimination process has generated a
c     row in L or U whose length is .gt.  n.)
c     ierr  = -2   --> The matrix L overflows the array al.
c     ierr  = -3   --> The matrix U overflows the array alu.
c     ierr  = -4   --> Illegal value for lfil.
c     ierr  = -5   --> zero row encountered.
c     
c     work arrays:
c=============
c     jw      = integer work array of length 2*n.
c     w       = real work array of length n+1.
c     
c----------------------------------------------------------------------
c     w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
c     jw(n+1:2n)  stores nonzero indicators
c     
c     Notes:
c     ------
c     The diagonal elements of the input matrix must be nonzero (at least
c     'structurally'). 
c     
c----------------------------------------------------------------------* 
c---- Dual drop strategy works as follows.                             *
c     *
c     1) Theresholding in L and U as set by droptol. Any element whose *
c     magnitude is less than some tolerance (relative to the abs       *
c     value of diagonal element in u) is dropped.                      *
c     *
c     2) Keeping only the largest lfil elements in the i-th row of L   * 
c     and the largest lfil elements in the i-th row of U (excluding    *
c     diagonal elements).                                              *
c     *
c     Flexibility: one can use droptol=0 to get a strategy based on    *
c     keeping the largest elements in each row of L and U. Taking      *
c     droptol .ne. 0 but lfil=n will give the usual threshold strategy *
c     (however, fill-in is then unpredictible).                        *
c----------------------------------------------------------------------*
c     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      real*8 tnorm, t, abs, s, fact 
      logical cinindex

      if (lfil .lt. 0) goto 998

      cinindex = ia(1) .eq. 0

      if (cinindex) then
         do i=1,n+1
            ia(i)=ia(i)+1
         end do
         do i = 1,ia(n+1)-1
            ja(i)=ja(i)+1
         end do
      end if

c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
c     
c     initialize nonzero indicator array. 
c     
      do 1 j=1,n
         jw(n+j) = 0
 1    continue
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
c     
c     unpack L-part and U-part of row of A in arrays w 
c     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c     
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
 170     continue
         jj = 0
         len = 0 
c     
c     eliminate previous rows
c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c     
c     determine smallest column index
c     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c     
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c     
c     zero out element in row by setting jw(n+jrow) to zero.
c     
         jw(n+jrow) = 0
c     
c     get the multiplier for row to be eliminated (jrow).
c     
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
c     
c     combine current row and row jrow
c     
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
c     
c     dealing with upper part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
c     
c     this is not a fill-in element 
c     
                  w(jpos) = w(jpos) - s

               endif
            else
c     
c     dealing  with lower part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
c     
c     this is not a fill-in element 
c     
                  w(jpos) = w(jpos) - s
               endif
            endif
 203     continue
c     
c     store this pivot element -- (from left to right -- no danger of
c     overlap with the working elements in L (pivots). 
c     
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150
 160     continue
c     
c     reset double-pointer to zero (U-part)
c     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
c     
c     update L-matrix
c     
         lenl = len 
         len = min0(lenl,lfil)
c     
c     sort by quick-split
c     
         call qsplit (w,jw,lenl,len)
c     
c     store L-part
c     
         do 204 k=1, len 
            if (ju0 .gt. iwk) goto 996
            alu(ju0) = w(k)
            jlu(ju0) = jw(k)
            ju0 = ju0+1
 204     continue
c     
c     save pointer to beginning of row ii of U
c     
         ju(ii) = ju0
c     
c     update U-matrix -- first apply dropping strategy 
c     
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
c     
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
c     
c     copy
c     
         t = abs(w(ii))
         if (len + ju0 .gt. iwk) goto 997
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
 302     continue
c     
c     store inverse of diagonal element of u
c     
         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
c     
         alu(ii) = 1.0d0 / w(ii) 
c     
c     update pointer to beginning of next row of U.
c     
         jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue

      nz=ju(n)-1

      if (cinindex) then 
         do i=1,nz
            jlu(i)=jlu(i)-1
         end do
      end if

      ierr = 0
      return
c     
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      write(*,*) 'input matrix may be wrong'
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      write(*,*),'insufficient storage in L'
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      write(*,*),'insufficient storage in L'
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      write(*,*),'illegal lfil entered'
      return
c     
c     zero row encountered
c     
 999  ierr = -5
      write(*,*),'zero row encountered'
      return
c----------------end-of-ilut--------------------------------------------
c-----------------------------------------------------------------------
      end
c     
c----------------------------------------------------------------------
      subroutine ilutp(n,a,ja,ia,lfil,droptol,permtol,mbloc,alu,
     *     jlu,iwk,ierr,nz)
c-----------------------------------------------------------------------
c     implicit none
      integer n,ja(1),ia(1),lfil,jlu(1),ju(n),jw(2*n),iwk,
     *     iperm(2*n),ierr
      real*8 a(1), alu(1), w(n+1), droptol
c----------------------------------------------------------------------*
c     *** ILUTP preconditioner -- ILUT with pivoting  ***              *
c     incomplete LU factorization with dual truncation mechanism       *
c----------------------------------------------------------------------*
c     author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996. *
c----------------------------------------------------------------------*
c     on entry:
c==========
c     n       = integer. The dimension of the matrix A.
c     
c     a,ja,ia = matrix stored in Compressed Sparse Row format.
c     ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR 
c     DETAILS. 
c     
c     lfil    = integer. The fill-in parameter. Each row of L and each row
c     of U will have a maximum of lfil elements (excluding the 
c     diagonal element). lfil must be .ge. 0.
c     ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
c     EARLIER VERSIONS. 
c     
c     droptol = real*8. Sets the threshold for dropping small terms in the
c     factorization. See below for details on dropping strategy.
c     
c     lfil    = integer. The fill-in parameter. Each row of L and
c     each row of U will have a maximum of lfil elements.
c     
c     permtol = tolerance ratio used to  determne whether or not to permute
c     two columns.  At step i columns i and j are permuted when 
c     
c     abs(a(i,j))*permtol .gt. abs(a(i,i))
c     
c     [0 --> never permute; good values 0.1 to 0.01]
c     
c     mbloc   = if desired, permuting can be done only within the diagonal
c     blocks of size mbloc. Useful for PDE problems with several
c     degrees of freedom.. If feature not wanted take mbloc=n.
c          
c     iwk     = integer. The lengths of arrays alu and jlu. If the arrays
c     are not big enough to store the ILU factorizations, ilut
c     will stop with an error message. 
c     
c     On return:
c===========
c     
c     alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c     the L and U factors together. The diagonal (stored in
c     alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c     contains the i-th row of L (excluding the diagonal entry=1)
c     followed by the i-th row of U.
c     
c     ju      = integer array of length n containing the pointers to
c     the beginning of each row of U in the matrix alu,jlu.
c     
c     iperm   = contains the permutation arrays. 
c     iperm(1:n) = old numbers of unknowns
c     iperm(n+1:2*n) = reverse permutation = new unknowns.
c     
c     ierr    = integer. Error message with the following meaning.
c     ierr  = 0    --> successful return.
c     ierr .gt. 0  --> zero pivot encountered at step number ierr.
c     ierr  = -1   --> Error. input matrix may be wrong.
c     (The elimination process has generated a
c     row in L or U whose length is .gt.  n.)
c     ierr  = -2   --> The matrix L overflows the array al.
c     ierr  = -3   --> The matrix U overflows the array alu.
c     ierr  = -4   --> Illegal value for lfil.
c     ierr  = -5   --> zero row encountered.
c     
c     work arrays:
c=============
c     jw      = integer work array of length 2*n.
c     w       = real work array of length n 
c     
c     IMPORTANR NOTE:
c     --------------
c     TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE, 
C     THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
c     changed]. SIMILARLY FOR THE U MATRIX. 
c     To permute the matrix back to its original state use the loop:
c     
c     do k=ia(1), ia(n+1)-1
c        ja(k) = iperm(ja(k)) 
c     enddo
c     
c-----------------------------------------------------------------------
c     local variables
c     
      integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,mbloc,
     *     icut
      real*8 s, tmp, tnorm,xmax,xmax0, fact, abs, t, permtol
      logical cinindex
      
      if (lfil .lt. 0) goto 998

c-----------------------------------------------------------------------
c     shift index for C routines
c-----------------------------------------------------------------------
      cinindex = ia(1) .eq. 0

      if (cinindex) then
         do i=1,n+1
            ia(i)=ia(i)+1
         end do
         do i = 1,ia(n+1)-1
            ja(i)=ja(i)+1
         end do
      end if

c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
c     
c     integer double pointer array.
c     
      do 1 j=1, n
         jw(n+j)  = 0
         iperm(j) = j
         iperm(n+j) = j
 1    continue
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1

         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/(j2-j1+1)
c     
c     unpack L-part and U-part of row of A in arrays  w  --
c     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c     
         do 170  j = j1, j2
            k = iperm(n+ja(j))
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
 170     continue
         jj = 0
         len = 0 
c     
c     eliminate previous rows
c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c     
c     determine smallest column index
c     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c     
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c     
c     zero out element in row by resetting jw(n+jrow) to zero.
c     
         jw(n+jrow) = 0
c     
c     get the multiplier for row to be eliminated: jrow
c     
         fact = w(jj)*alu(jrow)
c     
c     drop term if small
c     
         if (abs(fact) .le. droptol) goto 150
c     
c     combine current row and row jrow
c     
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
c     new column number
            j = iperm(n+jlu(k))
            jpos = jw(n+j)
            if (j .ge. ii) then
c     
c     dealing with upper part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenu = lenu+1
                  i = ii+lenu-1 
                  if (lenu .gt. n) goto 995
                  jw(i) = j
                  jw(n+j) = i 
                  w(i) = - s
               else
c     no fill-in element --
                  w(jpos) = w(jpos) - s
               endif
            else
c     
c     dealing with lower part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
c     
c     this is not a fill-in element
c     
                  w(jpos) = w(jpos) - s
               endif
            endif
 203     continue
c     
c     store this pivot element -- (from left to right -- no danger of
c     overlap with the working elements in L (pivots). 
c     
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150
 160     continue
c     
c     reset double-pointer to zero (U-part)
c     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
c     
c     update L-matrix
c     
         lenl = len 
         len = min0(lenl,lfil)
c     
c     sort by quick-split
c     
         call qsplit (w,jw,lenl,len)
c     
c     store L-part -- in original coordinates ..
c     
         do 204 k=1, len
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)  
            jlu(ju0) = iperm(jw(k))
            ju0 = ju0+1
 204     continue
c     
c     save pointer to beginning of row ii of U
c     
         ju(ii) = ju0
c     
c     update U-matrix -- first apply dropping strategy 
c     
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
c     
c     determine next pivot -- 
c     
         imax = ii
         xmax = abs(w(imax))
         xmax0 = xmax
         icut = ii - 1 + mbloc - mod(ii-1,mbloc)
         do k=ii+1,ii+len-1
            t = abs(w(k))
            if (t .gt. xmax .and. t*permtol .gt. xmax0 .and.
     *           jw(k) .le. icut) then
               imax = k
               xmax = t
            endif
         enddo
c     
c     exchange w's
c     
         tmp = w(ii)
         w(ii) = w(imax)
         w(imax) = tmp
c     
c     update iperm and reverse iperm
c     
         j = jw(imax)
         i = iperm(ii)
         iperm(ii) = iperm(j)
         iperm(j) = i
c     
c     reverse iperm
c     
         iperm(n+iperm(ii)) = ii
         iperm(n+iperm(j)) = j
c-----------------------------------------------------------------------
c     
         if (len + ju0 .gt. iwk) goto 997
c     
c     copy U-part in original coordinates
c     
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = iperm(jw(k))
            alu(ju0) = w(k)
            ju0 = ju0+1
 302     continue
c     
c     store inverse of diagonal element of u
c     
         if (w(ii) .eq. 0.0) w(ii) = (1.0D-4 + droptol)*tnorm
         alu(ii) = 1.0d0/ w(ii) 
c     
c     update pointer to beginning of next row of U.
c     
         jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue
c     
c     permute all column indices of LU ...
c     
      do k = jlu(1),jlu(n+1)-1
         jlu(k) = iperm(n+jlu(k))
      enddo
c     
c     ...and of A
c     
      do k=ia(1), ia(n+1)-1
         ja(k) = iperm(n+ja(k))
      enddo
c     
      nz=ju(n)-1

      if (cinindex) then
         do i=1,nz
            jlu(i)=jlu(i)-1
         end do
      end if

      ierr = 0
      return
c     
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      return
c     
c     zero row encountered
c     
 999  ierr = -5
      return
c----------------end-of-ilutp-------------------------------------------
      end
c     
c-----------------------------------------------------------------------
      subroutine ilus(n,a,colind,rwptr,
     1     s,relax,nzlu,luval,ijlu,maxstr,ierr)
c-----------------------------------------------------------------------
c     *
c     Incomplete LU factorization for s levels of fill-in.	       	*
c     *
c========================================================================
c     This routine performs incomplete LU factorization, with s levels	*
c     of fill-in. It is passed the matrix A in CSR format, and a real*8 *
c     parameter relax, which will indicate whether or not to perform	*
c     MILU factorization. If the parameter is zero, regular ILU is	*
c     performed instead. A single integer work array is passed in which *
c     to create the LU factors data structures, and for intermediate	*
c     storage in the computations. What is returned are three pointers	*
c     to the start locations in the work array, for the combined MSR	*
c     data structure containing the LU factors. These pointers will be 	*
c     aligned on double word boundaries, provided that the initial word	*
c     in the work array is aligned on a double word boundary.		*
c     *
c     If not enough work space is provided, an error code is returned	*
c     and the pointers are not set. An estimate of the additional       *
c     work space needed is returned in mneed, but this is only a guess.	*
c========================================================================
c     
c     parameters
c-----------
c     
c     on entry:
c==========
c     n       = integer. the dimension of the matrix a.
c     
c     a,colind,rwptr = matrix stored in CSR format.
c     
c     s       = integer; the fill-in parameter. Fill-in generated during the
c     factorization is kept only if it is caused by a nonzero element
c     in the matrix that was created on the s level or less.  Fill
c     generated by the original nonzeros are of level 1, those created
c     by such entries are level 2, etc. The original matrix entries
c     are defined as being of level 0. See symfac() for details.
c     
c     relax   = relaxation parameter. relax = 0 gives ilu(s), relax = 1
c     gives milu(s), and values between give the multiplier to use
c     before adding discarded fill to the diagonal.
c     
c     maxstr  = integer giving amount of integer word space available.
c     
c     intermediate variables:
c========================
c     
c     rowll   = pointer to temp vector of levels in a row
c     lastcol = pointer to temp vector of last column to update a row
c     levels  = pointer to temp vector keeping all levels during symfac.
c     colptrs = pointer to temp vector keeping column pointers in numfac.
c     
c     on return:
c===========
c     
c     lu, jlu = (pointers) matrix stored in modified sparse row (msr) 
c     format containing the l and u factors together. the diagonal (stored 
c     in lu(1:n) ) is inverted. each i-th row of the lu, jlu matrix 
c     contains the i-th row of l (excluding the diagonal entry=1) 
c     followed by the i-th row of u.
c     
c     ilup    = (pointer) integer array of length n containing the pointers to        
c     the beginning of each row of u in the matrix lu, jlu. 
c     
c     on return:
c===========
c     luval, ijlu = matrix stored in modified sparse row (msr) format containing
c     the l and u factors together. the diagonal (stored in
c     lu(1:n) ) is inverted. each i-th row of the luval, ijlu matrix 
c     contains the i-th row of l (excluding the diagonal entry=1) 
c     followed by the i-th row of u.  
c     ierr is an error flag:
c     ierr  = -i --> near zero pivot in step i
c     ierr  = 0  --> all's OK
c     ierr  = 1  --> not enough storage; check mneed.
c     ierr  = 2  --> illegal parameter
c     
c     mneed   = contains number of additional storage locations needed for lu
c     in case of failure, or total integer words used from work array
c     in case of success.
c     
c     Storage requirements:
c======================
c     
c     jlu:	nzlu  	          integers
c     ilup:	n 	   	  integers
c     lu:       nzlu 		  reals
c     luval     nzlu              reals
c     ijlu      nzlu              intergers
c     
c     temporary arrays
c     ----------------
c     rowll       n               integers
c     lastcol     n               integers
c     levels      nzlu            integers
c     colptrs     n               integers
c     
c     Warnings:
c===========
c     
c     1) A must have nonzero diagonal entries, at least. This does not assure
c     completion of the factorization, however.
c     2) The pointers ilup, jlu, lu are given as indices to the work array
c     that is passed in. So if the integer work array is passed in as
c     work(ptr), you will need to offset those pointers by the operations
c        jlu  = jlu  + ptr -1
c        ilup = ilup + ptr -1
c        lu   = lu   + ptr -1
c     on return.
c     
c========================================================================
      implicit none
c     
      integer      intdbl
      parameter  ( intdbl = 2 )
c     
      real*8  a(1), relax, luval(1)
      integer n, s, ierr, rwptr(1), colind(1), ijlu(1), i
      integer mneed, maxstr, nzlu, remain
c
c     Using iwork as a temp work space (double and int mixed). BAD!!! --zcs
c
      integer iwork(maxstr) 

      logical milu, cinindex

c     ----------------------
c     Pointers in work array
c     ----------------------
      integer lu, jlu, ilup
      integer lastcol, rowll, levels, colptrs

      ierr  =  0

c     ------------------------------------------------
c     transfer c index to fortran
c     ------------------------------------------------
      cinindex = rwptr(1) .eq. 0

      if (cinindex) then
         do i=1,n+1
            rwptr(i)=rwptr(i)+1
         end do
         do i = 1,rwptr(n+1)-1
            colind(i)=colind(i)+1
         end do
      end if

c     ------------------------------------------------
c     If requested level of fill is negative, set to 0
c     ------------------------------------------------
      s = max(s,0)

c     ----------------------------------------------
c     If relaxation parameter is zero, method is ILU
c     ----------------------------------------------
      milu = (relax .ne. 0.0d0)

c     ----------------------------------------------------------------
c     Compute pointers into work array.  This version uses two arrays
c     dependent on the (unknown in advance) size nzlu, jlu(1:nzlu) and
c     levels(1:nzlu).  To handle this, first the data space for the
c     work arrays in symfac is allocated, then half the remainder is
c     assigned to jlu, the rest to levels.  The partitioning of the 
c     integer work array is as:
c     
c     [ilup   jlu      ...      levels    rowll lastcol],
c     [-n-----nzlu-----...-------nzlu-------n------n---]
c     
c     The index of the last entry in the work array is maxstr, and
c     rowll and lastcol start at the far end of the work array.
c     Note that these are aligned on double word boundaries.
c     ----------------------------------------------------------------
      lastcol = maxstr - n + 1
      rowll = lastcol - n

c     ----------------------------------------------------------
c     Set ilup to the far left of the avaiable work array space.
c     ----------------------------------------------------------
      ilup = 1
      jlu  = ilup + n + mod(n,intdbl)

c     ----------------------------------------------------------------
c     Of space remaining, can allocate half to jlu and half to levels.
c     Compute the halfway dividing mark "remain".  Note that
c     growth available = total - lastcol  - rowll  - ilup 
c     ----------------------------------------------------------------
      remain = (maxstr - 3*n) / 2

c     -----------------------------------------------------
c     Now levels = storage for ilup + half remaining words.
c     -----------------------------------------------------
      levels = n + remain

c     ------------------------------------------------------------------
c     If remain is nonpositive, estimate 4n as needed storage, 2n for LU
c     and 2n for jlu: this will handle a diagonal preconditioner.
c     ------------------------------------------------------------------
      if (remain .le. 0) then
         ierr = 1
         mneed = 4*n + 1
         return
      end if
c     
c     -------------------------------
c     Perform symbolic factorization:
c     -------------------------------
      call symfac(n  , colind    ,rwptr     ,s      ,remain ,nzlu ,
     1     iwork(jlu) , iwork(ilup), iwork(lastcol), iwork(levels),
     2     iwork(rowll), ierr , mneed)

c     -----------------------------------------------------------------
c     If error encountered (probably because of lack of storage space),
c     return.
c     -----------------------------------------------------------------
      if (ierr .ne. 0) then
         return
      endif

c     -----------------------------------------------------------------
c     Check to see if enough space is available for the lu factorization.
c     The partitioning of the integer work array is as:
c     
c     [ilup   jlu           LU        colptrs],
c     [-n-----nzlu-----intdbl*nzlu-------n---]
c     
c     -----------------------------------------------------------------
      colptrs = maxstr - n + 1
      remain = maxstr - n - nzlu - n - 2 ! for safe, -2, c.s.
      if (remain .lt. intdbl*nzlu) then
         ierr = 1
         mneed = intdbl*nzlu - remain
         write(*,901) 'ierr=',ierr
         write(*,901) 'Error: need more memory',mneed
         return
      end if

c     ----------------------------------------
c     Set pointer to LU in integer work array:
c     ----------------------------------------
      lu = jlu + nzlu
      lu = lu + mod(lu-1,intdbl) + 2 ! for safe, +2, c.s.

c     --------------------------------
c     Perform numerical factorization.
c     --------------------------------
      call numfac(n, colind, rwptr, iwork(jlu), iwork(ilup), a, 
     &     iwork(lu), relax, ierr, iwork(colptrs), milu, nzlu, 
     &     ijlu, luval)

c     ----------------------------------------------------
c     Return actual integer words used in LU factors.
c     The last two terms are from forcing alignment of the
c     arrays jlu and alu on double word boundaries.
c     ----------------------------------------------------
      mneed = (1 + intdbl)*nzlu + n + mod(n,intdbl)
     &     + mod(nzlu,intdbl)
c     
c     transfer fortran index back to c
      if (cinindex) then 
         do i=1,nzlu
            ijlu(i)=ijlu(i)-1
         end do
      end if

 901  format(A,I15)

      return

c====================End of ILUS =======================================
      end
c     
c========================================================================
      subroutine symfac(n,colind,rwptr,levfill,nzmax,nzlu,
     1     ijlu,uptr,lastcol,levels,rowll,
     2     ierr,mneed)
c========================================================================
      implicit none
c========================================================================
c     *
c     Symbolic factorization of a matrix in compressed sparse row format*
c     with resulting factors stored in a single MSR data structure.     *
c     *
c     This routine uses the CSR data structure of A in two integer vectors*
c     colind, rwptr to set up the data structure for the ILU(levfill) 	*
c     factorization of A in the integer vectors ijlu and uptr.  Both L	*
c     and U are stored in the same structure, and uptr(i) is the pointer*
c     to the beginning of the i-th row of U in ijlu.			*
c     *
c========================================================================
c     *
c     Method Used                                                       *
c     ===========                                                      	*
c     *
c     The implementation assumes that the diagonal entries are		*
c     nonzero, and remain nonzero throughout the elimination		*
c     process.  The algorithm proceeds row by row.  When computing	*
c     the sparsity pattern of the i-th row, the effect of row		*
c     operations from previous rows is considered.  Only those		*
c     preceding rows j for which (i,j) is nonzero need be considered,	*
c     since otherwise we would not have formed a linear combination	*
c     of rows i and j.							*
c     *
c     The method used has some variations possible.  The definition	*
c     of ILU(s) is not well specified enough to get a factorization	*
c     that is uniquely defined, even in the sparsity pattern that	*
c     results.  For s = 0 or 1, there is not much variation, but for	*
c     higher levels of fill the problem is as follows:  Suppose		*
c     during the decomposition while computing the nonzero pattern	*
c     for row i the following principal submatrix is obtained:		*
c     _______________________						*
c     |          |           |					*
c     |          |           |					*
c     |  j,j     |    j,k    |					*
c     |          |           |					*
c     |__________|___________|					*
c     |          |           |					*
c     |          |           |					*
c     |  i,j     |    i,k    |					*
c     |          |           |					*
c     |__________|___________|					*
c     *
c     Furthermore, suppose that entry (i,j) resulted from an earlier	*
c     fill-in and has level s1, and (j,k) resulted from an earlier	*
c     fill-in and has level s2:						*
c     _______________________						*
c     |          |           |					*
c     |          |           |					*
c     | level 0  | level s2  |					*
c     |          |           |					*
c     |__________|___________|					*
c     |          |           |					*
c     |          |           |					*
c     | level s1 |           |					*
c     |          |           |					*
c     |__________|___________|					*
c     *
c     When using A(j,j) to annihilate A(i,j), fill-in will be incurred	*
c     in A(i,k).  How should its level be defined?  It would not be	*
c     operated on if A(i,j) or A(j,m) had not been filled in.  The 	*
c     version used here is to define its level as s1 + s2 + 1.  However,*
c     other reasonable choices would have been min(s1,s2) or max(s1,s2).*
c     Using the sum gives a more conservative strategy in terms of the	*
c     growth of the number of nonzeros as s increases.			*
c     *
c     levels(n+2:nzlu    ) stores the levels from previous rows,	*
c     that is, the s2's above.  levels(1:n) stores the fill-levels	*
c     of the current row (row i), which are the s1's above.		*
c     levels(n+1) is not used, so levels is conformant with MSR format.	*
c     *
c     Vectors used:							*
c     =============							*
c     *
c     lastcol(n):							*
c     The integer lastcol(k) is the row index of the last row		*
c     to have a nonzero in column k, including the current		*
c     row, and fill-in up to this point.  So for the matrix		*
c     *
c     |--------------------------|				*
c     | 11   12           15     |				*
c     | 21   22                26|				*
c     |      32  33   34         |				*
c     | 41       43   44         |				*
c     |      52       54  55   56|				*
c     |      62                66|				*
c     ---------------------------				*
c     *
c     after step 1, lastcol() = [1  0  0  0  1  0]		*
c     after step 2, lastcol() = [2  2  0  0  2  2]		*
c     after step 3, lastcol() = [2  3  3  3  2  3]		*
c     after step 4, lastcol() = [4  3  4  4  4  3]		*
c     after step 5, lastcol() = [4  5  4  5  5  5]		*
c     after step 6, lastcol() = [4  6  4  5  5  6]		*
c     *  
c     Note that on step 2, lastcol(5) = 2 because there is a	*
c     fillin position (2,5) in the matrix.  lastcol() is used	*
c     to determine if a nonzero occurs in column j because		*
c     it is a nonzero in the original matrix, or was a fill.		*
c     *  
c     rowll(n):								*
c     The integer vector rowll is used to keep a linked list of	        *
c     the nonzeros in the current row, allowing fill-in to be		*
c     introduced sensibly.  rowll is initialized with the		*
c     original nonzeros of the current row, and then sorted		*
c     using a shell sort.  A pointer called head         		*
c     (what ingenuity) is  initialized.  Note that at any		*
c     point rowll may contain garbage left over from previous		*
c     rows, which the linked list structure skips over.	        	*
c     For row 4 of the matrix above, first rowll is set to		*
c     rowll() = [3  1  2  5  -  -], where - indicates any integer.	*
c     Then the vector is sorted, which yields				*
c     rowll() = [1  2  3  5  -  -].  The vector is then expanded	*
c     to linked list form by setting head = 1  and         		*
c     rowll() = [2  3  5  -  7  -], where 7 indicates termination.	*
c     *  
c     ijlu(nzlu):							*
c     The returned nonzero structure for the LU factors.		*
c     This is built up row by row in MSR format, with both L		*
c     and U stored in the data structure.  Another vector, uptr(n),	*
c     is used to give pointers to the beginning of the upper		*
c     triangular part of the LU factors in ijlu.			*
c     *  
c     levels(n+2:nzlu):							*
c     This vector stores the fill level for each entry from		*
c     all the previous rows, used to compute if the current entry	*
c     will exceed the allowed levels of fill.  The value in		*
c     levels(m) is added to the level of fill for the element in	*
c     the current row that is being reduced, to figure if 		*
c     a column entry is to be accepted as fill, or rejected.		*
c     See the method explanation above.				        *
c     *  
c     levels(1:n):							*
c     This vector stores the fill level number for the current	        *
c     row's entries.  If they were created as fill elements		*
c     themselves, this number is added to the corresponding		*
c     entry in levels(n+2:nzlu) to see if a particular column		*
c     entry will							*
c     be created as new fill or not.  NOTE: in practice, the		*
c     value in levels(1:n) is one larger than the "fill" level of	*
c     the corresponding row entry, except for the diagonal		*
c     entry.  That is why the accept/reject test in the code		*
c     is "if (levels(j) + levels(m) .le. levfill + 1)".		        *
c     *  
c========================================================================
c     
c     on entry:
c==========
c     n       = The order of the matrix A.
c     ija     = Integer array. Matrix A stored in modified sparse row format.
c     levfill = Integer. Level of fill-in allowed.
c     nzmax   = Integer. The maximum number of nonzero entries in the
c     approximate factorization of a.  This is the amount of storage
c     allocated for ijlu.
c     
c     on return:
c===========
c     
c     nzlu   = The actual number of entries in the approximate factors, plus one.
c     ijlu   = Integer array of length nzlu containing pointers to 
c     delimit rows and specify column number for stored 
c     elements of the approximate factors of a.  the l 
c     and u factors are stored as one matrix.
c     uptr   = Integer array of length n containing the pointers to        
c     
c     ierr is an error flag:
c     ierr  = -i --> near zero pivot in step i
c     ierr  = 0  --> all's OK
c     ierr  = 1  --> not enough storage; check mneed.
c     ierr  = 2  --> illegal parameter
c     
c     mneed   = contains the actual number of elements in ldu, or the amount
c     of additional storage needed for ldu
c     
c     work arrays:
c=============
c     lastcol    = integer array of length n containing last update of the
c     corresponding column.
c     levels     = integer array of length n containing the level of
c     fill-in in current row in its first n entries, and
c     level of fill of previous rows of U in remaining part.
c     rowll      = integer array of length n containing pointers to implement a
c     linked list for the fill-in elements.
c     
c     
c     external functions:
c=============
c     ifix, float, min0, srtr
c     
c========================================================================

      integer n,colind(*),rwptr(*),ijlu(*),uptr(*),rowll(*), lastcol(*),
     1     levels(*), levfill,nzmax,nzlu
      integer  ierr,   mneed
      integer icolindj,ijlum,i,j,k,m,ibegin,iend,Ujbeg,Ujend
      integer head,prev,lm,actlev,lowct,k1,k2,levp1,lmk,nzi,rowct

c========================================================================
c     Beginning of Executable Statements
c========================================================================

c     --------------------------------------------------------------
c     Because the first row of the factor contains no strictly lower
c     triangular parts (parts of L), uptr(1) = ijlu(1) = n+2:
c     --------------------------------------------------------------
      ijlu(1)  =  n+2
      uptr(1)  =  n+2

c     --------------------------------------------------------
c     The storage for the nonzeros of LU must be at least n+1, 
c     for a diagonal matrix:
c     --------------------------------------------------------
      nzlu     =  n+1

c     --------------------------------------------------------------------
c     Number of allowed levels plus 1; used for the test of accept/reject.
c     See the notes about the methodology above.
c     --------------------------------------------------------------------
      levp1    =  levfill + 1
      
c     -------------------------------------------------------------
c     Initially, for all columns there were no nonzeros in the rows
c     above, because there are no rows above the first one.
c     -------------------------------------------------------------
      do i = 1,n
         lastcol(i) = 0
      end do
c     -------------------
c     Proceed row by row:
c     -------------------

      do 100 i = 1,n

c     ----------------------------------------------------------
c     Because the matrix diagonal entry is nonzero, the level of
c     fill for that diagonal entry is zero:
c     ----------------------------------------------------------
         levels(i) = 0

c     ----------------------------------------------------------
c     ibegin and iend are the beginning of rows i and i+1, resp.
c     ----------------------------------------------------------
         ibegin  =  rwptr(i)
         iend    =  rwptr(i+1)

c     -------------------------------------------------------------
c     Number of offdiagonal nonzeros in the original matrix's row i
c     -------------------------------------------------------------
         nzi   =  iend - ibegin
c     --------------------------------------------------------
c     If only the diagonal entry in row i is nonzero, skip the
c     fancy stuff; nothing need be done:
c     --------------------------------------------------------
         if (nzi .gt. 1) then

c     ----------------------------------------------------------
c     Decrement iend, so that it can be used as the ending index
c     in icolind of row i:
c     ----------------------------------------------------------
            iend          =  iend - 1

c     ---------------------------------------------------------
c     rowct keeps count of the number of nondiagonal entries in
c     the current row:
c     ---------------------------------------------------------
            rowct         =  0

c     ------------------------------------------------------------
c     For nonzeros in the current row from the original matrix A,
c     set lastcol to be the current row number, and the levels of
c     the entry to be 1.  Note that this is really the true level
c     of the element, plus 1.  At the same time, load up the work
c     array rowll with the column numbers for the original entries
c     from row i:
c     ------------------------------------------------------------
            do j = ibegin, iend
               icolindj           =  colind(j)
               lastcol(icolindj)  =  i
               if (icolindj .ne. i) then
                  levels(icolindj) = 1
                  rowct          =  rowct + 1
                  rowll(rowct)   =  icolindj
               end if
            end do
c     ---------------------------------------------------------
c     Sort the entries in rowll, so that the row has its column
c     entries in increasing order.
c     ---------------------------------------------------------
            call srtr(nzi-1,rowll)

c     ---------------------------------------------------------
c     Now set up rowll as a linked list containing the original
c     nonzero column numbers, as described in the methods section:
c     ---------------------------------------------------------
            head  =  rowll(1)
            k1    =  n+1
            do j = nzi-1, 1, -1
               k2        =  rowll(j)
               rowll(k2) =  k1
               k1        = k2
            end do
c     ------------------------------------------------------------
c     Increment count of nonzeros in the LU factors by the number
c     of nonzeros in the original matrix's row i.  Further
c     incrementing will be necessary if any fill-in actually occurs
c     ------------------------------------------------------------
            nzlu  =  nzlu + nzi - 1

c     ------------------------------------------------------------
c     The integer j will be used as a pointer to track through the
c     linked list rowll:
c     ------------------------------------------------------------
            j  =  head
c     
c     ------------------------------------------------------------
c     The integer lowct is used to keep count of the number of
c     nonzeros in the current row's strictly lower triangular part,
c     for setting uptr pointers to indicate where in ijlu the upperc
c     triangular part starts. 
c     ------------------------------------------------------------
            lowct =  0

c     ------------------------------------------------------------
c     Fill-in could only have resulted from rows preceding row i,
c     so we only need check those rows with index j < i.
c     Furthermore, if the current row has a zero in column j,
c     there is no need to check the preceding rows; there clearly
c     could not be any fill-in from those rows to this entry.
c     ------------------------------------------------------------
            do 80 while (j .lt. i)

c     ------------------------------------------------------------
c     Increment lower triangular part count, since in this case
c     (j<i) we got another entry in L:
c     ------------------------------------------------------------
               lowct = lowct  + 1

c     ---------------------------------------------------------
c     If the fill level is zero, there is no way to get fill in
c     occuring.  
c     ---------------------------------------------------------
               if (levfill .ne. 0) then

c     -----------------------------------------------------
c     Ujbeg is beginning index of strictly upper triangular
c     part of U's j-th row, and Ujend is the ending index
c     of it, in ijlu().
c     -----------------------------------------------------
                  Ujbeg = uptr(j)
                  Ujend = ijlu(j+1) - 1

c     -----------------------------------------------------
c     Need to set pointer to previous entry before working
c     segment of rowll, because if fill occurs that will be
c     a moving segment.
c     -----------------------------------------------------
                  prev  =  j

c     -----------------------------------------------------
c     lm is the next nonzero pointer in linked list rowll:
c     -----------------------------------------------------
                  lm    =  rowll(j)

c     -------------------------------------------------------
c     lmk is the fill level in this row, caused by
c     eliminating column entry j.  That is, level s1 from the
c     methodology explanation above.
c     -------------------------------------------------------
                  lmk   =  levels(j)

c     -------------------------------------------------------
c     Now proceed through the j-th row of U, because in the
c     elimination we add a multiple of it to row i to zero
c     out entry (i,j).  If a column entry in row j of U is
c     zero, there is no need to worry about fill, because it
c     cannot cause a fill in the corresponding entry of row i
c     -------------------------------------------------------
                  do 60  m = Ujbeg, Ujend

c     ----------------------------------------------------
c     ijlum is the column number of the current nonzero in
c     row j of U:
c     ----------------------------------------------------
                     ijlum =  ijlu(m)
c     
c     ---------------------------------------------------
c     actlev is the actual level (plus 1) of column entry
c     j in row i, from summing the level contributions
c     s1 and s2 as explained in the methods section.
c     Note that the next line could reasonably be
c     replaced by, e.g., actlev = max(lmk, levels(m)),
c     but this would cause greater fill-in:
c     ---------------------------------------------------
                     actlev = lmk + levels(m)

c     ---------------------------------------------------
c     If lastcol of the current column entry in U is not
c     equal to the current row number i, then the current
c     row has a zero in column j, and the earlier row j
c     in U has a nonzero, so possible fill can occur.
c     ---------------------------------------------------
                     if (lastcol(ijlum) .ne. i) then

c     --------------------------------------------------
c     If actlev < levfill + 1, then the new entry has an
c     acceptable fill level and needs to be added to the
c     data structure.
c     --------------------------------------------------
                        if (actlev .le. levp1) then

c     -------------------------------------------
c     Since the column entry ijlum in the current
c     row i is to be filled, we need to update
c     lastcol for that column number.  Also, the
c     level number of the current entry needs to be
c     set to actlev.  Note that when we finish 
c     processing this row, the n-vector levels(1:n)
c     will be copied over to the corresponding 
c     trailing part of levels, so that it can be
c     used in subsequent rows:
c     -------------------------------------------

                           lastcol(ijlum) = i
                           levels(ijlum) = actlev

c     -------------------------------------------
c     Now find location in the linked list rowll
c     where the fillin entry should be placed.
c     Chase through the linked list until the next
c     nonzero column is to the right of the fill
c     column number.
c     -------------------------------------------
                           do 50 while (lm .le. ijlum)
                              prev = lm
                              lm   = rowll(lm)
 50                        continue

c     -------------------------------------------
c     Insert new entry into the linked list for
c     row i, and increase the nonzero count for LU
c     -------------------------------------------
                           rowll(prev)  = ijlum
                           rowll(ijlum) = lm
                           prev       = ijlum
                           nzlu  =  nzlu  + 1
                        endif

c     -------------------------------------------------
c     Else clause is for when lastcol(ijlum) = i.  In
c     this case, the current column has a nonzero, but
c     it resulted from an earlier fill-in or from an
c     original matrix entry.  In this case, need to
c     update the level number for this column to be the
c     smaller of the two possible fill contributors,
c     the current fill number or the computed one from
c     updating this entry from a previous row.
c     -------------------------------------------------
                     else
                        levels(ijlum) = min0(levels(ijlum),actlev)
                     endif

c     -------------------------------------------------
c     Now go and pick up the next column entry from row
c     j of U:
c     -------------------------------------------------
 60               continue

c     -------------------------------------------
c     End if clause for levfill not equal to zero
c     -------------------------------------------
               endif

c     ------------------------------------------------------
c     Pick up next nonzero column index from the linked
c     list, and continue processing the i-th row's nonzeros.
c     This ends the first while loop (j < i).
c     ------------------------------------------------------
               j = rowll(j)

 80         continue

c     ---------------------------------------------------------
c     Check to see if we have exceeded the allowed memory
c     storage before storing the results of computing row i's
c     sparsity pattern into the ijlu and uptr data structures.
c     ---------------------------------------------------------
            if (nzlu .gt. nzmax) then
               mneed = ifix((float(n-i)/float(2*i))*3*nzlu)
               ierr  = 1
               return
            endif

c     ---------------------------------------------------------
c     Storage is adequate, so update ijlu data structure.
c     Row i ends at nzlu + 1:
c     ---------------------------------------------------------
            ijlu(i+1)   =  nzlu + 1

c     ---------------------------------------------------------
c     ... and the upper triangular part of LU begins at
c     lowct entries to right of where row i begins.
c     ---------------------------------------------------------
            uptr(i)     =  ijlu(i)  + lowct

c     -----------------------------------------------------
c     Now chase through linked list for row i, recording
c     information into ijlu.  At same time, put level data
c     into the levels array for use on later rows:
c     -----------------------------------------------------
            j  =  head
            k1 =  ijlu(i)
            do k  =  k1, nzlu
               ijlu(k)    =  j
               levels(k)  =  levels(j)
               j          =  rowll(j)
            end do

         else

c     ---------------------------------------------------------
c     This else clause ends the (nzi > 1) if.  If nzi = 1, then
c     the update of ijlu and uptr is trivial:
c     ---------------------------------------------------------
            ijlu(i+1)   =  nzlu + 1
            uptr(i)     =  ijlu(i)

         endif

c     ----------------------------------------------
c     And you thought we would never get through....
c     ----------------------------------------------

 100  continue
      
      ierr = 0
      return

c======================== End of symfac ==============================
      end
c     
      subroutine srtr(num,q)
      implicit none
c========================================================================
c========================================================================
c     
c     Implement shell sort, with hardwired increments.  The algorithm for
c     sorting entries in A(0:n-1) is as follows:
c----------------------------------------------------------------
c     inc = initialinc(n)
c     while inc >= 1
c     for i = inc to n-1
c     j = i
c     x = A(i)
c     while j >= inc and A(j-inc) > x
c     A(j) = A(j-inc)
c     j    = j-inc
c     end while 
c     A(j) = x
c     end for
c     inc = nextinc(inc,n)
c     end while
c----------------------------------------------------------------
c     
c     The increments here are 1, 4, 13, 40, 121, ..., (3**i - 1)/2, ...
c     In this case, nextinc(inc,n) = (inc-1)/3.  Usually shellsort
c     would have the largest increment the largest integer of the form
c     (3**i - 1)/2 that is less than n, but here it is fixed at 121
c     because most sparse matrices have 121 or fewer nonzero entries
c     per row.  If this routine is expanded for a complete sparse
c     factorization routine, or if a large number of levels of fill is
c     allowed, then possibly it should be replaced with more efficient
c     sorting.
c     
c     Any set of increments with 1 as the first one will result in a
c     true sorting algorithm.
c     
c========================================================================
c     
      integer num, q(num)
c     
      integer iinc(5), key, icn, ih, ii, i, j, jj
      data iinc/1,4,13,40,121/
c     
      if   (num .eq. 0) then
         icn   =  0
      elseif   (num .lt. 14) then
         icn   =  1
      elseif   (num .lt. 41) then
         icn   =  2
      elseif   (num .lt. 122) then
         icn   =  3
      elseif   (num .lt. 365) then
         icn   =  4
      else
         icn   =  5
      end if
      do 40 ii = 1, icn
         ih = iinc(icn + 1 - ii)
         do 30  j = ih+1, num
            i = j-ih
            key = q(j)
            do 10 jj = 1, j-ih, ih
               if (key .ge. q(i)) then
                  go to 20
               else
                  q(i + ih) = q(i)
                  i = i - ih
               end if
 10         continue
 20         continue
            q(i + ih) = key
 30      continue
 40   continue
c     
      return
      end
c     
c========================================================================
      subroutine  numfac(n, colind,rwptr,  jlu, uptr, a, lu,
     *     relax, ierr, colptrs, milu,nzlu,ijlu,luval)
c========================================================================
c     
c     Numerical factorization with given sparsity pattern identified in
c     the MSR data structure  jlu.
c     
c     Because of the MSR storage format, the method used is a deferred update
c     version.  At the k-th step, the updates from the previous rows m,
c     m < k and m corresponding to a nonzero in row k of L/U, are applied
c     to row k.  A relaxation parameter is used.  When it is 0.0, ILU
c     factorization is performed.  Otherwise the parameter is multiplied
c     times the discarded fillin and added to the diagonal entry.  So when 
c     the parameter is equal to 1.0, classical MILU factorization is performed.
c     
c     This routine implicitly assumes that the sparsity pattern specified
c     by the MSR data structure contains that of A.  Note the comments in
c     the code if you want to remove this restriction.
c     
c     Matrix A is input in CSR format, and the output MSR array contains
c     both L and U, with L assumed to be unit lower triangular.  The diagonal
c     entries of U are stored in inverted form to avoid divisions when
c     applying the preconditioner.
c     
c     Finally, if a small pivot element is found in the factorization, it
c     is replaced with a small quantity based on machine epsilon, the sign
c     of the pivot element, and the norm of the current row.
c     
c     
c     The process can be viewed best in the dense case, shown here
c     for k = 4 and n = 6.  Rows 1 through 3 have been computed and
c     are complete at the beginning of this stage (they need no
c     further updates).  The  corresponding dense factorization is
c     usually viewed as consisting of a triangular solve to find the
c     entries L(k,1:k-1), followed by a vector*matrix update to get
c     the entries U(k,k:n).  When a saxpy form is used for both the
c     triangular solve and the vector*matrix update, both phases have
c     the same form and can be computed with the same loop.
c     
c     Stage 0: Copy row 1 of A over to row 1 of LU.
c     
c     Stage 1: Compute LU(4,1), as LU(4,1)/U(1,1).
c     
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     | slv |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stage 2: Update LU(4,2:n) using multiplier in LU(4,1) and row 1 of LU:
c     
c     -------------------------------------
c     |     | U12 | U13 | U14 | U15 | U16 |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     | M   | upd | upd | upd | upd | upd |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stage 3: Compute LU(4,2), as LU(4,2)/U(2,2).
c     
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     | slv |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stage 4: Update LU(4,3:n) using multiplier in LU(4,2) and row 2 of LU:
c     
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     | U23 | U24 | U25 | U26 |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     | Ml  | upd | upd | upd | upd |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stage 5: Compute LU(4,3), as LU(4,3)/U(3,3).
c     
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     | slv |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stage 6: Update LU(4,4:n) using multiplier in LU(4,3) and row 3 of LU:
c     
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     | U34 | U35 | U36 |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     | M   | upd | upd | upd |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stage 7: Invert diagonal entry U(4,4), which is now complete:
c     
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     | Inv |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     |     |     |     |     |     |     |
c     |     |     |     |     |     |     |
c     -------------------------------------
c     
c     Stages 8 onwards: Go on to the next row, repeating above process.
c     
c     Note:  When applying updates from previous rows to the
c     current one, we only need apply those for which there is a
c     corresponding nonzero in the current row of LU; otherwise the
c     multiplier would be zero.
c     
c========================================================================
c     a, colind, rwptr    :  Matrix A in CSR format
c     lu, jlu, uptr: Data structure for L/U matrix in MSR format.
c     lu contains the nonzero values of L/U
c     lu(1:n) diagonal entries listed in order
c     lu(n+1) is unused,
c     lu(n+2:nzlu+1) are the off-diagonal entries in
c     column increasing order.
c     jlu contains pointers and column indices.
c     jlu(1:n) pointers to the beginning of each row
c     in arrays lu and jlu
c     jlu(n+1) is a pointer to one beyond the last
c     entry in lu and jlu
c     jlu(n+2:nnz+1) are the column indices
c     uptr(1:n) contains pointers the first entry of U in
c     each row.
c     n             : order of the matrices.
c     relax         : relaxation parameter for MILU methods.
c     mult          : temporary scalar for holding multiplier L(k,j).
c     ierr          : return error code.
c     ierr = 0  -> all's OK.
c     ierr < 0  -> row -ierr had a small pivot
c     k             : loop index for current row number
c     j             : loop index for current row number
c     indj          : index for traversing row k; index of j-th entry in
c     the lu and jlu arrays.
c     indja         : index for traversing row k of A.
c     inds          : index for traversing updating row of L/U when 
c     processing row k by applying sparse saxpy ops from
c     previous (updating) rows to row k (updated row).
c     jluj          : index j of current column number in row k
c     jlus          : index s of current column number in row j
c     ijaj          : column index j for entries of A.
c     
c     colptrs(1:n)     : used as an indicator array, to perform updates only
c     on columns that corr to allowed nonzeros.   If
c     column j is an allowed nonzero entry of the
c     current row, then colptrs(j) is the
c     index in arrays lu(:) and jlu(:) of
c     the corresponding entry (k,j) of L/U
c     milu             : logical indicating whether to relax diagonal entries
c     or not
c     rwnrm            : row norm of current row; used for determining small
c     diagonal element replacement.
c     nrw              : number of nonzeros in current row of A.
c     
c========================================================================

      implicit none

c     ------------------------
c     Declaration of arguments
c     ------------------------
      integer  n, colind(*), rwptr(*), jlu(*), uptr(*),ijlu(*)
      real*8   a(*), lu(*), relax, pilu,luval(*)
      integer  colptrs(n), ierr, nrw,nzlu

c     ---------------
c     Functions used
c     ---------------
c     real*8 dlamch
c     external dlamch

c     ---------------
c     Local variables
c     ---------------
      integer  k, indj, inds, indja
      integer  jluj, jlus, ijaj
      logical  milu
      real*8  SMALL
      real*8  rwnrm, mult

c========================================================================
c     Beginning of Executable Statements
c========================================================================

c     SMALL = sqrt(dlamch('E'))
      SMALL = 1.49d-8
      pilu=1.0
c     pilu=0.95
c     write(*,*) 'in numfac 1'
c     ------------------------------------------------------------
c     colptrs is used to hold the indices of entries in LU of 
c     row k.  It is initialized to zero here, and then reset after 
c     each row's work.
c     ------------------------------------------------------------
      do k =  1, n
         colptrs( k ) = 0
      end do
c     write(*,*) 'in numfac 2'
c     ------------------
c     Proceed row by row
c     ------------------
      do 200 k = 1, n

c     --------------------------------------------------------------
c     Set up colptrs with indices in lu of allowed nonzeros of row k
c     --------------------------------------------------------------
         do indj = jlu(k), jlu(k+1)-1
            colptrs(jlu(indj)) = indj
            lu(indj) = 0.0d0
         end do

c     ---------------------------------------------------------
c     Set the diagonal entry (not needed for CSR format of ILU)
c     ---------------------------------------------------------
         colptrs(k) =  k

c     ----------------------------------------------------------------
c     Copy row k of A over to LU.  Note that the "if" test in the loop
c     can be removed if it is known that the sparsity pattern of A
c     is contained in that of LU, which is the case for (M)ILU(s).
c     ----------------------------------------------------------------
         rwnrm = 0.0d0
         nrw = 0
         do indja = rwptr(k), rwptr(k+1)-1
            ijaj = colind(indja)
c     if (colptrs(ijaj) .ne. 0) then
            nrw = nrw + 1
            rwnrm =  rwnrm + abs(a(indja))
            lu(colptrs(ijaj)) = a(indja)
c     end if
         end do

c     -------------------------------------------------------------------
c     The first segment of the next loop on indj effectively solves
c     the transposed upper triangular system
c     U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'
c     via sparse saxpy operations, throwing away disallowed fill.
c     When the loop index indj reaches the k-th column (i.e., the
c     diagonal entry), then the innermost sparse saxpy operation 
c     effectively is applying the previous updates to the corresponding 
c     part of U via sparse vector*matrix, discarding disallowed fill-in
c     entries.  That operation is 
c     U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
c     -------------------------------------------------------------------

         do 100 indj = jlu(k), uptr(k)-1


c     -----------------------------------------------------------
c     jluj is the col number of current entry in row k of L,
c     and index of diagonal entry of U in same column.  For LU in
c     CSR format, that diag entry will require fancier indexing.
c     -----------------------------------------------------------
            jluj = jlu(indj)

c     -----------------------------------------------------------
c     Solve for next unknown in row k of L: L_kj = L_kj/U_jj ...
c     -----------------------------------------------------------
            lu(indj) = lu(indj)*lu(jluj)
            mult = lu(indj)

c     -------------------------------------------------------------
c     ... and use it as a multiplier to update the entries s in row
c     k of L, s = j+1, ... , k-1, and the entries s in row k of U,
c     s = k, ..., n.
c     -------------------------------------------------------------
            if (milu) then
               do inds = uptr(jluj), jlu(jluj+1)-1
                  jlus = jlu(inds)
                  if (colptrs(jlus) .ne. 0) then
                     lu(colptrs(jlus)) = lu(colptrs(jlus))
     &                    - mult*lu(inds)
                  else
                     lu(k) = lu(k) - relax*mult*lu(inds)
                  end if
               end do
            else
               do inds = uptr(jluj), jlu(jluj+1)-1
                  jlus = jlu(inds)
                  if (colptrs(jlus) .ne. 0) then
                     lu(colptrs(jlus)) = lu(colptrs(jlus))
     &                    - mult*lu(inds)*pilu
                  end if
               end do
            end if

 100     continue

c     ----------------------------------------------------------
c     Finished with row k of LU; reset colptrs indices to zero 
c     for next row, and invert diagonal entry of U for this row.
c     ----------------------------------------------------------
         do indj = jlu(k), jlu(k+1)-1
            colptrs(jlu(indj)) = 0
         end do
c     
         colptrs(k) =  0
c     
         if (abs(lu(k)) .le. SMALL*rwnrm/nrw) then
            lu(k) = sign(SMALL*rwnrm/nrw, lu(k))
         end if
         lu(k) = 1.0d0/lu(k)

 200  continue

      ierr  = 0
      
      do k=1,nzlu
         luval(k) = lu(k)
         ijlu(k) = jlu(k)
      end do

      return
c==================End of numfac =====================================
      end
c     
      subroutine qsplit(a,ind,n,ncut)
      real*8 a(n)
      integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c     
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c     
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
      real*8 tmp, abskey
      integer itmp, first, last
c-----
      first = 1
      last = n
      if (ncut .lt. first .or. ncut .gt. last) return
c     
c     outer loop -- while mid .ne. ncut do
c     
 1    mid = first
      abskey = abs(a(mid))
      do 2 j=first+1, last
         if (abs(a(j)) .gt. abskey) then
            mid = mid+1
c     interchange
            tmp = a(mid)
            itmp = ind(mid)
            a(mid) = a(j)
            ind(mid) = ind(j)
            a(j)  = tmp
            ind(j) = itmp
         endif
 2    continue
c     
c     interchange
c     
      tmp = a(mid)
      a(mid) = a(first)
      a(first)  = tmp
c     
      itmp = ind(mid)
      ind(mid) = ind(first)
      ind(first) = itmp
c     
c     test for while loop
c     
      if (mid .eq. ncut) return
      if (mid .gt. ncut) then
         last = mid-1
      else
         first = mid+1
      endif
      goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
      end
      
c======================================================================== 
      subroutine symbfactor(n,colind,rwptr,levfill,nzmax,nzlu,
     *           ijlu,uptr,ierr)
      implicit none

c======================================================================== 
c======================================================================== 
c                                                                      	*
c  Symbolic factorization of a matrix in compressed sparse row format, 	*
c    with resulting factors stored in a single MSR data structure.     	*
c                                                                      	*
c  This routine uses the CSR data structure of A in two integer vectors	*
c    colind, rwptr to set up the data structure for the ILU(levfill) 	*
c    factorization of A in the integer vectors ijlu and uptr.  Both L	*
c    and U are stored in the same structure, and uptr(i) is the pointer	*
c    to the beginning of the i-th row of U in ijlu.			*
c                                                                      	*
c======================================================================== 
c                                                                      	*
c    Method Used                                                       	*
c    ===========                                                      	*
c                                                                      	*
c  The implementation assumes that the diagonal entries are		*
c  nonzero, and remain nonzero throughout the elimination		*
c  process.  The algorithm proceeds row by row.  When computing		*
c  the sparsity pattern of the i-th row, the effect of row		*
c  operations from previous rows is considered.  Only those		*
c  preceding rows j for which (i,j) is nonzero need be considered,	*
c  since otherwise we would not have formed a linear combination	*
c  of rows i and j.							*
c                                                                      	*
c  The method used has some variations possible.  The definition	*
c  of ILU(s) is not well specified enough to get a factorization	*
c  that is uniquely defined, even in the sparsity pattern that		*
c  results.  For s = 0 or 1, there is not much variation, but for	*
c  higher levels of fill the problem is as follows:  Suppose		*
c  during the decomposition while computing the nonzero pattern		*
c  for row i the following principal submatrix is obtained:		*
c       _______________________						*
c       |          |           |					*
c       |          |           |					*
c       |  j,j     |    j,k    |					*
c       |          |           |					*
c       |__________|___________|					*
c       |          |           |					*
c       |          |           |					*
c       |  i,j     |    i,k    |					*
c       |          |           |					*
c       |__________|___________|					*
c  									*
c  Furthermore, suppose that entry (i,j) resulted from an earlier	*
c  fill-in and has level s1, and (j,k) resulted from an earlier		*
c  fill-in and has level s2:						*
c       _______________________						*
c       |          |           |					*
c       |          |           |					*
c       | level 0  | level s2  |					*
c       |          |           |					*
c       |__________|___________|					*
c       |          |           |					*
c       |          |           |					*
c       | level s1 |           |					*
c       |          |           |					*
c       |__________|___________|					*
c  									*
c  When using A(j,j) to annihilate A(i,j), fill-in will be incurred	*
c  in A(i,k).  How should its level be defined?  It would not be	*
c  operated on if A(i,j) or A(j,m) had not been filled in.  The 	*
c  version used here is to define its level as s1 + s2 + 1.  However,	*
c  other reasonable choices would have been min(s1,s2) or max(s1,s2).	*
c  Using the sum gives a more conservative strategy in terms of the	*
c  growth of the number of nonzeros as s increases.			*
c  									*
c  levels(n+2:nzlu    ) stores the levels from previous rows,		*
c  that is, the s2's above.  levels(1:n) stores the fill-levels		*
c  of the current row (row i), which are the s1's above.		*
c  levels(n+1) is not used, so levels is conformant with MSR format.	*
c  									*
c  Vectors used:							*
c  =============							*
c  									*
c  lastcol(n):								*
c  	The integer lastcol(k) is the row index of the last row		*
c  	to have a nonzero in column k, including the current		*
c  	row, and fill-in up to this point.  So for the matrix		*
c  									*
c             |--------------------------|				*
c             | 11   12           15     |				*
c             | 21   22                26|				*
c             |      32  33   34         |				*
c             | 41       43   44         |				*
c             |      52       54  55   56|				*
c             |      62                66|				*
c             ---------------------------				*
c  									*
c             after step 1, lastcol() = [1  0  0  0  1  0]		*
c             after step 2, lastcol() = [2  2  0  0  2  2]		*
c             after step 3, lastcol() = [2  3  3  3  2  3]		*
c             after step 4, lastcol() = [4  3  4  4  4  3]		*
c             after step 5, lastcol() = [4  5  4  5  5  5]		*
c             after step 6, lastcol() = [4  6  4  5  5  6]		*
c									*  
c          Note that on step 2, lastcol(5) = 2 because there is a	*
c          fillin position (2,5) in the matrix.  lastcol() is used	*
c   	to determine if a nonzero occurs in column j because		*
c   	it is a nonzero in the original matrix, or was a fill.		*
c									*  
c  rowll(n):								*
c  	The integer vector rowll is used to keep a linked list of	*
c  	the nonzeros in the current row, allowing fill-in to be		*
c   	introduced sensibly.  rowll is initialized with the		*
c  	original nonzeros of the current row, and then sorted		*
c  	using a shell sort.  A pointer called head         		*
c  	(what ingenuity) is  initialized.  Note that at any		*
c  	point rowll may contain garbage left over from previous		*
c  	rows, which the linked list structure skips over.		*
c  	For row 4 of the matrix above, first rowll is set to		*
c   	rowll() = [3  1  2  5  -  -], where - indicates any integer.	*
c   	Then the vector is sorted, which yields				*
c   	rowll() = [1  2  3  5  -  -].  The vector is then expanded	*
c  	to linked list form by setting head = 1  and         		*
c   	rowll() = [2  3  5  -  7  -], where 7 indicates termination.	*
c									*  
c  ijlu(nzlu):								*
c  	The returned nonzero structure for the LU factors.		*
c  	This is built up row by row in MSR format, with both L		*
c  	and U stored in the data structure.  Another vector, uptr(n),	*
c  	is used to give pointers to the beginning of the upper		*
c  	triangular part of the LU factors in ijlu.			*
c									*  
c  levels(n+2:nzlu):							*
c  	This vector stores the fill level for each entry from		*
c  	all the previous rows, used to compute if the current entry	*
c  	will exceed the allowed levels of fill.  The value in		*
c  	levels(m) is added to the level of fill for the element in	*
c   	the current row that is being reduced, to figure if 		*
c  	a column entry is to be accepted as fill, or rejected.		*
c  	See the method explanation above.				*
c									*  
c  levels(1:n):								*
c  	This vector stores the fill level number for the current	*
c  	row's entries.  If they were created as fill elements		*
c  	themselves, this number is added to the corresponding		*
c  	entry in levels(n+2:nzlu) to see if a particular column		*
c       entry will							*
c  	be created as new fill or not.  NOTE: in practice, the		*
c  	value in levels(1:n) is one larger than the "fill" level of	*
c  	the corresponding row entry, except for the diagonal		*
c  	entry.  That is why the accept/reject test in the code		*
c  	is "if (levels(j) + levels(m) .le. levfill + 1)".		*
c									*  
c======================================================================== 
c
c on entry:
c========== 
c  n       = The order of the matrix A.
c  ija     = Integer array. Matrix A stored in modified sparse row format.
c  levfill = Integer. Level of fill-in allowed.
c  nzmax   = Integer. The maximum number of nonzero entries in the
c           approximate factorization of a.  This is the amount of storage
c           allocated for ijlu.
c
c on return:
c=========== 
c
c nzlu   = The actual number of entries in the approximate factors, plus one.
c ijlu   = Integer array of length nzlu containing pointers to 
c           delimit rows and specify column number for stored 
c           elements of the approximate factors of a.  the l 
c           and u factors are stored as one matrix.
c uptr   = Integer array of length n containing the pointers to        
c
c ierr is an error flag:
c        ierr  = -i --> near zero pivot in step i
c        ierr  = 0  --> all's OK
c        ierr  = 1  --> not enough storage; check mneed.
c        ierr  = 2  --> illegal parameter
c
c mneed   = contains the actual number of elements in ldu, or the amount
c		of additional storage needed for ldu
c
c work arrays:
c=============
c lastcol    = integer array of length n containing last update of the
c              corresponding column.
c levels     = integer array of length n containing the level of
c              fill-in in current row in its first n entries, and
c              level of fill of previous rows of U in remaining part.
c rowll      = integer array of length n containing pointers to implement a
c              linked list for the fill-in elements.
c
c
c external functions:
c====================
c ifix, float, min0, srtr
c
c======================================================================== 


      integer n,colind(*),rwptr(*),ijlu(*),uptr(*),rowll(n), lastcol(n),
     1        levfill,nzmax,nzlu, levels(nzmax)
      integer  ierr,   mneed
      integer icolindj,ijlum,i,j,k,m,ibegin,iend,Ujbeg,Ujend
      integer head,prev,lm,actlev,lowct,k1,k2,levp1,lmk,nzi,rowct
       logical cinindex

c======================================================================== 
c       Beginning of Executable Statements
c======================================================================== 

       cinindex = ( rwptr(1) .eq. 0 )

      if (cinindex) then 
         do i=1,n+1
            rwptr(i)=rwptr(i)+1
         end do
         do i = 1,rwptr(n+1)-1
            colind(i)=colind(i)+1
         end do
      end if

c     --------------------------------------------------------------
c     Because the first row of the factor contains no strictly lower
c     triangular parts (parts of L), uptr(1) = ijlu(1) = n+2:
c     --------------------------------------------------------------
      ijlu(1)  =  n+2
      uptr(1)  =  n+2

c     --------------------------------------------------------
c     The storage for the nonzeros of LU must be at least n+1, 
c     for a diagonal matrix:
c     --------------------------------------------------------
      nzlu     =  n+1

c     --------------------------------------------------------------------
c     Number of allowed levels plus 1; used for the test of accept/reject.
c     See the notes about the methodology above.
c     --------------------------------------------------------------------
      levp1    =  levfill + 1

c     -------------------------------------------------------------
c     Initially, for all columns there were no nonzeros in the rows
c     above, because there are no rows above the first one.
c     -------------------------------------------------------------
      do i = 1,n
      	lastcol(i) = 0
      end do

c     -------------------
c     Proceed row by row:
c     -------------------

      do 100 i = 1,n

c       ----------------------------------------------------------
c       Because the matrix diagonal entry is nonzero, the level of
c       fill for that diagonal entry is zero:
c       ----------------------------------------------------------
      	levels(i) = 0

c       ----------------------------------------------------------
c       ibegin and iend are the beginning of rows i and i+1, resp.
c       ----------------------------------------------------------
      	ibegin    =  rwptr(i)
      	iend    =  rwptr(i+1)

c       -------------------------------------------------------------
c       Number of offdiagonal nonzeros in the original matrix's row i
c       -------------------------------------------------------------
      	nzi   =  iend - ibegin

c       --------------------------------------------------------
c       If only the diagonal entry in row i is nonzero, skip the
c       fancy stuff; nothing need be done:
c       --------------------------------------------------------
      	if (nzi .gt. 1) then

c           ----------------------------------------------------------
c           Decrement iend, so that it can be used as the ending index
c           in icolind of row i:
c           ----------------------------------------------------------
            iend          =  iend - 1

c           ---------------------------------------------------------
c           rowct keeps count of the number of nondiagonal entries in
c           the current row:
c           ---------------------------------------------------------
            rowct          =  0

c           ------------------------------------------------------------
c           For nonzeros in the current row from the original matrix A,
c           set lastcol to be the current row number, and the levels of
c           the entry to be 1.  Note that this is really the true level
c           of the element, plus 1.  At the same time, load up the work
c           array rowll with the column numbers for the original entries
c           from row i:
c           ------------------------------------------------------------
            do j = ibegin, iend
               icolindj           =  colind(j)
               lastcol(icolindj)  =  i
               if (icolindj .ne. i) then
                  levels(icolindj)   =  1
                  rowct          =  rowct + 1
                  rowll(rowct)   =  icolindj
               end if
            end do

c           ---------------------------------------------------------
c           Sort the entries in rowll, so that the row has its column
c           entries in increasing order.
c           ---------------------------------------------------------
            call srtr(nzi-1,rowll)

c           ---------------------------------------------------------
c           Now set up rowll as a linked list containing the original
c           nonzero column numbers, as described in the methods section:
c           ---------------------------------------------------------
            head  =  rowll(1)
            k1    =  n+1
            do j = nzi-1, 1, -1
               k2        =  rowll(j)
               rowll(k2) =  k1
               k1        = k2
            end do

c           ------------------------------------------------------------
c           Increment count of nonzeros in the LU factors by the number
c           of nonzeros in the original matrix's row i.  Further
c           incrementing will be necessary if any fill-in actually occurs
c           ------------------------------------------------------------
            nzlu  =  nzlu + nzi - 1

c           ------------------------------------------------------------
c           The integer j will be used as a pointer to track through the
c           linked list rowll:
c           ------------------------------------------------------------
            j  =  head
c
c           ------------------------------------------------------------
c           The integer lowct is used to keep count of the number of
c           nonzeros in the current row's strictly lower triangular part,
c           for setting uptr pointers to indicate where in ijlu the upperc
c           triangular part starts. 
c           ------------------------------------------------------------
            lowct =  0

c           ------------------------------------------------------------
c           Fill-in could only have resulted from rows preceding row i,
c           so we only need check those rows with index j < i.
c           Furthermore, if the current row has a zero in column j,
c           there is no need to check the preceding rows; there clearly
c           could not be any fill-in from those rows to this entry.
c           ------------------------------------------------------------
            do 80 while (j .lt. i)

c              ------------------------------------------------------------
c              Increment lower triangular part count, since in this case
c              (j<i) we got another entry in L:
c              ------------------------------------------------------------
               lowct = lowct  + 1

c              ---------------------------------------------------------
c              If the fill level is zero, there is no way to get fill in
c              occuring.  
c              ---------------------------------------------------------
               if (levfill .ne. 0) then

c                 -----------------------------------------------------
c                 Ujbeg is beginning index of strictly upper triangular
c                 part of U's j-th row, and Ujend is the ending index
c                 of it, in ijlu().
c                 -----------------------------------------------------
                  Ujbeg = uptr(j)
                  Ujend = ijlu(j+1) - 1

c                 -----------------------------------------------------
c                 Need to set pointer to previous entry before working
c                 segment of rowll, because if fill occurs that will be
c                 a moving segment.
c                 -----------------------------------------------------
                  prev  =  j

c                 -----------------------------------------------------
c                 lm is the next nonzero pointer in linked list rowll:
c                 -----------------------------------------------------
                  lm    =  rowll(j)

c                 -------------------------------------------------------
c                 lmk is the fill level in this row, caused by
c                 eliminating column entry j.  That is, level s1 from the
c                 methodology explanation above.
c                 -------------------------------------------------------
                  lmk   =  levels(j)

c                 -------------------------------------------------------
c                 Now proceed through the j-th row of U, because in the
c                 elimination we add a multiple of it to row i to zero
c                 out entry (i,j).  If a column entry in row j of U is
c                 zero, there is no need to worry about fill, because it
c                 cannot cause a fill in the corresponding entry of row i
c                 -------------------------------------------------------
                  do 60  m = Ujbeg, Ujend

c                    ----------------------------------------------------
c                    ijlum is the column number of the current nonzero in
c                    row j of U:
c                    ----------------------------------------------------
                     ijlum =  ijlu(m)
c
c                    ---------------------------------------------------
c                    actlev is the actual level (plus 1) of column entry
c                    j in row i, from summing the level contributions
c                    s1 and s2 as explained in the methods section.
c                    Note that the next line could reasonably be
c                    replaced by, e.g., actlev = max(lmk, levels(m)),
c                    but this would cause greater fill-in:
c                    ---------------------------------------------------
                     actlev = lmk + levels(m)

c                    ---------------------------------------------------
c                    If lastcol of the current column entry in U is not
c                    equal to the current row number i, then the current
c                    row has a zero in column j, and the earlier row j
c                    in U has a nonzero, so possible fill can occur.
c                    ---------------------------------------------------
                     if (lastcol(ijlum) .ne. i) then

c                    --------------------------------------------------
c                    If actlev < levfill + 1, then the new entry has an
c                    acceptable fill level and needs to be added to the
c                    data structure.
c                    --------------------------------------------------
                        if (actlev .le. levp1) then

c                          -------------------------------------------
c                          Since the column entry ijlum in the current
c                          row i is to be filled, we need to update
c                          lastcol for that column number.  Also, the
c                          level number of the current entry needs to be
c                          set to actlev.  Note that when we finish 
c                          processing this row, the n-vector levels(1:n)
c                          will be copied over to the corresponding 
c                          trailing part of levels, so that it can be
c                          used in subsequent rows:
c                          -------------------------------------------

                           lastcol(ijlum) = i
                           levels(ijlum) = actlev

c                          -------------------------------------------
c                          Now find location in the linked list rowll
c                          where the fillin entry should be placed.
c                          Chase through the linked list until the next
c                          nonzero column is to the right of the fill
c                          column number.
c                          -------------------------------------------
                           do 50 while (lm .le. ijlum)
                              prev = lm
                              lm   = rowll(lm)
 50                        continue

c                          -------------------------------------------
c                          Insert new entry into the linked list for
c                          row i, and increase the nonzero count for LU
c                          -------------------------------------------
                           rowll(prev)  = ijlum
                           rowll(ijlum) = lm
                           prev       = ijlum
                           nzlu  =  nzlu  + 1
                        endif

c                    -------------------------------------------------
c                    Else clause is for when lastcol(ijlum) = i.  In
c                    this case, the current column has a nonzero, but
c                    it resulted from an earlier fill-in or from an
c                    original matrix entry.  In this case, need to
c                    update the level number for this column to be the
c                    smaller of the two possible fill contributors,
c                    the current fill number or the computed one from
c                    updating this entry from a previous row.
c                    -------------------------------------------------
                     else
                        levels(ijlum) = min0(levels(ijlum),actlev)
                     endif

c                  -------------------------------------------------
c                  Now go and pick up the next column entry from row
c                  j of U:
c                  -------------------------------------------------
 60                continue

c              -------------------------------------------
c              End if clause for levfill not equal to zero
c              -------------------------------------------
               endif

c              ------------------------------------------------------
c              Pick up next nonzero column index from the linked
c              list, and continue processing the i-th row's nonzeros.
c              This ends the first while loop (j < i).
c              ------------------------------------------------------
               j = rowll(j)

 80         continue

c           ---------------------------------------------------------
c           Check to see if we have exceeded the allowed memory
c           storage before storing the results of computing row i's
c           sparsity pattern into the ijlu and uptr data structures.
c           ---------------------------------------------------------
            if (nzlu .gt. nzmax) then
               mneed = ifix((float(n-i)/float(2*i))*3*nzlu)
               ierr  = 1
               return
            endif

c           ---------------------------------------------------------
c           Storage is adequate, so update ijlu data structure.
c           Row i ends at nzlu + 1:
c           ---------------------------------------------------------
            ijlu(i+1)   =  nzlu + 1

c           ---------------------------------------------------------
c           ... and the upper triangular part of LU begins at
c           lowct entries to right of where row i begins.
c           ---------------------------------------------------------
            uptr(i)     =  ijlu(i)  + lowct

c           -----------------------------------------------------
c           Now chase through linked list for row i, recording
c           information into ijlu.  At same time, put level data
c           into the levels array for use on later rows:
c           -----------------------------------------------------
            j  =  head
            k1 =  ijlu(i)
            do k  =  k1, nzlu
               ijlu(k)    =  j
               levels(k)  =  levels(j)
               j          =  rowll(j)
            end do

         else

c           ---------------------------------------------------------
c           This else clause ends the (nzi > 1) if.  If nzi = 1, then
c           the update of ijlu and uptr is trivial:
c           ---------------------------------------------------------
            ijlu(i+1)   =  nzlu + 1
            uptr(i)     =  ijlu(i)

         endif

c        ----------------------------------------------
c        And you thought we would never get through....
c        ----------------------------------------------

 100  continue
c      write(*,*) 'nzlu ==',nzlu
      if (cinindex) then
         do i=1,nzlu
            ijlu(i)=ijlu(i)-1
         end do
         do i=1,n
            uptr(i)=uptr(i)-1
         end do
         do i = 1,rwptr(n+1)-1
            colind(i)=colind(i)-1
         end do
        do i=1,n+1
            rwptr(i)=rwptr(i)-1
         end do
      end if 

      ierr = 0
      return

c======================== End of symbfac ==============================

      end
      
