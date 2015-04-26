c     
c     ilu.f
c     
c     WARNING: The subroutines contained in this file are modified from 
c     SPARSEKIT 2.0 and they are NOT part of the FASP package.
c
c     SPARSKIT is a package of FORTRAN subroutines for working with
c     sparse matrices. It includes general sparse matrix manipulation
c     routines as well as a few iterative solvers.
c
c     SPARSKIT is free software; you can redistribute it and/or modify 
c     it under the terms of the  GNU Lesser General Public License as 
c     published by the Free Software Foundation [version 2.1 of the 
c     License, or any later version.]
c
c-----------------------------------------------------------------------
c     Created  by Shiquan Zhang  on 12/27/2009.
c     Modified by Chensong Zhang on 05/25/2010.
c     Modified by Chensong Zhang on 04/21/2012.
c     Modified by Chunsheng Feng on 04/26/2015.
c-----------------------------------------------------------------------

!> \file ilu.f
!> \brief ILU routines for preconditioning adapted from SPARSEKIT
!>
!> \note Incomplete Factorization Methods: ILUk, ILUt, ILUtp

c-----------------------------------------------------------------------
      subroutine iluk(n,a,ja,ia,lfil,alu,jlu,iwk,ierr,nzlu)
c-----------------------------------------------------------------------
      implicit none 
      integer n,nzlu,lfil,iwk,ierr
      real*8 a(1),alu(1)
c     !,w(n)
      integer ja(1),ia(1),jlu(1)
c     !,ju(n),levs(iwk),jw(3*n)
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
c     jlu     = integer array of length n containing the pointers to
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
      real*8, allocatable:: w(:) 
      integer,allocatable:: ju(:),jw(:),levs(:)

      if (lfil .lt. 0) goto 998
      
      allocate(w(n))
      allocate(ju(n))
      allocate(jw(3*n))
      allocate(levs(iwk))
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

      deallocate(w)
      deallocate(ju)
      deallocate(jw)
      deallocate(levs)

      ierr = 0
      return
c     
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      deallocate(w)
      deallocate(ju)
      deallocate(jw)
      deallocate(levs)
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      deallocate(w)
      deallocate(ju)
      deallocate(jw)
      deallocate(levs)
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      deallocate(w)
      deallocate(ju)
      deallocate(jw)
      deallocate(levs)
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      deallocate(w)
      deallocate(ju)
      deallocate(jw)
      deallocate(levs)
      return
c     
c     zero row encountered in A or U. 
c     
 999  ierr = -5
      deallocate(w)
      deallocate(ju)
      deallocate(jw)
      deallocate(levs)
      return
c----------------end-of-iluk--------------------------------------------
c-----------------------------------------------------------------------
      end
c     
c-----------------------------------------------------------------------
      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,iwk,ierr,nz)
c-----------------------------------------------------------------------
!$      use omp_lib
      implicit none 
      integer n 
      real*8 a(1),alu(1),droptol
c     !,w(n+1)
      integer ja(1),ia(1),jlu(1),lfil,iwk,ierr,nz
c      !,ju(n),jw(2*n)
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
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,NE,len
      real*8 t, abs, s, fact, tmp 
c      !tnorm(n)
      logical cinindex
      real*8, allocatable:: w(:),tnorm(:) 
      integer,allocatable:: ju(:),jw(:)

      if (lfil .lt. 0) goto 998

      allocate(ju(n))
      allocate(jw(2*n))
      allocate(w(n+1))
      allocate(tnorm(n))

      cinindex = ia(1) .eq. 0

      if (cinindex) then
         NE = n+1       !modify by chunsheng 2012,Sep,1
!$OMP PARALLEL DO PRIVATE(I)              
         do i=1,NE
            ia(i)=ia(i)+1
         end do
!$OMP END PARALLEL DO
         NE = ia(n+1)-1
!$OMP PARALLEL DO PRIVATE(I)              
         do i = 1,NE
            ja(i)=ja(i)+1
         end do
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(J)              
      do 1 j=1,n
         jw(n+j) = 0
 1    continue
!$OMP END PARALLEL DO
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(II,j1,j2,k,tmp)
      do ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tmp = 0.0d0
         do k=j1,j2
            tmp = tmp+abs(a(k))
         enddo
!         if (tnorm(ii) .eq. 0.0) goto 999
         tmp = tmp/real(j2-j1+1)
         tnorm(ii) = tmp*droptol;
      enddo
!$OMP END PARALLEL DO

      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
!         tnorm = 0.0d0
!         do 501 k=j1,j2
!            tnorm = tnorm+abs(a(k))
! 501     continue
!         if (tnorm(ii) .eq. 0.0) goto 999
!         if (tnorm .eq. 0.0) goto 999
!         tnorm = tnorm/real(j2-j1+1)
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
!            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
            if (abs(w(ii+k)) .gt. tnorm(ii)) then 
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
c         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
         if (w(ii) .eq. 0.0) w(ii) = tnorm(ii)
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
!$OMP PARALLEL DO PRIVATE(I)              
         do i=1,nz
            jlu(i)=jlu(i)-1
         end do
!$OMP END PARALLEL DO               
      end if

      deallocate(ju)
      deallocate(jw)
      deallocate(w)
      deallocate(tnorm)

      ierr = 0
      return
c     
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      deallocate(ju)
      deallocate(jw)
      deallocate(w)
      deallocate(tnorm)
      write(*,*) 'input matrix may be wrong'
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      deallocate(ju)
      deallocate(jw)
      deallocate(w)
      deallocate(tnorm)
      write(*,*),'insufficient storage in L'
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      deallocate(ju)
      deallocate(jw)
      deallocate(w)
      deallocate(tnorm)
      write(*,*),'insufficient storage in U'
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      deallocate(ju)
      deallocate(jw)
      deallocate(w)
      deallocate(tnorm)
      write(*,*),'illegal lfil entered'
      return
c     
c     zero row encountered
c     
 999  ierr = -5
      deallocate(ju)
      deallocate(jw)
      deallocate(w)
      deallocate(tnorm)
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
      integer n,ja(1),ia(1),lfil,jlu(1),iwk,ierr
c     !ju(n),jw(2*n),iperm(2*n) 
      real*8 a(1), alu(1), droptol
c      !w(n+1)
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
      real*8, allocatable:: w(:) 
      integer,allocatable:: ju(:),jw(:),iperm(:)
      
      if (lfil .lt. 0) goto 998
      
      allocate(ju(n))
      allocate(jw(2*n))
      allocate(iperm(2*n))
      allocate(w(n+1))

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
      deallocate(ju)
      deallocate(jw)
      deallocate(iperm)
      deallocate(w)
      return
c     
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      deallocate(ju)
      deallocate(jw)
      deallocate(iperm)
      deallocate(w)
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      deallocate(ju)
      deallocate(jw)
      deallocate(iperm)
      deallocate(w)
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      deallocate(ju)
      deallocate(jw)
      deallocate(iperm)
      deallocate(w)
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      deallocate(ju)
      deallocate(jw)
      deallocate(iperm)
      deallocate(w)
      return
c     
c     zero row encountered
c     
 999  ierr = -5
      deallocate(ju)
      deallocate(jw)
      deallocate(iperm)
      deallocate(w)
      return
c----------------end-of-ilutp-------------------------------------------
      end
c     
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
c----------------end-of-srtr-------------------------------------------
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
c uptr   = Integer array of length n containing the pointers to upper trig
c           matrix       
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
      
