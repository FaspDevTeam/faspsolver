! These routines are part of the matching MG method
! Copyright (C) Ludmil Zikatanov 2006/2007
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! DOES NOT WORK FOR NONSYM PATTERN

!=====================================================================
      subroutine cut0(n,ia,ja,a,iaw,jaw,
     >     jblk,iblk,nblk,lwork1,lwork2,lwork3,msize)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension jblk(n),ia(*),ja(*),iaw(*),jaw(*),iblk(n),
     >     lwork1(n),lwork2(n),lwork3(n)
      dimension a(*)
C---------------------------------------------------------------------
C...  INPUT is ia and ja; output is "jblk()" which is the permutation"
C... iblk is the same here as iblk; iord is just jblk.
C... msize is input and is the maximal block size
C---------------------------------------------------------------------
      tol = 1.d-12
      nnz = ia(n+1)-1
C
      call icopyv(ia,iaw,n+1)
      call icopyv(ja,jaw,nnz)
C
      do 30 i = 2, n
         iaa = ia(i)
         iab = ia(i+1)-1
         do 20 ii = iaa,iab
            j = ja(ii)
            if (i .le. j) go to 20
            iac = ia(j)
            iad = ia(j+1)-1
            do 10 iik = iac, iad
               iii = iik
               if (ja(iii) .eq. i) go to 15
 10         continue
c            write(*,*) ' Error j has no i, (i,j) = (',i,j,
c     >           '): STOPPING because of Nonsymmetric pattern'
c            stop 11
            iii=0
            aa=a(ii)
            bb=0d0
            go to 100
 15         continue
            aa = a(ii)
            bb = a(iii)
            if(i .ne. j .and. aa .gt. tol) aa = 0.0d00 ! do not consider positive...
            if(i .ne. j .and. bb .gt. tol) bb = 0.0d00 ! do not consider positive...
            if(dabs(aa) .lt. tol .and. dabs(bb) .lt. tol) then
               jaw(ii) = 0
               jaw(iii) = 0
               go to 20
            end if
 100        continue
            call chsize(aa,bb,tol,imin)
            if (imin .eq. 0) go to 20
            if (imin .eq. 1) then
               jaw(ii) = 0
            else
               if(iii .ne. 0) jaw(iii) = 0
            end if
 20      continue
 30   continue
C
      call shift(iaw,jaw,n)
C
C      call lpri(iaw,jaw,n)
      call dfs(n,iaw,jaw,nblk,iblk,jblk,
     >         lwork1,lwork2,lwork3)
C
      mxyz = 0
      do k = 1 , nblk
         n0 = iblk(k)
         n1 = iblk(k+1)
         mxyz = max0(mxyz,n1-n0)
c         if(n1-n0 .eq. 121) then
c            write(*,*) k
c            read(*,*)
c c        end if
      end do
      write(*,*) ' Ordering:::', nblk, ' blocks.', 
     >     ' Max size = ', mxyz
ccc      return
!     max block size
C
C... put the blocks now in iaw, jaw, 
C... if there are any that are big. 
C... this below works for nonoverlapping blocks. 
      do k = 1 , nblk
         n0 = iblk(k)
         n1 = iblk(k+1) - 1
         do jk = n0,n1
            inode = jblk(k)
            lwork1(inode) = k
         end do
      end do
      iaw(1) = 1
      iptr = 1
      nblk1 = 0
      do k = 1,nblk
         n0 = iblk(k)
         n1 = iblk(k+1) - 1
         if(n1 - n0 + 1 .gt. msize) then
            kloc = 0
            do k01 = n0 , n1
               i = jblk(k01)
               numi = lwork1(i)
               iaa = ia(i)
               iab = ia(i+1)-1
               do jk = iaa,iab
                  j = ja(jk)
                  numj = lwork1(j)
                  if(numi .eq. numj) then
                     jaw(iptr) = j
                     iptr = iptr + 1
                  end if
               end do
               nblk1 = nblk1 + 1
               iaw(nblk1+1) = iptr
            end do 
         else
            do k01 = n0 , n1
               jaw(iptr) =  jblk(k01)
               iptr = iptr + 1
            end do
            nblk1 = nblk1 + 1
            iaw(nblk1+1) = iptr
         end if
      end do
      nblk = nblk1
      do k = 1 , nblk+1
         iblk(k) = iaw(k)
      end do
      do k = 1 , iblk(nblk+1)-1
         jblk(k) = jaw(k)
      end do
CCCCCCCCCCC
cccccc
C
      mxyz = 0
      do k = 1 , nblk
         n0 = iblk(k)
         n1 = iblk(k+1)
         mxyz = max0(mxyz,n1-n0)
c         if(n1-n0 .eq. 121) then
c            write(*,*) k
c            read(*,*)
c c        end if
      end do
      write(*,*) ' Second Pass:::', nblk, ' blocks.', 
     >     ' Max size = ', mxyz
cccccccccccccccccccccccccccccccccc
      return
      end
C======================================================================
      subroutine chsize(a,b,tol,imin)
C======================================================================
      implicit real*8(a-h,o-z)
C---------------------------------------------------------------------
C...
      parameter (big0 = 10d+00)
C---------------------------------------------------------------------
      imin = 0
      if (dabs(a) .gt. 1.d-13 .and. dabs(b) .gt. 1.d-13) then
         ra = dabs(abs(a)/dabs(b)-1.d00)
         if(ra .lt. tol) then
            return
         else
            go to 10
         end if
      end if
C
      if(dabs(a) .gt. 1.d-13 .or. dabs(b) .gt. 1.d-13) go to 10
      return
C
 10   continue
      if(dabs(a) .gt. big0*dabs(b)) then
         imin = 2
      else if(dabs(b) .gt. big0*dabs(a)) then
         imin = 1
      else
         imin = 0
      end if
C
      return
      end
C====================================================================
C====================================================================
      subroutine shift(nxadj,nadj,n)
C====================================================================
      integer nxadj(*),nadj(*),n
C---------------------------------------------------------------------
C...
C---------------------------------------------------------------------
      kend = nxadj(n+1) - 1
      istrt = nxadj(1)
      do k = 1, n
         iend = nxadj(k+1) - 1
         klngth = iend - istrt + 1
         do j = istrt,iend
            if(nadj(j) .eq. 0) klngth = klngth - 1
         end do
         nxadj(k+1) = nxadj(k) + klngth
         istrt = iend + 1
      end do
      l = 0
      do k = 1, kend
         if (nadj(k) .ne. 0) then
            l = l + 1
            nadj(l) = nadj(k)
         end if
      end do
      return
      end
C======================================================================
      subroutine dfs(
     I     n,ia,ja,
     O     nblk,iblk,jblk,
     W     lowlink,iedge,numb)
C======================================================================
      integer ia(*),ja(*),jblk(*),iblk(*)
      integer n,n1,nb,nblk,count,k1,sp,vp,v,wp,w,v1,i
      integer lowlink(*),iedge(*),numb(*)
C---------------------------------------------------------------------
C...  I--input, O--output, W--working
C...  This subroutine performs the Tarjan's algorithm for finding the
C...  strongly connected components of a diggraph.  This is a
C...  modification of the FRED GUSTAVSON'S implementation of the
C...  algorithm.  There is no warranty that such version would work
C...  properly.  The new numbering is in CSR format the pointers are in
C...  IBLK, and the actual numbering is in JBLK.
C...  
C---------------------------------------------------------------------
C...  Initialization.
C
      nblk = 0
      nb = 0
C
      count  = nb+1
      k1 = count
      n1 = n+1
      vp = n1
      sp = n1
C
      do i = 1 , n
         iedge(i) = ia(i)
         numb(i) = 0
         lowlink(i) = 0
      end do
C
      numb(n1) = 0
C
C...  Get out when the  renumbering is done;
C...  Otherwise begin the at vertex K1 
C...
C
 10   if(count .eq. n1) then
C
         iblk(nblk+1) = n+1
         return
      end if
C
      do i = k1,n
         if(numb(i) .eq. 0) go to 30
      end do
      write(*,*) ' There is an error in DEPTH FIRST SEARCH.'
      stop 1
C
 30   continue
      v = i
      k1 = v + 1
      go to 50
C     
C...  :::
 40   continue
      vp = vp - 1
      iblk(vp) = v
      v=w
C
C...  
 50   continue
      nb = nb + 1
      numb(v) = nb
      lowlink(v) = numb(v)
      sp = sp - 1
      jblk(sp) = v
      v1 = v+1
C
C...  
 60   continue
      wp = iedge(v)
      w = ja(wp)
C
      iedge(v) = wp+1
C
      if(numb(w) .ge. numb(v)) go to 70
      if(numb(w) .eq. 0) go to 40
      lowlink(v) = min0(lowlink(v),numb(w))
 70   continue
      if(iedge(v) .lt. ia(v1)) go to 60
C
C...  
C  
      if(lowlink(v) .lt. numb(v)) go to 90
C
      nblk = nblk + 1
      iblk(nblk) = count
C
 80   continue
      w = jblk(sp)
      numb(w) = n1
      sp = sp + 1
      jblk(count) = w
      count = count + 1
      if(v .ne. w) go to 80
C
C...  
C
      if(sp .eq. n1) go to 10
 90   continue
      w = v
      v = iblk(vp)
      vp = vp + 1
      v1 = v + 1
      lowlink(v) = min0(lowlink(v),lowlink(w))
      go to 70
C
      end
C====================================================================
      subroutine permat(iord,ia,ja,an,n,m,iat,jat,ant)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*),ia(*),ja(*),iat(*),jat(*),n,m
      dimension an(*),ant(*)
C---------------------------------------------------------------------
C...  Permutes the matrix (P^t*A*P)
C---------------------------------------------------------------------
cc      call outmat(ia,ja,an,n)
cc      write(*,*) ' iord :' 
cc      write(*,*) (iord(k),k=1,n)
cc      read(*,*)
C
      call perm0(iord,ia,ja,an,n,m,iat,jat,ant)
C
cc      call aat(iat,jat,ant,n,m,ia,ja,an)
cc      call outmat(ia,ja,an,n)
cc      read(*,*)
C
      call perm0(iord,iat,jat,ant,m,n,ia,ja,an)
C
      return
      end
C====================================================================
      subroutine pervec(iord,u1,u2,n)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*)
      dimension u1(*),u2(*)
C---------------------------------------------------------------------
C...  Permutes the vector (u1 = P^t*u1)
C---------------------------------------------------------------------
      do k = 1, n
         u2(k) = u1(k)
      end do
      do k = 1, n
         u1(k) = u2(iord(k))
      end do
      return
      end
C====================================================================
      subroutine perback(iord,u1,u2,n)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*)
      dimension u1(*),u2(*)
C---------------------------------------------------------------------
C...  Permutes to the original vector (u1 = P*u1)
C---------------------------------------------------------------------
      do k = 1, n
         u2(k) = u1(k)
      end do
      do k = 1 , n
         u1(iord(k)) = u2(k)
      end do
      return
      end
C====================================================================
      subroutine perm0(iord,ia,ja,an,n,m,iat,jat,ant)
C====================================================================
      implicit real*8(a-h,o-z)
      integer*4 iord(*),ia(*),ja(*),iat(*),jat(*),n,m
      dimension an(*),ant(*)
C--------------------------------------------------------------------
C...  Permutation of the columns of a general sparse matrix A.
C...
C...  Input:
C...    IA, JA, AN    - given matrix A in RRCU.
C...    IORD          - the permutation. 
C...    N             - number of rows of the matrix.
C...    M             - number of columns of the matrix.
C...
C...  Output:
C...    IAT, JAT, ANT - transposed matrix A^t in RRCO.
C...
C...  Note:
C...    N+1 is the dimension of IA.
C...    M+1 is the dimension of IAT.
C--------------------------------------------------------------------
      mh = m + 1
      nh = n + 1
      do 10 i = 2, mh
         iat(i) = 0
 10   continue
      iab = ia(nh) - 1
      do 20 i = 1, iab
         j = ja(i) + 2
         if(j .le. mh) iat(j) = iat(j) + 1
 20   continue
      iat(1) = 1
      iat(2) = 1
      if (m .ne. 1) then
         do 30 i = 3, mh
            iat(i) = iat(i) + iat(i-1)
 30      continue
      end if
      do 50 i = 1, n
         iaa = ia(iord(i))
         iab = ia(iord(i)+1) - 1
         if(iab .lt. iaa) go to 50
         do 40 jp = iaa, iab
            j = ja(jp) + 1
            k = iat(j) 
            jat(k) = (i)
            ant(k) = an(jp)
            iat(j) = k + 1
 40      continue
 50   continue
C
      return
      end
C====================================================================
      subroutine icopyv(iu,iv,n)
C====================================================================
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension iu(*),iv(*)
C--------------------------------------------------------------------
C...  IV <--- IU : integer
C--------------------------------------------------------------------
      do i = 1 , n
         iv(i) = iu(i)
      end do
      return
      end
C=================================================================
      subroutine mxfrm2(n,ia,ja,nblk,iblock,jblock,
     >     mask,maxa,memt,maxbs)
      implicit real*8(a-h,o-z), integer (i-n)
      dimension mask(*),ia(*),ja(*),maxa(*),
     >     iblock(*),jblock(*)
C
C...  Goes over the blocks and forms the integer pointers for envelope
C.... (or skyline) format
C...  In another word computes the memory
C...  needed to store the decomposed blocks.  The total number of
C.... integers will be IBLOCK(NBLK+1)-1 + NBLK
C...  The total number of reals is MEMT
C...  The max block size is computed in MAXBS
      maxbs = 0
      iii = 0
      memt = 0
      i0 = 0
C..   memt on output is the total storage requirement for all the
C..   decomposed blocks
      do is = 1 , nblk
         ibl0 = iblock(is)
         ibl1 = iblock(is+1)-1
         nloc = ibl1-ibl0 + 1
         maxbs = max0(maxbs,nloc)
C... mask all entries in the block with the local numbering. 
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = i
         end do
         do i = 1 , nloc
            iblk = ibl0 + i - 1
c            write(*,*) iblk,i, ibl0
            ki = jblock(iblk)
            iaa = ia(ki) 
            iab = ia(ki+1)-1
            ihght = 0
            if (iaa .le. iab) then
               do kij = iaa, iab
                  kj = ja(kij)
                  j = mask(kj)
                  if(j .ne. 0 .and. j .le. i) then
CCCCCCCCCCCCCCCCCCthis is the i,j entry
                     ihght = max0(ihght,i-j)
                  end if
               end do
            end if
 !           write(*,*) iii,is,n
            maxa(iii+i+1) = ihght
         end do
         maxa(iii+1) = 1
         do i = 1 , nloc
            i0 = iii + i
            maxa(i0+1) = maxa(i0)+maxa(i0+1)+1
         end do
         memt = memt + maxa(i0+1)-1
C...  compute the address of the next block. 
         iii = iii + nloc + 1
C...  zero the mask so that everyting is as it was
         do i = 1, nloc
            iblk = ibl0 + i - 1
cc            write(*,*) iblk, ibl0, i
            ki = jblock(iblk)
            mask(ki) = 0
         end do
      end do
C...  
c      write(*,*) memt,maxbs
c      read(*,*) iii
      return
      end
C=================================================================
      subroutine sky2ns(n,ia,ja,a,nblk,iblock,jblock,mask,maxa,au,al)
      implicit real*8(a-h,o-z)
      parameter (zero=0d+0, one=1d+0)
      dimension iblock(*),jblock(*),mask(*),ia(*),ja(*),a(*)
      dimension maxa(*),au(*),al(*)
C
C...  Extracts block from ia,ja,a, and puts it in a skyline format (or
C...  envelope format). As input the integer pointer arrays should be
C...  stored in "maxa". Does an LU decomposition on each block and
C...  stores the lower triangles in au(). Each row of the lower triangle
C...  is stored between m1(ib)+maxa(i)--m1(ib)+maxa(i+1)-1, where m1 is
C...  the appropriate shift to get to the place relevant to the block
C...  "ib" m1(ib) is a shift, that is computed for each block -- it is
C...  not an array and is computed on fly.
      iii = 0
      m1 = 0
      i0 = 0
      do is = 1 , nblk
         ibl0 = iblock(is)
         ibl1 = iblock(is+1)-1
         nloc = ibl1-ibl0 + 1
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = i
         end do
C...  here (i,i), should be at position au(maxa(i))
C...  (i,i-1) is at au(maxa(i) + 1)
C...  (i,i-k) is at au(maxa(i) + k), for k < i
         do i = 1 , nloc
            i0 = iii + i
            ih0=maxa(i0) + m1 
            ih1=maxa(i0+1)-1 + m1
            do j0 = ih0,ih1
c     write(*,*) j0
               au(j0) = zero
               al(j0) = zero
            end do
         end do
C
         do i = 1 , nloc
            i0 = iii + i
C
C... Compute the addresses in the big array.
C
            ih0=maxa(i0) + m1 
            ih1=maxa(i0+1)-1 + m1
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            iaa = ia(ki) 
            iab = ia(ki+1)-1
            if (iaa .le. iab) then
               do kij = iaa, iab
                  kj = ja(kij)
                  j = mask(kj)
                  if(j .ne. 0 .and. i .le. j) then
                     j0   = iii + j
                     jh0  = maxa(j0) + m1
                     kdst = j-i
                     au(jh0 + kdst) = a(kij)
ccccc                     write(*,*) i,j,a(kij),jh0+kdst,'u',nloc,n
                  end if
                  if(j .ne. 0 .and. j .le. i) then
                     kdst = i-j
                     al(ih0 + kdst) = a(kij)
c                     write(*,*) i,j,a(kij), ih0+kdst,'l'
                  end if
               end do
            end if
         end do
C...  The first address of the next maxa. 
C... Now decompose the thing
C         
         i0 = iii + 1
         ih0 = m1 + 1
         ih1 = m1 + maxa(i0+nloc)-1
cc         write(*,*) 'm1,maxa(i0), maxa(i0+nloc)-1'
cc         write(*,*) m1,maxa(i0), maxa(i0+nloc)-1
C         write(*,*) 'MAXA=[',
C     >        (maxa(kk),kk=i0,i0+nloc),'];'
C         write(*,*) 'AU=[',
C     >        (au(ih0+kk),kk=0,maxa(i0+nloc)-2),'];'
C         write(*,*) 'AL=[',
C     >        (al(ih0+kk),kk=0,maxa(i0+nloc)-2),'];'
C         read(*,*) kk
         call doluns(au(ih0),al(ih0),maxa(i0),nloc)
cc         call dolu(au(ih0),maxa(i0),nloc)
CC         write(*,*) 'AU:',(au(ih0+kk),kk=0,maxa(i0+nloc)-2)
CC         read(*,*) kk
cc         write(*,*) 'AL:',(al(ih0+kk),kk=0,maxa(i0+nloc)-2)
cc         read(*,*) kk
C
         iii = iii + nloc + 1
C...  zero the mask so that everyting is as it was
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = 0
         end do
C... shift to the beginning of the next block.
         m1 = ih1
      end do
C...
      return
      end
C====================================================================
C UP TO HERE
C===================================================================
      subroutine fbgs2ns(n,ia,ja,a,x,b,
     >     nblk,iblock,jblock,mask,maxa,au,al,rhsloc,memt)
C===================================================================
      implicit real*8(a-h,o-z)
      parameter (zero=0d+0, one=1d+0)
      dimension ia(*),ja(*),a(*),x(*),b(*),
     >     iblock(*),jblock(*),mask(*),rhsloc(*)
      dimension au(*),al(*),maxa(*)
C
C... Block Gauss-Seidel, a.k.a. multiplicative Schwarz iteration
C... Blocks are in iblock, jblock
C... Decomposed diagonal blocks are in au()
C.... Mask must be 0 on entry
      iii = 0
      m1 = 0
      i0 = 0
      do is = 1 , nblk
         ibl0 = iblock(is)
         ibl1 = iblock(is+1)-1
         nloc = ibl1-ibl0 + 1
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = i
            rhsloc(i) = b(ki)
         end do
         do i = 1 , nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            iaa = ia(ki) 
            iab = ia(ki+1)-1
            if (iaa .le. iab) then
               do kij = iaa, iab
                  kj = ja(kij)
                  j = mask(kj)
                  if(j .eq. 0) 
     >                 rhsloc(i) = rhsloc(i) - a(kij)*x(kj)
               end do
            end if
         end do
         i0 = iii + 1
         ih0=m1 + 1
         ih1=maxa(i0 + nloc) - 1 + m1
         call sluns(au(ih0),al(ih0),rhsloc,maxa(i0),nloc)
cccccccccccccccccc         call slvlu(au(ih0),rhsloc,maxa(i0),nloc)
         iii = iii + nloc + 1
C...  zero the mask so that everyting is as it was
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = 0
            x(ki) = rhsloc(i)
         end do
C... shift to the beginning of the next block.
         m1 = ih1
      end do
C...  DONE; Here m1 should be = memt and iii = iblock(nblk+1)-1+nblk
      return
      end
C===================================================================
      subroutine bbgs2ns(n,ia,ja,a,x,b,
     >     nblk,iblock,jblock,mask,maxa,au,al,rhsloc,memt)
C===================================================================
      implicit real*8(a-h,o-z)
      parameter (zero=0d+0, one=1d+0)
      dimension ia(*),ja(*),a(*),x(*),b(*),
     >     iblock(*),jblock(*),mask(*),rhsloc(*)
      dimension au(*),al(*),maxa(*)
C
C...  backward Block Gauss-Seidel, a.k.a. multiplicative Schwarz
C...  iteration Blocks are in iblock, jblock Decomposed diagonal blocks
C...  are in au(). The difference with the forward one is only in
C...  computing the addresses for the L and U factors 
C.... Mask must be 0 on entry and is 0 on exit
      ibl0 = iblock(nblk)
      ibl1 = iblock(nblk+1)-1

      iii = ibl1 + nblk
      m1 = memt
      i0 = 0
      do is = nblk,1,-1
         ibl0 = iblock(is)
         ibl1 = iblock(is+1)-1
         nloc = ibl1-ibl0 + 1
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = i
            rhsloc(i) = b(ki)
         end do
         do i = 1 , nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            iaa = ia(ki) 
            iab = ia(ki+1)-1
            if (iaa .le. iab) then
               do kij = iaa, iab
                  kj = ja(kij)
                  j = mask(kj)
                  if(j .eq. 0)
     >                 rhsloc(i) = rhsloc(i) - a(kij)*x(kj)
               end do
            end if
         end do
         i0 = iii - nloc
         nzloc = maxa(iii) - 1
         ih0 = m1 - nzloc + 1
         ih1 = ih0 - 1 
cc         write(*,*) memt,maxa(i0),ih0,nloc,nzloc
         call sluns(au(ih0),al(ih0),rhsloc,maxa(i0),nloc)
         iii = iii - nloc - 1
C...  zero the mask so that everyting is as it was
         do i = 1, nloc
            iblk = ibl0 + i - 1
            ki = jblock(iblk)
            mask(ki) = 0
            x(ki) = rhsloc(i)
         end do
C... shift to the beginning of the next block.
         m1 = ih1
      end do
C... DONE
      return
      end
C=====================================================================
      subroutine doluns(au,al,maxa,nn)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension au(1),al(1),maxa(1)
C---------------------------------------------------------------------
C...  To calculate (l)*(u) factorization of nonsymmetric
C...  stiffness matrix using upper and lower envelopes
C---------------------------------------------------------------------
      if (nn.eq.1) return
      nlow=0
      ish = 0
      do 100 n=1,nn
      kn=maxa(n)
      ku=maxa(n+1)-1
      kh=ku-kn
C
      al(kn) = 1.d0
      if (kh) 100,100,110
C
 110  k = n - kh
      ki=maxa(k)
      al(ku) = al(ku)/au(ki)
cc    au(ku) = au(ku) !!!!!!
      kh1 = kh - 1
C
      if (kh1) 120,120,130
 130  k=n-kh1
      ic=0
      klt=ku
      do 140 j=1,kh1
         ic=ic+1
         klt=klt-1
         ki=maxa(k)
         nd=maxa(k+1)-ki-1
C
         if (nd) 141,145,145
 145     kk=min(ic,nd)
         cl=0.0d00
         cu=0.0d00
         do 150 l=1,kk
            kind = ki  + l
            lind = klt + l
            cl=cl + au(kind)*al(lind)
            cu=cu + al(kind)*au(lind)
 150     continue
         al(klt)=(al(klt)-cl)/au(ki)
         au(klt)= au(klt)-cu
 141     k=k+1
 140  continue
C
 120  b=0.d0
      kl = kn + 1
      do 160 kk=kl,ku
         cu=au(kk)
         cl=al(kk)
         b=b+cl*cu
 160  continue
      au(kn)=au(kn) - b
      if(dabs(au(kn)) .lt. 1.d-13) then
         ish = ish + 1
      end if
 100  continue
C
      if(ish .ne. 0) then
         write(*,'(a,i10,a)') 'there were ',ish,
     >        'NEGATIVE diagonal elements in the U factor'
      end if
 1111 format(//10x,'          there were  '/
     >     10x,'***Zero/negative diagonal entry in stiffnes matrix***'/
     >     /10x,'      i = ',i10,/10x,'     a(i,i) =',e15.7/)
      return
      end
C=====================================================================
      subroutine sluns(au,al,v,maxa,nn)
C=====================================================================
      implicit real*8 (a-h,o-z)
      dimension au(1),al(1),v(1),maxa(1)
C---------------------------------------------------------------------
C...  To reduce and back-substitute iteration vectors
C...  for nonsymmetrical matrix  
C---------------------------------------------------------------------
      do 400 n = 1 , nn
         kl = maxa(n) + 1
         ku = maxa(n+1) - 1
         if (ku - kl) 400,410,410
 410     k = n
         c = 0.d0
         do 420 kk = kl , ku
            k = k - 1
            c = c + al(kk) * v(k)
 420     continue
         v(n) = v(n) - c
 400  continue
C
cc      read (mt1) (a(j),j=1,nwk)
C
      k = maxa(nn)
      v(nn) = v(nn) / au(k)
C
      n = nn
      do 500 l = 2 , nn
         kl = maxa(n) + 1
         ku = maxa(n+1) - 1
         if (ku - kl) 530,510,510
 510     k = n
         do 520 kk = kl , ku
            k = k - 1
            v(k) = v(k) - au(kk) * v(n)
 520     continue
 530     nl = n - 1
         ln = maxa(nl)
         v(nl) = v(nl) / au(ln)
         n = n - 1
 500  continue
C
      return
      end
C=====================================================================
      subroutine dolu(a,maxa,nn)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension a(*),maxa(*)
C---------------------------------------------------------------------
C...  To calculate (l)*(d)*(l)(i) factorization of a matrix
C---------------------------------------------------------------------
      if (nn.eq.1) return
      nlow=0
      do 200 n=1,nn
         kn=maxa(n)
         kl=kn+1
         ku=maxa(n+1)-1
         kh=ku-kl
C
         if (kh) 304,240,210
 210     k=n-kh
         ic=0
         klt=ku
         do 260 j=1,kh
            ic=ic+1
            klt=klt-1
            ki=maxa(k)
            nd=maxa(k+1)-ki-1
            if (nd) 260,270,270
 270        kk=min0(ic,nd)
            c=0.0d00
            do 280 l=1,kk
               c=c+a(ki+l)*a(klt+l)
 280        continue
            a(klt)=a(klt)-c
            k=k+1
 260     continue
C
 240     k=n
         b=0.0d00
         do 300 kk=kl,ku
            k=k-1
            ki=maxa(k)
            c=a(kk)/a(ki)
c            if (dabs(c).lt.1.e07) go to 290
c            stop
 290        b=b+c*a(kk)
            a(kk)=c
 300     continue
         a(kn)=a(kn)-b
C
 304     if (a(kn)) 310,310,200
 310     nlow=nlow+1
         if (a(kn).eq.0) a(kn)=-1.e-16
         go to 200
 320     write(*,2000) n,a(kn)
         stop
 200  continue
C
 2000 format(//'   stop - matrix not positive definite'  ,//
     >         '   nonpositive diag element for equation' ,
     >     i4,//)
      return
      end
C=====================================================================
      subroutine slvlu(a,v,maxa,nn)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension a(*),v(*),maxa(1)
C---------------------------------------------------------------------
C... Solve with LUx = v. RHS is in V result is also in V.
C---------------------------------------------------------------------
      do 400 n=1,nn
         kl=maxa(n)+1
         ku=maxa(n+1)-1
         if (ku-kl) 400,410,410
 410     k=n
         c=0.
         do 420 kk=kl,ku
            k=k-1
            c=c+a(kk)*v(k)
 420     continue
         v(n)=v(n) - c
 400  continue
C
      do 480 n=1,nn
         k=maxa(n)
         v(n)=v(n)/a(k)
 480  continue
      if (nn.eq.1) return
      n=nn
      do  l=2,nn
         kl=maxa(n)+1
         ku=maxa(n+1)-1
         if (ku-kl) 530,510,510
 510     k=n
         do 520 kk=kl,ku
            k=k-1
            v(k)=v(k)-a(kk)*v(n)
 520     continue
 530     n=n-1
      end do
C
      return
      end
C====================================================================
      subroutine ijacrs(ln,ia,ja,a,n,nnz,ir,ic,aij)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(*),ja(*),ir(*),ic(*)
      dimension a(*),aij(*)
C--------------------------------------------------------------------
C (c) Hwanho, ltz 1999
C...  This subroutine converts the structures of a given matrix from
C...  (IR,IC,AIJ) to (IA,JA,A).
C...
C...  Input:
C...    IR, IC - the row and column indices of a given matrix.
C...    AIJ    - numerical values of the matrix corresponding to IR
C...             and IC.
C...  Output:
C...    IA, JA - structure of matrix A in RRCU. The order of A is N.
C...    A      - numerical values of nonzeros of A in RRCU.
C...    N      - the size of the matrix is N by N.
C...    NNZ    - number of nonzeros of A.
C--------------------------------------------------------------------
C
      ni   = 0
      nj   = 0
      nnz  = 0
C
c      write(*,'(e28.18)') (aij(k), k=1,5)
      do k = 1 , ln
         ni     = max(ni,ir(k))
         nj     = max(nj,ic(k))
      end do
C
      k = ln
      n = ni
      if (ni .ne. nj) then
         write(*,*) '**Warning: of rows and columns do not match:',
     >        ni, nj
         write(*,*)
      end if
C
      nnz = k
      do k = 1 , n+1
         ia(k) = 0
      end do
      do 30 k = 1, nnz
         irk = ir(k) + 1
         ia(irk) = ia(irk) + 1
 30   continue
C
      nzk = 1
      do 40 k = 2, n+1
         if (ia(k) .ne. 0) then
            nzk = k - 1
            go to 50
         end if
 40   continue
C
 50   continue
C
      do 60 k = 1, nzk
         ia(k) = 1
 60   continue
C
      do 70 k = nzk+1, n+1
         ia(k) = ia(k-1) + ia(k)
 70   continue
C
      do 80 k = 1, n
         ica = ia(k)
         icb = ia(k+1)
         if (icb .gt. ica)  ja(ia(k+1)-1) = ia(k)
 80   continue
C
      do 90 k = 1, nnz
         ik     = ir(k)
         iend   = ia(ik+1) - 1
         jp     = ja(iend)
         ja(jp) = ic(k)
         a(jp)  = aij(k)
         if (iend .ne. jp)  ja(iend) = ja(iend) + 1
 90   continue
      do i = 1 , n
         id = 0
         iai = ia(i)
         do jk = iai,ia(i+1)-1
            j = ja(jk)
            if(j .eq. i) then
               id = jk
               go to 100
            end if
         end do
         write(*,*) ' No diagonal element found in row :', i
         stop
 100     continue
         iw = ja(iai)
         ja(id) = iw
         ja(iai) = i
         asave = a(iai)
         a(iai) = a(id)
         a(id) = asave
      end do
C
      return
      end
C====================================================================
      subroutine sympat(ln,ia,ja,n,ir,ic,aij)
C====================================================================
      implicit real*8 (a-h,o-z)
      dimension ia(*),ja(*),ir(*),ic(*)
      dimension aij(*)
C--------------------------------------------------------------------
C ltz (2011) -- symmetrize the pattern of a matrix. Uses both csr and 
C     (i,j,v) structures
C--------------------------------------------------------------------
C
      lnnew=ln 
      do 30 i = 1, n
         iaa = ia(i)
         iab = ia(i+1)-1
         do 20 ii = iaa,iab
            j = ja(ii)
            if (i .eq. j) go to 20
            iac = ia(j)
            iad = ia(j+1)-1
            do 10 iik = iac, iad
               iii = iik
               if (ja(iii) .eq. i) go to 20
 10         continue
cccccc            write(*,*) ' Error j has no i, (i,j) = (',i,j
            lnnew=lnnew+1
            ir(lnnew)=j
            ic(lnnew)=i
            aij(lnnew)=0d0
 20      continue
 30   continue
C
      write (*,*) ' old-ln=',ln, 'new-ln=',lnnew
C      read(*,*) kk
      ln=lnnew
C
      return
      end
C====================================================================
      subroutine levels(inroot,ia,ja,mask,nlvl,iblock,jblock,maxlev)

C... Ludmil Zikatanov (1995); Following: Alan George, Joseph
C... C Liu. Computer Solution of Large Sparse Positive Definite Systems.
C... C Prentice Hall, 1981,
C------------------------------------------------------------------
      implicit real*8(a-h,o-z), integer(i-n)
      dimension ia(*),ja(*),jblock(*),mask(*),iblock(*)
c     
C If mask is non zero than we are at avisited node (not exactly as in the
C reference). Also an additional parameter maxlev is introduced
      if(mask(inroot) .ne. 0) return
C... This is a diagonal. 
      if(ia(inroot+1)-ia(inroot) .le. 1) then
         nlvl=1
         iblock(nlvl)=1
         jblock(iblock(nlvl))=inroot
         iblock(nlvl+1)=2
         return
      end if
      nlvl = 0
      jblock(1) = inroot
      lvlend = 0
      nsize  = 1
C     Initialization
      mask(inroot) = 1
 200  continue
      lbegin = lvlend + 1
      lvlend = nsize
      nlvl = nlvl + 1
      iblock(nlvl) = lbegin
      do i = lbegin,lvlend
         node = jblock(i)
         jstrt = ia(node)
         jstop = ia(node+1) -1
         if(jstop .ge. jstrt) then
            do  j = jstrt, jstop
               nbr = ja(j)
               if((mask(nbr) .le. 0)) then
                  nsize = nsize + 1
                  jblock(nsize) = nbr
                  mask(nbr) = nlvl
               end if
            end do
         end if
      end do
c     
      lvsize = nsize - lvlend
      if(lvsize .gt. 0 .and. nlvl .lt. maxlev) go to 200
      iblock(nlvl+1) = lvlend + 1
      do i = 1 , nsize
         node = jblock(i)
         mask(node) = 0
      end do
      nsize = iblock(nlvl+1) - iblock(nlvl)
C C     write(*,*)
C C     write(*,*) ' level structure rooted at node:', inroot
C C     do i=1,nlvl
C C       write(*,*)' Level:', i 
C C       write(*,*)(jblock(jk),jk=iblock(i),iblock(i+1)-1)
C C     end do
C C     read(*,*)
      return
      end
