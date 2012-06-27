!> \file   factor.f
!> \brief  LU factoraization for CSR matrix
!> \author Ludmil Zikatanov
!> \date   01/01/2002

C=====================================================================
      subroutine sfactr(ia,ja,n,iu,ju,ip,nwku)
C=====================================================================
      dimension ia(1),ja(1),iu(1),ju(1),ip(n)
C--------------------------------------------------------------------
C...  SYMBOLIC FACTORIZATION OF A SYMMETRIC SPARSE MATRIX.
C...  Gives the upper triagle (U) of the LU factorization;
C--------------------------------------------------------------------
      nm = n - 1
      nh = n + 1
      do i = 1 , n
         iu(i) = 0
         ip(i) = 0
      end do
C
      jp = 1
      do i = 1 , nm
         jpi = jp
         jpp = n + jp - i
         min = nh
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .ge. iaa) then
            do j = iaa,iab
               jj = ja(j)
               ju(jp) = jj
               jp = jp + 1
               if(jj .lt. min) min = jj
               iu(jj) = i
            end do
         end if
          last = ip(i)
         if(last .eq. 0) go to 60
         l = last
 40      l = ip(l)
         lh = l + 1
         iua = iu(l)
         iub = iu(lh) - 1
         if(lh .eq. i) iub = jpi - 1
         iu(i) = i
         do j = iua,iub
            jj = ju(j)
            if(iu(jj) .ne. i) then
               ju(jp) = jj
               jp = jp + 1
               iu(jj) = i
               if(jj .lt. min)min = jj
            end if
         end do
         if(jp .eq. jpp) go to 70
         if(l.ne. last) go to 40
 60      if (min .eq. nh) go to 90
 70      l = ip(min)
         if(l .eq. 0) go to 80
         ip(i) = ip(l)
         ip(l) = i
         go to 90
 80      ip(min) = i
         ip(i) = i
 90      continue
         iu(i) = jpi
      end do
      iu(n) = jp
      iu(nh) = jp
      nwku = iu(n+1)
      return
      end

C=====================================================================
      subroutine sfactr_new(ia,ja,n,iu,ju,ip,nwku,mem_chk)
C=====================================================================
      dimension ia(1),ja(1),iu(1),ju(1),ip(n)
C--------------------------------------------------------------------
C...  SYMBOLIC FACTORIZATION OF A SYMMETRIC SPARSE MATRIX.
C--------------------------------------------------------------------
      nm = n - 1
      nh = n + 1
      iu(nh) = 0
      do i = 1 , n
         iu(i) = 0
         ip(i) = 0
      end do
C
      jp = 1
      do i = 1 , nm
cc         if(mod(i,10) .eq. 1.and. nm .lt. 7000) write(*,*) i
         jpi = jp
         jpp = n + jp - i
         min = nh
         iaa = ia(i)
         iab = ia(i+1)-1
         if(iab .ge. iaa) then
            do j = iaa,iab
               jj = ja(j)
C...           memory check, 
               if(jp .le. mem_chk) then
                  ju(jp) = jj
               else
                  iu(n+1) = mem_chk+3
                  return
               end if
               jp = jp + 1
               if(jj .lt. min) min = jj
               iu(jj) = i
            end do
         end if
         last = ip(i)
         if(last .eq. 0) go to 60
         l = last
 40      l = ip(l)
         lh = l + 1
         iua = iu(l)
         iub = iu(lh) - 1
         if(lh .eq. i) iub = jpi - 1
         iu(i) = i
         do j = iua,iub
C...        memory check, 
            if(j .le. mem_chk) then
               jj = ju(j)
            else
               iu(n+1) = mem_chk+3
               return
            end if
            if(iu(jj) .ne. i) then
C...           memory check, 
               if(jp .le. mem_chk) then
                  ju(jp) = jj
               else
                  iu(n+1) = mem_chk+3
                  return
               end if
               jp = jp + 1
               iu(jj) = i
               if(jj .lt. min) min = jj
            end if
         end do
         if(jp .eq. jpp) go to 70
         if(l.ne. last) go to 40
 60      if (min .eq. nh) go to 90
 70      l = ip(min)
         if(l .eq. 0) go to 80
         ip(i) = ip(l)
         ip(l) = i
         go to 90
 80      ip(min) = i
         ip(i) = i
 90      continue
         iu(i) = jpi
      end do
      iu(n) = jp
      iu(nh) = jp
      nwku = iu(n+1)
      return
      end

C=====================================================================
      subroutine factor(ia,ja,n,iu,ju,ip,iup,an,ad,un,di)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension ia(1),ja(1),iu(1),ju(1),ip(1),iup(1)
      dimension an(1),ad(1),di(1),un(1)
C--------------------------------------------------------------------
C...  Factorization of symmetric sparse matrix.
C--------------------------------------------------------------------
      do 10 j = 1 , n
         iup(j) = 0
         ip(j) = 0
 10   continue
C
      do 130 i = 1 , n
cc       if(mod(i,10).eq.1 .and. n.lt.7000) write(*,*) ' factoring',i
         ih = i + 1
         iua = iu(i)
         iub = iu(ih) - 1
         if(iub .lt. iua) go to 40
         do 20 j = iua,iub
            di(ju(j)) = 0.0d00
 20      continue
         iaa = ia(i)
         iab = ia(ih) - 1
         if(iab .lt. iaa) go to 40
         do 30 j = iaa,iab
            di(ja(j)) = an(j)
 30      continue
 40      di(i) = ad(i)
         last = ip(i)
         if(last .eq. 0) go to 90
         ln = ip(last)
 50      l = ln
         ln = ip(l)
         iuc = iup(l)
         iud = iu(l + 1) -1
         um = un(iuc) * di(l)
         do 60 j = iuc,iud
            jj = ju(j)
            di(jj) = di(jj) - un(j) * um
 60      continue
         un(iuc) = um
         iup(l) = iuc + 1
         if(iuc .eq. iud) go to 80
         j = ju(iuc+1)
         jj = ip(j)
         if(jj .eq. 0) go to 70
         ip(l) = ip(jj)
         ip(jj) = l
         go to 80
 70      ip(j) = l
         ip(l) = l
 80      if(l .ne. last) go to 50
C
 90      di(i) = 1.d00/di(i)
         if(iub .lt. iua) go to 120
         do 100 j = iua,iub
            un(j) = di(ju(j))
 100     continue
         j = ju(iua)
         jj = ip(j)
         if(jj .eq. 0) go to 110
         ip(i) = ip(jj)
         ip(jj) = i
         go to 120
 110     ip(j) = i
         ip(i) = i
 120     iup(i) = iua
 130  continue
C
      return
      end

C=====================================================================
      subroutine forbac(iu,ju,un,di,n,x)
C=====================================================================
      implicit real*8(a-h,o-z)
      dimension iu(1),ju(1),un(1),x(1),di(1)
C--------------------------------------------------------------------
C...  REDUCTION AND BACK SUBSTITUTION FOR SOLVING THE
C...  SYMMETRIC SYSTEM WITH DECOMPOSED (FACTORIZED MATRIX)
C--------------------------------------------------------------------
      nm = n-1
c     do i = 1 , n
c        x(i) = b(i)
c     end do
      do k = 1 , nm
         iua = iu(k)
         iub = iu(k+1) - 1
         xx = x(k)
         if(iub .ge. iua) then
            do i = iua,iub
               x(ju(i)) = x(ju(i))-un(i)*xx
            end do
         end if
         x(k) = xx*di(k)
      end do
      x(n) = x(n)*di(n)
      k = nm
 50   iua = iu(k)
      iub = iu(k+1) - 1
      if(iub .lt. iua) go to 70
      xx = x(k)
      do i = iua,iub
         xx = xx - un(i) *x(ju(i))
      end do
      x(k) = xx
 70   k = k - 1
      if(k .gt. 0) go to 50
C
      return
      end
C===================================================================== 
