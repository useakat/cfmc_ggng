      subroutine coldel(np,a,q0, aa,den)

      implicit none

      integer i
      integer np,a(np),ng,flag,aa(np,2),pos0,q0
      real*8 den

      den = 1d0

      pos0 = 1
      do i = pos0,np

         aa(i,1) = a(i)
         if (i.eq.pos0) then
            aa(i,2) = q0
         else
            aa(i,2) = a(i-1)
         endif

         if (a(i).eq.q0) then
            pos0 = i+1
            goto 100
         endif

      enddo

 100  do i = pos0,np   
         aa(i,1) = a(i)
         aa(i,2) = a(i)
         den = den*3
      enddo
      
      return
      end
