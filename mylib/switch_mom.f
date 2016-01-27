      subroutine switch_mom(p1,p2)

      implicit none

      integer i
      real*8 p1(0:3),p2(0:3),ptmp(0:3)

      do i = 0,3
         ptmp(i) = p1(i)
         p1(i) = p2(i)
         p2(i) = ptmp(i)
      enddo

      return
      end
