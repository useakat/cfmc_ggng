      subroutine etaorder(n,p, iflg)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK
      implicit none

      integer i
      integer n,iflg    
      real*8 p(0:3,n)

      real*8 eta(n)
      integer ip(n)

      real*8 peta
      external peta
  
      iflg = 0

      do i = 1,n
         eta(i) = peta(p(0,i))
      enddo

      call rsort_asc(n,eta,ip)
      call porder(n,p,ip)

      return
      end
