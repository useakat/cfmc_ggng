      subroutine ptorder_chk(n,p, flg)
c by Yoshitaro Takaesu -Aug/03/2011 @KEK
      implicit none

      integer i
      integer n,flg
      real*8 p(0:3,n)

      real*8 pt(n)

      real*8 ppt
      external ppt
    
      flg = 0

      do i = 1,n
         pt(i) = ppt(p(0,i))
      enddo

      do i = 1,n-1
         if( pt(i+1).lt.pt(i) ) then
            flg = 1
            return
         endif
      enddo
   
      end

      
