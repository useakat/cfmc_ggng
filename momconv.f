      subroutine momconv(q)
c by Yoshitaro Takaesu JAN 28 2013 @KIAS
      implicit none
      include 'cparam.inc'
      integer i
      real*8 q(0:3,next)
      external transfer

      do i = 1,next
         call transfer(i, q(0,i))
      enddo

      return
      end

      subroutine momconv4(q)
c by Yoshitaro Takaesu -Aug/03/2011 @KEK
      implicit none
      
      include 'cparam.inc'

      real*8 q(0:3,next)
      REAL*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      external transfer

      call transfer(1, p1)
      call transfer(2, p2)
      call transfer(3, p3)
      call transfer(4, p4)
      call transfer4(q,p1,p2,p3,p4)

      end


      subroutine momconv5(q)
c by Yoshitaro Takaesu -Aug/03/2011 @KEK
      implicit none
      
      include 'cparam.inc'

      real*8 q(0:3,next)
      REAL*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      external transfer

      call transfer(1, p1)
      call transfer(2, p2)
      call transfer(3, p3)
      call transfer(4, p4)
      call transfer(5, p5)
      call transfer5(q,p1,p2,p3,p4,p5)

      end


      subroutine momconv6(q)
c by Yoshitaro Takaesu -Aug/03/2011 @KEK
      implicit none
      
      include 'cparam.inc'

      real*8 q(0:3,next)
      REAL*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      external transfer

      call transfer(1, p1)
      call transfer(2, p2)
      call transfer(3, p3)
      call transfer(4, p4)
      call transfer(5, p5)
      call transfer(6, p6)
      call transfer6(q,p1,p2,p3,p4,p5,p6)

      end

