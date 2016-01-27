      real*8 function fvdot(p)
c by Yoshitaro Takaesu -Aug/03/2011 @KEK
      implicit none
      
      real*8 p(4)

      fvdot = p(1)**2 -p(2)**2 -p(3)**2 -p(4)**2 

      end
