      real*8 function peta(p)

      implicit none

      real*8 p(0:3)

      peta = 0.5d0*dlog((p(0) +p(3))/(p(0) -p(3)))

      return
      end
