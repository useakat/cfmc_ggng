      real*8 function ppt(p)

      implicit none

      real*8 p(0:3)

      ppt = dsqrt( p(1)**2 +p(2)**2 )

      return
      end
