      subroutine bfunc_ini

      include 'cparam.inc'
      include 'pshk.inc'
      include 'bfunc.inc'
      include 'lhe.inc'

      do i = 1,next
         pmass(i) = 0d0
      enddo
      do i = 1,next
         M(1,i) = pmass(i)
         M(2,i) = M(1,i)**2
         PUP(5,i) = pmass(i)
      enddo

      return
      end
