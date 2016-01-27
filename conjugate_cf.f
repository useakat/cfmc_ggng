      subroutine conjugate_cf(n,a,b)
      implicitnone

      integer n
      integer i,a(n),b(n)

      do i = 1,n
         if (i.ne.n) then
            b(i) = a(i+1)
         else
            b(i) = a(1)
         endif
      enddo

      return
      end
