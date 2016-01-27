      subroutine mom_perm(n,a,P)

      implicitnone

      integer n,i,j
      integer a(n)
      real*8 P(0:3,n),Ptemp(0:3,n)

      do i = 1,n
         do j = 0,3
            Ptemp(j,i) = P(j,i)
         enddo
      enddo

      do i = 1,n
         do j = 0,3
            P(j,i) = Ptemp(j,a(i))
         enddo
      enddo

      return
      end
