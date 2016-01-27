      subroutine wavsort(n,order,W)
      implicitnone

      integer n,i,j
      integer order(n)
      complex*16 W(18,n),WW(18,n)
     
      do i = 1,n
         do j = 1,18
            WW(j,i) = W(j,i)
         enddo
      enddo

      do i = 1,n
         do j = 1,18
            W(j,i) = WW(j,order(i))
         enddo
      enddo

      return
      end
      
