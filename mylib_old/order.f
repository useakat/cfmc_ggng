      subroutine rorder(n,a,m)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK
      implicit none
      
      integer i,j
      integer n
      integer m(n)
      real*8 a(n)

      real*8 aa(n)

      do i = 1,n
         aa(i) = a(i)
      enddo

      do i = 1,n
         a(i) = aa(m(i))
      enddo

      end


      subroutine porder(n,a,m)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK
      implicit none
      
      integer i,j
      integer n
      integer m(n)
      real*8 a(0:3,n)

      real*8 aa(0:3,n)

      do i = 1,n
         do j = 0,3
            aa(j,i) = a(j,i)
         enddo
      enddo

      do i = 1,n
         do j = 0,3
            a(j,i) = aa(j,m(i))
         enddo
      enddo

      end
