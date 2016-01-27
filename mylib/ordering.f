      subroutine order_pt(n,P,Q)
      implicitnone

      integer n,order(n),used(n),imax,i,j
      real*8 P(0:3,n),Q(0:3,n),pt(n),maxpt
      real*8 ppt
      external ppt
      
      do i = 1,n
         pt(i) = ppt(P(0,i))
      enddo
      do i = 1,n
         order(i) = i
         used(i) = 0
      enddo

      do j = 1,n
         maxpt = 0d0
         do i = 1,n
            if ((used(i).eq.0).and.(pt(i).gt.maxpt)) then
               maxpt = pt(i)
               imax = i
            endif
         enddo
         order(j) = imax
         used(imax) = 1
      enddo

      do i = 1,n
         do j = 0,3
            Q(j,i) = P(j,order(i))
         enddo
      enddo

      return
      end
