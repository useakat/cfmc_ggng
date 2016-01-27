      subroutine isort_asc(n, a,m)
c by Yoshitaro Takaesu -Nov/18/2012 @KIAS

      implicit none
      
      integer i,j
      integer n
      integer m(n)
      integer a(n)

      integer mini,tmp
      integer min

      do i = 1,n
         m(i) = i
      enddo
      do i = 1,n-1
         mini = i
         min = a(i)
         do j = i+1,n
            if( a(j) < min ) then
               min = a(j)
               mini = j
            endif
         enddo
         if( mini.ne.i ) then
            a(mini) = a(i)
            a(i) = min
            tmp = m(mini)
            m(mini) = m(i)
            m(i) = tmp
         endif
      enddo

      end


      subroutine isort_des(n, a,m)
c by Yoshitaro Takaesu -Nov/18/2012 @KIAS

      implicit none
      
      integer i,j
      integer n
      integer m(n)
      integer a(n)

      integer maxi,tmp
      integer max

      do i = 1,n
         m(i) = i
      enddo
      do i = 1,n-1
         maxi = i
         max = a(i)
         do j = i+1,n
            if( a(j) > max ) then
               max = a(j)
               maxi = j
            endif
         enddo
         if( maxi.ne.i ) then
            a(maxi) = a(i)
            a(i) = max
            tmp = m(maxi)
            m(maxi) = m(i)
            m(i) = tmp
         endif
      enddo

      end


      subroutine rsort_asc(n, a,m)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK

      implicit none
      
      integer i,j
      integer n
      integer m(n)
      real*8 a(n)

      integer mini,tmp
      real*8 min

      do i = 1,n
         m(i) = i
      enddo
      do i = 1,n-1
         mini = i
         min = a(i)
         do j = i+1,n
            if( a(j) < min ) then
               min = a(j)
               mini = j
            endif
         enddo
         if( mini.ne.i ) then
            a(mini) = a(i)
            a(i) = min
            tmp = m(mini)
            m(mini) = m(i)
            m(i) = tmp
         endif
      enddo

      end


      subroutine rsort_des(n, a,m)
c by Yoshitaro Takaesu -Jul/28/2011 @KEK

      implicit none
      
      integer i,j
      integer n
      integer m(n)
      real*8 a(n)

      integer maxi,tmp
      real*8 max

      do i = 1,n
         m(i) = i
      enddo
      do i = 1,n-1
         maxi = i
         max = a(i)
         do j = i+1,n
            if( a(j) > max ) then
               max = a(j)
               maxi = j
            endif
         enddo
         if( maxi.ne.i ) then
            a(maxi) = a(i)
            a(i) = max
            tmp = m(maxi)
            m(maxi) = m(i)
            m(i) = tmp
         endif
      enddo

      end
