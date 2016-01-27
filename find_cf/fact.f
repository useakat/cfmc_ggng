      integer function fact(n)
      implicit none
      integer n,i      
      fact = 1
      if (n.eq.0) then
         fact = 1
      else
         do i = 1, n
            fact = fact*i
         enddo
      endif   
     
      return
      end
