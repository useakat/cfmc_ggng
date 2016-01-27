      subroutine nperm(np,n, a)

      implicitnone

      integer np
      integer n,a(np),flag,i

      do i = 1,np
         a(i) = i
      enddo
      if (n.gt.1) then 
         do i = 1,n-1
            call ipnext(a,np,flag)
         enddo
      endif
   
      return
      end
      
