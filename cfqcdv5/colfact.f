      subroutine colfact(n,aa,bb, cf)

      implicit none
      
      integer n,i,nused(n),cfp,npos1,npos0,npos2
      integer cf,aa(n,2),bb(n,2)

      do i = 1,n
         nused(i) = 0
      enddo
      cfp = 0
      npos1 = 1

      npos0 = npos1
 100  continue 
c      do i = 1,n
c         if (aa(i,1).eq.npos1) then
            nused(aa(i,1)) = 1
            npos2 = aa(i,2)
c         endif
c      enddo
    
      do i = 1,n
         if (npos2.eq.bb(i,2)) then
            npos1 = bb(i,1)
         endif
      enddo
      if (npos1.eq.npos0) then
         cfp = cfp +1
         do i = 1,n
            if (nused(i).ne.1) then
               npos1 = i
               goto 200
            endif
         enddo
         goto 300
 200     npos0 = npos1
         goto 100
      endif
      goto 100

 300  cf = 3**cfp

      return
      end
