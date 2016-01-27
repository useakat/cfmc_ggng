      subroutine gen_basiccf(ibcf,next,cf)
      implicitnone
      integer next,ibcf
      integer cf(next)
      integer i,pos

      cf(1) = 1
      pos = 3
      do i = 2,next
         if (i.ne.ibcf+1) then
            cf(i) = pos
            pos = pos +1
         else
            cf(i) = 2
         endif
      enddo
      
      return
      end
