      subroutine chcf(iicf)
      implicitnone
      include 'cparam.inc'
      integer iicf,pos,i

      floworder(1,1) = 1
      pos = 3
      do i = 2,next
         if (i.ne.iicf+1) then
            floworder(i,1) = pos
            pos = pos +1
         else
            floworder(i,1) = 2
         endif
      enddo

      return
      end
