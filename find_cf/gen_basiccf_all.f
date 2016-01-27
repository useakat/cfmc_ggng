      subroutine gen_basiccf(next,cf)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 15 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer next
      integer cf(next,next-1)
C     LOCAL VARIABLES 
      integer i,j
      integer pos
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do i = 1,next-1
         pos = 3
         cf(1,i) = 1
         do j = 2,next
            if (j.ne.i+1) then
               cf(j,i) = pos
               pos = pos +1
            elseif (j.eq.i+1) then
               cf(j,i) = 2
            endif
         enddo
      enddo

      return
      end
