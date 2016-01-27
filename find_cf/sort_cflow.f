      subroutine sort_cflow(atm,ia,natm,next,cflow)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C
C     
C     ****************************************************
      implicit none

C     GLOBAL VARIABLES

C     CONSTANTS

C     ARGUMENTS 
      integer natm,next
      integer ia(natm),cflow(next),atm(2,natm)
C     LOCAL VARIABLES 
      integer i,j,cfpos
C     EXTERNAL FUNCTIONS

C     ----------
C     BEGIN CODE
C     ----------
      cfpos = 0
      do i = 1,natm
         do j = 1,2
            if (atm(j,ia(i)).ne.0) then
               cfpos = cfpos +1
               cflow(cfpos) = atm(j,ia(i))
            endif
         enddo
      enddo

      return
      end
