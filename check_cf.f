      subroutine check_cf(next,cf1,cf2,iflag)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Dec.8 2012
C     ****************************************************
      implicitnone
C     
C     ARGUMENTS
C     
      integer next
      integer cf1(next),cf2(next),iflag
C     
C     LOCAL VARIABLES 
C     
      integer i
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 0
      do i = 1,next
         if (cf1(i).ne.cf2(i)) return
      enddo
      iflag = 1
      
      return
      end
      
