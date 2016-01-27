      subroutine make_get_amp2
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'opt_info.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      character*4 cproc
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      open(1,file='get_amp2.f',status='replace')
      write(1,*) "      subroutine get_amp2(Q,amp2)"
      write(1,*) "      implicit none"
      write(1,*) "      include 'cparam.inc'"
      write(1,*) "      real*8 Q(0:3,next),amp2"
      write(1,*) "      integer i"
      write(1,*) "      real*8 matele"
      write(1,*) "      "
      write(1,*) "      if (icf.eq.0) then"
      write(1,*) "      call smatrix_mg_per(Q,amp2)"


      return
      end
