      real*8 function esbeta(a,b)
C     ****************************************************     
C     
C     The function for calculating beta(a,b).
C
C     Input:
C           real*8 a,b
C     Output:
C     Amplitude squared and summed
C
C     By Yoshitaro Takaesu @KEK Jun. 26 2011
C     ****************************************************
      implicit none
C     
C     ARGUMENTS 
C     
      real*8 a,b
C
C     ----------
C     BEGIN CODE
C     ----------
      esbeta = dsqrt( 1d0 +a**2 +b**2 -2*( a +b +a*b ) )

      return
      end
