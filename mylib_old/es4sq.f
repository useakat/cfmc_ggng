      real*8 function es4sq(P1)
C     ****************************************************     
C     
C     The function for square of a lorentz four-vector.
C
C     Input:
C           real*8 P1: Lorentz four-vector
C     Output:
C            P1^2
C
C     By Yoshitaro Takaesu @KEK Jun. 26 2011
C     ****************************************************
      implicit none
C     
C     ARGUMENTS 
C     
      real*8 P1(4)
C
C     EXTERNAL FUNCTIONS
C
      real*8 es4dot
      external es4dot
C
C     ----------
C     BEGIN CODE
C     ----------
      es4sq = es4dot(P1,P1)

      return
      end
