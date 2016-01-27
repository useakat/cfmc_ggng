      real*8 function es4dot(P1,P2)
C     ****************************************************     
C     
C     The function for inner product of lorentz four-vectors.
C
C     Input:
C           real*8 P1,P2: Lorentz four-vector
C     Output:
C            P1*P2
C
C     By Yoshitaro Takaesu @KEK Jun. 26 2011
C     ****************************************************
      implicit none
C     
C     ARGUMENTS 
C     
      real*8 P1(4),P2(4)
C     ----------
C     BEGIN CODE
C     ----------
      es4dot = P1(1)*P2(1) -P1(2)*P2(2) -P1(3)*P2(3) -P1(4)*P2(4)

      return
      end
