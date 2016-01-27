      subroutine get_chjacob_diagram_6ch(P)
C     ****************************************************
C     
C     The driver for making include and dat files for QCD
C     processes in color-flow-sampling option of MG. 
C     Input:
C     pp    4 momentum of external particles
C     Output:
C     Amplitude squared and summed
C
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicitnone

      INCLUDE 'coupl.inc'
      include 'cparam.inc'
C     
C     CONSTANTS
C     
      INTEGER    NWAVEFUNCS
      PARAMETER (NWAVEFUNCS=80)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I
      COMPLEX*16 W(18,NWAVEFUNCS),AMP(27)
      real*8 P(0:3,next)
      integer NCOMB
      parameter (NCOMB=32)
      INTEGER NHEL(NEXT,NCOMB)
      INTEGER IHEL
      INTEGER IC(NEXT)
      DATA (NHEL(I,   1),I=1,5) /-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,5) /-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,5) /-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,5) /-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,5) /-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,5) /-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,5) /-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,5) /-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,5) /-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,5) /-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,5) /-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,5) /-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,5) /-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,5) /-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,5) /-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,5) /-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,5) / 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,5) / 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,5) / 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,5) / 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,5) / 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,5) / 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,5) / 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,5) / 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,5) / 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,5) / 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,5) / 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,5) / 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,5) / 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,5) / 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,5) / 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,5) / 1, 1, 1, 1, 1/
C
C     ----------
C     BEGIN CODE
C     ----------
      DO i=1,NEXT
        IC(i) = +1
      ENDDO
      do i=1,nwgt
         AA(i) = 0d0
      enddo
      DO IHEL=1,NCOMB
         CALL VXXXXX(P(0,1),ZERO,NHEL(1,IHEL),-1*IC(1),W(1,1))
         CALL VXXXXX(P(0,2),ZERO,NHEL(2,IHEL),-1*IC(2),W(1,2))
         CALL VXXXXX(P(0,3),ZERO,NHEL(3,IHEL),+1*IC(3),W(1,3))
         CALL VXXXXX(P(0,4),ZERO,NHEL(4,IHEL),+1*IC(4),W(1,4))
         CALL VXXXXX(P(0,5),ZERO,NHEL(5,IHEL),+1*IC(5),W(1,5))
         CALL VVV1_1(W(1,1),W(1,2),GC_4,ZERO, ZERO, W(1,6))
         CALL VVV1_1(W(1,3),W(1,4),GC_4,ZERO, ZERO, W(1,7))
         CALL VVV1_1(W(1,3),W(1,5),GC_4,ZERO, ZERO, W(1,8))
         CALL VVV1_1(W(1,4),W(1,5),GC_4,ZERO, ZERO, W(1,9))
         CALL VVV1_1(W(1,1),W(1,3),GC_4,ZERO, ZERO, W(1,10))
         CALL VVV1_1(W(1,2),W(1,4),GC_4,ZERO, ZERO, W(1,11))
C     Amplitude(s) for diagram number 1
         CALL VVV1_0(W(1,6),W(1,7),W(1,5),GC_4,AMP(1))
         CALL VVV1_1(W(1,3),W(1,5),GC_4,ZERO, ZERO, W(1,8))
C     Amplitude(s) for diagram number 2
         CALL VVV1_0(W(1,6),W(1,8),W(1,4),GC_4,AMP(2))
         CALL VVV1_1(W(1,4),W(1,5),GC_4,ZERO, ZERO, W(1,9))
C     Amplitude(s) for diagram number 3
         CALL VVV1_0(W(1,6),W(1,3),W(1,9),GC_4,AMP(3))
C     Amplitude(s) for diagram number 5
         CALL VVV1_0(W(1,10),W(1,11),W(1,5),GC_4,AMP(7))
         CALL VVV1_1(W(1,2),W(1,5),GC_4,ZERO, ZERO, W(1,12))
C     Amplitude(s) for diagram number 6
         CALL VVV1_0(W(1,10),W(1,12),W(1,4),GC_4,AMP(8))
C     Amplitude(s) for diagram number 7
         CALL VVV1_0(W(1,10),W(1,2),W(1,9),GC_4,AMP(9))
         CALL VVV1_1(W(1,1),W(1,4),GC_4,ZERO, ZERO, W(1,13))
         CALL VVV1_1(W(1,2),W(1,3),GC_4,ZERO, ZERO, W(1,14))
C     Amplitude(s) for diagram number 9
         CALL VVV1_0(W(1,13),W(1,14),W(1,5),GC_4,AMP(13))
C     Amplitude(s) for diagram number 10
         CALL VVV1_0(W(1,13),W(1,12),W(1,3),GC_4,AMP(14))
C     Amplitude(s) for diagram number 11
         CALL VVV1_0(W(1,13),W(1,2),W(1,8),GC_4,AMP(15))
         CALL VVV1_1(W(1,1),W(1,5),GC_4,ZERO, ZERO, W(1,15))
C     Amplitude(s) for diagram number 13
         CALL VVV1_0(W(1,15),W(1,14),W(1,4),GC_4,AMP(19))
C     Amplitude(s) for diagram number 14
         CALL VVV1_0(W(1,15),W(1,11),W(1,3),GC_4,AMP(20))
C     Amplitude(s) for diagram number 15
         CALL VVV1_0(W(1,15),W(1,2),W(1,7),GC_4,AMP(21))
C     Amplitude(s) for diagram number 17
         CALL VVV1_0(W(1,1),W(1,14),W(1,9),GC_4,AMP(25))
C     Amplitude(s) for diagram number 18
         CALL VVV1_0(W(1,1),W(1,11),W(1,8),GC_4,AMP(26))
C     Amplitude(s) for diagram number 19
         CALL VVV1_0(W(1,1),W(1,12),W(1,7),GC_4,AMP(27))
         AA(1) = AA(1) +AMP(8)*DCONJG(AMP(8))     ! 1: (3,4,5)
         AA(2) = AA(2) +AMP(7)*DCONJG(AMP(7))     ! 2: (3,5,4)
         AA(3) = AA(3) +AMP(14)*DCONJG(AMP(14))   ! 3: (4,3,5)
         AA(4) = AA(4) +AMP(13)*DCONJG(AMP(13))   ! 4: (4,5,3)
         AA(5) = AA(5) +AMP(20)*DCONJG(AMP(20))   ! 5: (5,3,4)
         AA(6) = AA(6) +AMP(19)*DCONJG(AMP(19))   ! 6: (5,4,3)
         AA(7) = AA(7) +AMP(9)*DCONJG(AMP(9))     ! 7: (3,(4,5))
         AA(8) = AA(8) +AMP(15)*DCONJG(AMP(15))   ! 8: (4,(3,5))
         AA(9) = AA(9) +AMP(21)*DCONJG(AMP(21))   ! 9: (5,(3,4))
         AA(10) = AA(10) +AMP(25)*DCONJG(AMP(25)) ! 10: ((4,5),3)
         AA(11) = AA(11) +AMP(26)*DCONJG(AMP(26)) ! 11: ((3,5),4)
         AA(12) = AA(12) +AMP(27)*DCONJG(AMP(27)) ! 12: ((3,4),5)
      ENDDO

      AAA = 0d0
      do i = 1,nwgt
         AAA = AAA +AA(i)
      enddo

      AA(1) = AA(1) +AA(6)
      AA(2) = AA(2) +AA(4)
      AA(3) = AA(3) +AA(5)
      AA(4) = AA(7) +AA(10)
      AA(5) = AA(8) +AA(11)
      AA(6) = AA(9) +AA(12)

      return
      end
