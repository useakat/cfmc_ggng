      SUBROUTINE SMATRIX_cfqcd(P,ANS)
C     
C     Generated by MadGraph 5 v. 1.3.9, 2011-08-01
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: g21 g32 > g43 g54 g65 g16 WEIGHTED=4
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
C     
C     LOCAL VARIABLES 
C     
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 MATRIX_cfqcd_gg4g
      INTEGER IHEL,IDEN, I
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./
      DATA (NHEL(I,   1),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   2),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   3),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   4),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   5),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   6),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   7),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   8),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   9),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  10),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  11),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  12),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  13),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  14),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  15),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  16),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  17),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  18),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  19),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  20),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  21),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  22),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  23),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  24),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  25),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  26),I=1,6) /-1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  27),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  28),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  29),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  30),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  31),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  32),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  33),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  34),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  35),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  36),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  37),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  38),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  39),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  40),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  41),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  42),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  43),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  44),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  45),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  46),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  47),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  48),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  49),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  50),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  51),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  52),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  53),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  54),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  55),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  56),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  57),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  58),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  59),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  60),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  61),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  62),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  63),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  64),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA IDEN/ 6144/
C     ----------
C     BEGIN CODE
C     ----------
      NTRY=NTRY+1
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      ANS = 0D0
      DO IHEL=1,NCOMB
        IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
          T=MATRIX_cfqcd_gg4g(P ,NHEL(1,IHEL),JC(1))
          ANS=ANS+T
          IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
            GOODHEL(IHEL)=.TRUE.
          ENDIF
        ENDIF
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END


      REAL*8 FUNCTION MATRIX_cfqcd_gg4g(P,NHEL,IC)
C     
C     Generated by MadGraph 5 v. 1.3.9, 2011-08-01
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g21 g32 > g43 g54 g65 g16 WEIGHTED=4
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=1)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=6, NCOLOR=1)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
      INCLUDE 'cparam.inc'
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,1),I=1,1) /1/
C     ----------
C     BEGIN CODE
C     ----------
      CALL VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL VXXXXX(P(0,3),ZERO,NHEL(3),+1*IC(3),W(1,3))
      CALL VXXXXX(P(0,4),ZERO,NHEL(4),+1*IC(4),W(1,4))
      CALL VXXXXX(P(0,5),ZERO,NHEL(5),+1*IC(5),W(1,5))
      CALL VXXXXX(P(0,6),ZERO,NHEL(6),+1*IC(6),W(1,6))
      CALL wavsort(next,floworder(1,1),W)
      CALL GLUON6_0(W(1,1),W(1,2),W(1,3),W(1,4),W(1,5),W(1,6),GC_4
     $ ,AMP(1))
      JAMP(1)=+AMP(1)

      MATRIX_cfqcd_gg4g = 3**(next-2)*(3**2-1d0)*JAMP(1)*DCONJG(JAMP(1))
     &     /(2d0**4)

      END
