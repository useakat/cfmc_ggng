C     -------------------------
      SUBROUTINE SMATRIX_co1(ICF,P,ANS)
C     -------------------------
C     
C     Generated by MadGraph 5 v. 1.4.1, 2012-06-02
C     By the MadGraph Development Team
C     Please visit us at https://launchpad.net/madgraph5
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
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
      REAL*8 MATRIX_co1_gg4g
      INTEGER IHEL,IDEN, I,ICF
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
      DATA IDEN/6144/
C     ----------
C     BEGIN CODE
C     ----------
      NTRY=NTRY+1
      ANS = 0D0
      DO IHEL=1,NCOMB
        IF (GOODHEL(IHEL) .OR. NTRY .LT. 2) THEN
          T=MATRIX_co1_gg4g(ICF,P ,NHEL(1,IHEL))
          ANS=ANS+T
          IF (T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
            GOODHEL(IHEL)=.TRUE.
          ENDIF
        ENDIF
      ENDDO
      ANS=ANS/DBLE(IDEN)
      END

C     ------------------------------
      REAL*8 FUNCTION MATRIX_co1_gg4g(ICF,P,NHEL)
C     ------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external momenta P(0:3,NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NPERMS
      PARAMETER (NPERMS=120)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,PERM(NEXTERNAL),ICF
      COMPLEX*16 ZTEMP
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 ONEPERM_co1_gg4g
      EXTERNAL ONEPERM_co1_gg4g

c      ZTEMP = (0.D0,0.D0)
c      DO I=1,NPERMS
        CALL GETPERM_co1_gg4g(ICF,PERM)
        ZTEMP=ONEPERM_co1_gg4g(P,NHEL,PERM)
c        ZTEMP=ZTEMP+ONEPERM(P,NHEL,PERM)
c      ENDDO
      MATRIX_co1_gg4g=DBLE(ZTEMP)

      RETURN
      END

C     --------------------------------------
      COMPLEX*16 FUNCTION ONEPERM_co1_gg4g(P,NHEL,PM)
C     --------------------------------------
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NFLOWS
      PARAMETER (NFLOWS=1)
      INTEGER    NPERMS
      PARAMETER (NPERMS=1)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL)
      INTEGER PM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J
      COMPLEX*16 ZTEMP
      COMPLEX*16 JAMP(NFLOWS)
      INTEGER PERMS(NEXTERNAL,NPERMS),IFERM(NPERMS),PERM(NEXTERNAL)
      DATA (PERMS(I,   1),I=1,6) / 1, 2, 3, 4, 5, 6/
      DATA IFERM/ 1/
C     
C     EXTERNAL FUNCTIONS
C     
      COMPLEX*16 FLOW1
      EXTERNAL FLOW1
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=PM(PERMS(I,1))
      ENDDO
      JAMP(1)=IFERM(1)*FLOW1(P,NHEL,PERM)
      ZTEMP = (0.D0,0.D0)
      ZTEMP = ZTEMP+1/8D0*JAMP(1)*DCONJG(81D0*(JAMP(1)))
      ONEPERM_co1_gg4g=ZTEMP

      RETURN
      END

C     ------------------------------------
      SUBROUTINE GETPERM_co1_gg4g(IPERM,PERM)
C     ------------------------------------
C     
C     Gives permutation number IPERM. 
C     Return value is the fermion factor due to PERM
C     
C     Process: g g > g g g g QCD=4 QED=0 WEIGHTED=4 singlet_QCD=0
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
C     
C     ARGUMENTS 
C     
      INTEGER IPERM,PERM(NEXTERNAL)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,IFLAG
      LOGICAL OK
      INTEGER COMP(NEXTERNAL)
      DATA COMP/1,1,1,1,1,1/
C     ----------
C     BEGIN CODE
C     ----------
      DO I=1,NEXTERNAL
        PERM(I)=I
      ENDDO
      I=1
      DO WHILE(I.LT.IPERM)
        CALL IPNEXT(PERM,NEXTERNAL,IFLAG)
        OK=.TRUE.
        DO J=1,NEXTERNAL
          IF(COMP(PERM(J)).NE.COMP(J))THEN
            OK=.FALSE.
            EXIT
          ENDIF
        ENDDO
        IF(OK) I=I+1
      ENDDO
      END

