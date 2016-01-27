      SUBROUTINE SPINIT(SEED)

c      INCLUDE 'crndm.inc'

      REAL POL(4)
      integer SEED


c      IF(ISEED.NE.0) THEN
c
c         PRINT *,'RM48 initialized: ',ISEED

c         NRNDM1 = 0
c         NRNDM2 = 0

c         CALL RM48IN(ISEED,NRNDM1,NRNDM2)
         NANDM1 = rand(SEED)
         NANDM2 = rand(SEED+10)

c      END IF

C     Initialize TAUOLA

C      CALL INITDK
C      CALL DEXAY(-1,POL)

      RETURN
      END
