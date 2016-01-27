************************************************************************
*    ==============================                                    *
      SUBROUTINE XHFILL(ID, DX, FX )
*    ==============================                                    *
* ((Function))                                                         *
*     To fill histograms.                                              *
*   This routine identifies the bin number which is to be updated      *
*   with weight FX*WGT.  Up to five points per histogram are able      *
*   to be stacked before calling BHUPDT or SHUPDT.                     *
* ((Input))                                                            *
*   ID    : Histogram identification number                            *
*   DX    : Input value                                                *
*   FX    : Input value of the function                                *
* ((Author))                                                           *
*   S.Kawabata         June '90 at KEK                                 *
*                                                                      *
************************************************************************

      IMPLICIT NONE

      INTEGER ID
      DOUBLE PRECISION DX,FX

      INCLUDE 'base0.h'
      INCLUDE 'base3.h'

      INCLUDE 'plot.h'

      INTEGER IHIST,IX
      INTEGER IP1,IP2,IP3

      REAL X,XMIN,XMAX,DEV,FXWGT
      INTEGER NXBIN
      
      INTEGER I,K
      

      IF( NHIST .GT. 0 ) THEN

          I  = IABS(MOD( ID, 13 )) + 1
          IF( XHASH(1, I) .EQ. 1 ) THEN
            IF( ID .EQ. MAPL( 1, XHASH(2,I))) THEN
                IHIST = XHASH(2,I)
                GO TO 200
            ENDIF
          ELSEIF( XHASH(1, I) .GT. 1 ) THEN
            DO 100 K = 2, XHASH(1,I)+1
               IF( ID .EQ. MAPL( 1, XHASH(K,I))) THEN
                   IHIST = XHASH(K,I)
                   GO TO 200
               ENDIF
  100       CONTINUE
          ENDIF
      ENDIF
C     IF( LU .GT. 0 ) THEN
C         WRITE(LU,9000) ID
C     ENDIF
C9000 FORMAT(1X,'No Histogram corresponds to ID =',I5,
C    .      /1X,' This call is neglected.')
      RETURN
C

  200 X     = DX*1.0

          IX    = -1
          IP1   = MAPL(2,IHIST)
          XMIN  = BUFF(IP1)
          XMAX  = BUFF(IP1+1)
          NXBIN = IBUF(IP1+2)
          DEV   = BUFF(IP1+3)
          IF(     X .LT. XMIN ) THEN
                  IX   = 0
          ELSEIF( X .GT. XMAX ) THEN
                 IX   = NXBIN + 1
          ELSE
                 IX   = INT((X - XMIN)/DEV + 1.0)
                 IF( IX .GT. NXBIN ) IX = NXBIN
          ENDIF
C        PRINT*,'ID, IHIST, IFBASE =',ID,IHIST,(IFBASE(I),I=1,NHIST)

      IF( IBASES .EQ. 1 ) THEN

          IP2       = MAPL(3,IHIST) + IX
          IBUF(IP2) = IBUF(IP2) + 1
          FXWGT     = FX*WGT
          IP2       = IP2 + 52
          BUFF(IP2) = BUFF(IP2) + FXWGT
          IP2       = IP2 + 52
          BUFF(IP2) = BUFF(IP2) + FXWGT*FXWGT
*   Add March 1994
          IFBASE(IHIST) = 1

      ELSE
C        PRINT*,'ID, IHIST, IFBASE =',ID,IHIST,(IFBASE(I),I=1,NHIST)

         IP3        =  MAPL(4,IHIST)
         IBUF(IP3)  = IX

      ENDIF

C
      RETURN
      END
