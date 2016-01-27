C***********************************************************************
C*                                                                     *
C*========================                                             *
C*    SUBROUTINE BSETGU                                                *
C*========================                                             *
C*((Function))                                                         *
C*     Initialization of Bases progam                                  *
C*     This is called only when IFLAG=0.                               *
C*     ( IFLAG = 0 ; First Trial of Defining Grid step )               *
C*                                                                     *
C*    Changed by S.Kawabata    Aug. 1984 at Nagoya Univ.               *
C*    Last update              Oct. 1985 at KEK                        *
C*                                                                     *
C***********************************************************************
      SUBROUTINE BSETGU

      IMPLICIT NONE
C     IMPLICIT REAL*8 (A-H,O-Z)

      INCLUDE 'base1.h'
      INCLUDE 'base4.h'
      INCLUDE 'base6.h'

      INTEGER M,NDM
      DOUBLE PRECISION RC,XN,DR,XO
      DOUBLE PRECISION XIN(NDMX)
      INTEGER NSP,IPSAVE

      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0)

      INTEGER I,J,K

C
C---------------------------------------------------------------
C           Define the number of grids and sub-regions
C---------------------------------------------------------------
C==> Determine NG : Number of grids
      write(99,*) "BASES"
      write(99,*) "NCALL",NCALL
          NG    = (NCALL/2.)**(1./NWILD)
         IF(NG .GT. 25) NG  = 25
  100    IF(NG .LT.  2) NG  =  1
         IF(NG**NWILD .GT. LENG) THEN
            NG  = NG - 1
            GO TO 100
         ENDIF
C
C==> Determine ND : Number of sub-regions
          M     = NDMX/NG
          ND    = M*NG
C
C==> Determine NPG: Number of sampling points per subhypercube
          NSP   = NG**NWILD
          NPG   = NCALL/NSP

          XI(1,1)= ONE
          MA(1)  = 1
          DX(1)  = XU(1)-XL(1)

          IF( NDIM .GT. 1 ) THEN
             DO 130 J = 2,NDIM
                 XI(1,J)= ONE
                 DX(J)  = XU(J)-XL(J)
                 IF( J .LE. NWILD ) THEN
                    MA(J)  = NG*MA(J-1)
                 ENDIF
  130         CONTINUE
          ENDIF
C
C---------------------------------------------------------------
C           Set size of subregions uniform
C---------------------------------------------------------------
          NDM   = ND-1
          RC    = ONE/ND
          DO 155 J =1,NDIM
             K     = 0
             XN    = 0.D0
             DR    = XN
             I     = K
  140        K     = K+1
             DR    = DR+ONE
             XO    = XN
             XN    = XI(K,J)
  145       IF(RC .GT. DR) GO TO 140
             I     = I+1
             DR    = DR-RC
             XIN(I)= XN-(XN-XO)*DR
            IF(I .LT. NDM) GO TO 145
             DO 150 I  = 1,NDM
                XI(I,J)= XIN(I)
  150        CONTINUE
             XI(ND,J)  = ONE
  155     CONTINUE
********************************************* Updated Feb.08 '94
          IF( ITSX .GT. 0 ) THEN
              IPSAVE = 1
              XACC    = 1.0D70
              XTI     = 0.0D0
              XTSI    = XACC
              ITSX    = 1
              DO 200 J = 1, NDIM
              DO 200 I = 1, ND
                 XSAVE(I,J) = XI(I,J)
  200         CONTINUE
          ENDIF
C
      RETURN
      END
