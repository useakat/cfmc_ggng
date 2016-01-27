************************************************************************
*                                                                      *
*    ==========================                                        *
      SUBROUTINE BSWRIT( LUN,LUN2 )
*    =====================                                             *
* ((Purpose))                                                          *
*     Read temporary result from disk file.                            *
* ((Auther))                                                           *
*     S.Kawabata  June '90 at KEK                                      *
*                                                                      *
************************************************************************

      IMPLICIT NONE
C     IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER LUN,LUN2

      INCLUDE 'base1.h'
      INTEGER ND1(5*MXDIM+3)
      EQUIVALENCE (ND1(1),XL(1))
      INCLUDE 'base3.h'
      INTEGER ND3(11)
      EQUIVALENCE (ND3(1),SCALLS)
      INCLUDE 'base4.h'
      INTEGER ND4(2*MXDIM*(NDMX+1)+4*LENG+MXDIM+3)
      EQUIVALENCE (ND4(1),XI(1,1))
      INCLUDE 'base5.h'
      INTEGER ND5(22*ITM)
      EQUIVALENCE (ND5(1),ITRAT(1,0))

      INCLUDE 'randm.h'

      INCLUDE 'plot.h'
      INTEGER NPH(18*(NHS+NSC)+29)
      EQUIVALENCE (NPH(1),XHASH(1,1))

      INCLUDE 'ninfo.h'

      INTEGER I


      IF( NODEID .NE. 0 ) RETURN

      REWIND LUN
      WRITE(LUN) ND1,ND3,ND4,ND5,ND6,NPH
C     WRITE(LUN) ND1,ND2,ND3,ND4,ND5,ND6,NPH
      IF(NW .EQ. 0 ) NW = 281
      WRITE(LUN) NW,(IBUF(I),I=1,NW)
      WRITE(LUN2,*) XI

      do i = 1,10
         write(99,*) "HeperCube",I
         write(99,*) "DXP",DXP(I)
         write(99,*) "DXD",DXD(I)
      enddo
      write(99,*) "NG",NG
      write(99,*) "ND",ND
C
      write(28,*) DXP
      RETURN
      END
