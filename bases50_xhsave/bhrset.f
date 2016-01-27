************************************************************************
*    ====================                                              *
      SUBROUTINE BHRSET
*    ====================                                              *
* ((Purpose))                                                          *
*     To reset contents of histograms and scatter plots.               *
* ((Author))                                                           *
*     S.Kawabata  June '90 at KEK                                      *
************************************************************************

      IMPLICIT NONE

      INCLUDE 'plot.h'

      INTEGER IP2

      INTEGER I,J
*                                                                      *
*--------------------------- Entry point ------------------------------*
*                                                                      *
*-------------------------- Clear Histograms --------------------------*
*                                                                      *
         DO 200 J    = 1, NHIST
           IP2       = MAPL(3,J)
           DO 100 I  = IP2,IP2+259
             IBUF(I) = 0
  100      CONTINUE
           IFBASE(J) = 0
  200    CONTINUE
*                                                                      *
*-------------------------- Clear Scat. Plots -------------------------*
*                                                                      *
         DO 500  J   = 1, NSCAT
           IP2       = MAPD(3,J)
           DO 400  I = IP2,IP2+2500
             IBUF(I) = 0.0
  400      CONTINUE
  500    CONTINUE
*                                                                      *
      RETURN
      END
