************************************************************************
*    ===========================                                       *
      SUBROUTINE SHFILL( NTRY )
*    ===========================                                       *
* ((Function))                                                         *
*     To fill the number of trials for a event generation              *
* ((Input))                                                            *
*    NTYR : the number of trials for the current event                 *
* ((Author))                                                           *
*    S.Kawabata    April 1994                                          *
*                                                                      *
************************************************************************

      IMPLICIT NONE

      INTEGER NTRY

      INCLUDE 'plotsp.h'


      IF( NTRY .LE. NBIN ) THEN
          IBUFSP( NTRY ) = IBUFSP( NTRY ) + 1
      ELSE
          IBUFSP( NBIN+1 ) = IBUFSP( NBIN+1 ) + 1
      ENDIF

      RETURN
      END
