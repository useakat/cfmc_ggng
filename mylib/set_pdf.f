      subroutine set_pdf(lpp,pdg,x,Q, pdf)
C     ****************************************************
C     The subroutine for setting PDFs
C
C     Input:
C           lpp(2): beam types; 0: NO PDF, 1: proton, -1: anti-proton       
C           pdg(2): PDG code of initial partons
C           x(2): longitudinal momentum fractions of initial partons
C           Q(2): factorization scales
C
C     Output:
C            pdf(2): PDF of initial partons
C
C     By Yoshitaro Takaesu @KEK Dec.19 2011
C     ****************************************************
      implicitnone
C     
C     ARGUMENTS 
C     
      integer lpp(2),pdg(2)
      real*8 pdf(2),x(2),Q(2)
C     
C     LOCAL VARIABLES 
C     
      integer lp
      real*8 xpq(-7:7)
C     ----------
C     BEGIN CODE
C     ----------
      call pftopdg(abs(lpp(1)),x(1),Q(1),xpq)
      lp = sign(1,lpp(1))
      pdf(1) = xpq(pdg(1)*lp)
      call pftopdg(abs(lpp(2)),x(2),Q(2),xpq)
      lp = sign(1,lpp(2))
      pdf(2) = xpq(pdg(2)*lp)

      end
