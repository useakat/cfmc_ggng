      subroutine set_histo(xl,xu)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'cparam.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i
      character*80 ctit
      integer maxdim
      parameter (maxdim=50)
      real*8 xl(maxdim),xu(maxdim)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      call xhinit(20,0d0,1d0,50,"dsigma/dx1")
      call xhinit(21,0d0,1d0,50,"dsigma/dx2")
      call xhinit(22,0d0,2d0,50,"dsigma/M2/SumM2")
c      call dhinit(22,0d0,1d0,30,0d0,1d0,30,
c     &     'd sigma/ d x1/d x2')
      do i=1,ndim
         write(ctit,1010) i,xl(i),xu(i)
 1010    format('x(',i2,') ; ',e12.6,' - ',e12.6)
         call xhinit(i+30,xl(i),xu(i),50,ctit)
      enddo
      call xhinit(61,0d0,100d0,50,"dsigma/pt_g1")
      call xhinit(62,0d0,100d0,50,"dsigma/pt_g2")
      call xhinit(63,0d0,100d0,50,"dsigma/pt_g3")
      call xhinit(64,0d0,100d0,50,"dsigma/pt_g4")
      call xhinit(65,0d0,10d0,50,"dsigma/y_g1")
      call xhinit(66,0d0,10d0,50,"dsigma/y_g2")
      call xhinit(67,0d0,10d0,50,"dsigma/y_g3")
      call xhinit(68,0d0,10d0,50,"dsigma/y_g4")
      call xhinit(69,0d0,6.5d0,50,"dsigma/phi_g1")
      call xhinit(70,0d0,6.5d0,50,"dsigma/phi_g2")
      call xhinit(71,0d0,6.5d0,50,"dsigma/phi_g3")
      call xhinit(72,0d0,6.5d0,50,"dsigma/phi_g4")
      call xhinit(73,0d0,5d0,50,"dsigma/dr_g1g2")
      call xhinit(74,0d0,5d0,50,"dsigma/dr_g1g3")
      call xhinit(75,0d0,5d0,50,"dsigma/dr_g1g4")
      call xhinit(76,0d0,5d0,50,"dsigma/dr_g2g3")
      call xhinit(77,0d0,5d0,50,"dsigma/dr_g2g4")
      call xhinit(78,0d0,5d0,50,"dsigma/dr_g3g4")
      call xhinit(79,0d0,100d0,50,"dsigma/m_g1g2")
      call xhinit(80,0d0,100d0,50,"dsigma/m_g1g3")
      call xhinit(81,0d0,100d0,50,"dsigma/m_g1g4")
      call xhinit(82,0d0,100d0,50,"dsigma/m_g2g3")
      call xhinit(83,0d0,100d0,50,"dsigma/m_g2g4")
      call xhinit(84,0d0,100d0,50,"dsigma/m_g3g4")
      call xhinit(85,0d0,100d0,50,"dsigma/kt_g1g2")
      call xhinit(86,0d0,100d0,50,"dsigma/kt_g1g3")
      call xhinit(87,0d0,100d0,50,"dsigma/kt_g1g4")
      call xhinit(88,0d0,100d0,50,"dsigma/kt_g2g3")
      call xhinit(89,0d0,100d0,50,"dsigma/kt_g2g4")
      call xhinit(90,0d0,100d0,50,"dsigma/kt_g3g4")

      return
      end
