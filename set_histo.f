      subroutine set_histo(xl,xu)
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
      call xhinit(31,0d0,1d0,50,'dsigma/dx1')
      call xhinit(32,0d0,1d0,50,'dsigma/dx2')
      call xhinit(33,0d0,2d0,50,'dsigma/M2/SumM2')
c      call dhinit(22,0d0,1d0,30,0d0,1d0,30,
c     &     'd sigma/ d x1/d x2')
      do i=1,ndim
         write(ctit,1010) i,xl(i),xu(i)
 1010    format('x(',i2,') ; ',e12.6,' - ',e12.6)
         call xhinit(i+33,xl(i),xu(i),50,ctit)
      enddo
      call xhinit(41,0d0,100d0,50,'dsigma/pt_g1')
      call xhinit(42,0d0,100d0,50,'dsigma/pt_g2')
      call xhinit(43,0d0,100d0,50,'dsigma/pt_g3')
      call xhinit(44,0d0,10d0,50,'dsigma/y_g1')
      call xhinit(45,0d0,10d0,50,'dsigma/y_g2')
      call xhinit(46,0d0,10d0,50,'dsigma/y_g3')
      call xhinit(47,0d0,5d0,50,'dsigma/dr_g1g2')
      call xhinit(48,0d0,5d0,50,'dsigma/dr_g1g3')
      call xhinit(49,0d0,5d0,50,'dsigma/dr_g2g3')
      call xhinit(50,0d0,100d0,50,'dsigma/m_g1g2')
      call xhinit(51,0d0,100d0,50,'dsigma/m_g1g3')
      call xhinit(52,0d0,100d0,50,'dsigma/m_g2g3')

      return
      end
