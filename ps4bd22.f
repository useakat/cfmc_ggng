      subroutine ps4bd22(rnd,P1,mass,smin1,smax1,smin2,smax2,P,jacob
     &     ,wgt,ierr)
C     ****************************************************     
C     
C     The subroutine for
C     processes in color-flow-sampling option of MG. 
C     Input:
C     pp    4 momentum of external particles
C     Output:
C     Amplitude squared and summed
C
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     
C     CONSTANTS
C     
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     
C     ARGUMENTS 
C     
      real*8 rnd(8),P1(0:3),mass(4),smin1,smax1,smin2,smax2
      real*8 P(0:3,4),jacob,wgt
      integer ierr
C     
C     GLOBAL VARIABLES
C     
C     
C     LOCAL VARIABLES 
C     
      real*8 sqrts,costh,phi,K(0:3,2),KK(0:3,2),smass(2)
      real*8 jacob0,jacob1,jacob2,jacob3,wgt0,wgt1,wgt2,wgt3
C     
C     EXTERNAL FUNCTIONS
C     
      real*8 esbeta,es4dot,es4sq
      external esbeta,es4dot,es4sq
C     ----------
C     BEGIN CODE
C     ----------
      ierr = 0

      if (es4sq(P1).ge.0d0) then
         sqrts = dsqrt(es4sq(P1))
      else
         ierr = 1
         return
      endif

      smass(1) = dsqrt( smin1**2 +( smax1**2 -smin1**2 )*rnd(1) )
      smass(2) = dsqrt( smin2**2 +( smax2**2 -smin2**2 )*rnd(2) )
      wgt0 = ( smax1**2 -smin1**2 )/(2*pi)*( smax2**2 -smin2**2 )/(2*pi)
      jacob0 = 1d0

c      smass(1) = dsqrt( dexp( dlog(smin1**2) +( dlog(smax1**2) 
c     &     -dlog(smin1**2) )*rnd(1) ) )
c      smass(2) = dsqrt( dexp( dlog(smin2**2) +( dlog(smax2**2) 
c     &     -dlog(smin2**2) )*rnd(2) ) )
c      wgt0 = ( dlog(smax1**2) -dlog(smin1**2) )/(2*pi)
c     &     *( dlog(smax2**2) -dlog(smin2**2) )/(2*pi)
c      jacob0 = smass(1)**2*smass(2)**2

c      smass(1) = dsqrt( dexp( dlog(smin1**2) +( dlog(smax1**2) 
c     &     -dlog(smin1**2) )*rnd(1) ) )
c      smass(2) = dsqrt( dexp( dlog(smin2**2) +( dlog(smax2**2) 
c     &     -dlog(smin2**2) )*rnd(2) ) )
c      wgt = ( dlog(smax1**2) -dlog(smin1**2) )/(2*pi)
c     &     *( dlog(smax2**2) -dlog(smin2**2) )/(2*pi)
c      jacob = smass(1)**3*smass(2)**3

c      smass(1) = dsqrt( 1d0/( 1d0/smax1**2 +( 1d0/smin1**2 
c     &     -1d0/smax1**2 )*rnd(1) ) )
c      smass(2) = dsqrt( 1d0/( 1d0/(smax2**2) +( 1d0/(smin2**2) 
c     &     -1d0/(smax2**2) )*rnd(2) ) )
c      wgt0 = ( 1d0/(smin1**2) -1d0/(smax1**2) )/(2*pi)
c     &     *( 1d0/(smin2**2) -1d0/(smax2**2) )/(2*pi)
c      jacob0 = smass(1)**4*smass(2)**4

c      smass(1) = dsqrt( smin1**2 +( smax1**2 
c     &     -smin1**2 )*rnd(1) )
c      smass(2) = dsqrt( smin2**2 +( smax2**2
c     &     -smin2**2 )*rnd(2) )
c      wgt = ( smax1**2 -smin1**2 )/(2*pi)
c     &     *( smax2**2 -smin2**2 )/(2*pi)
c      jacob = 1d0

      call ps2bdt(rnd(3),P1,smass,K,jacob1,wgt1,ierr)
c      call ps2bd(rnd(3),P1,smass,K,jacob1,wgt1,ierr)
      call ps2bd(rnd(5),K(0,1),mass(1),P(0,1),jacob2,wgt2,ierr)
      call ps2bd(rnd(7),K(0,2),mass(3),P(0,3),jacob3,wgt3,ierr)

      wgt = wgt0*wgt1*wgt2*wgt3
c      jacob = jacob*jacob1*jacob2*jacob3
      jacob = jacob0*jacob1*jacob2*jacob3


      return
      end


      subroutine ps4bd22_2(rnd,P1,mass,smin1,smax1,smin2,smax2,P,jacob
     &     ,wgt,ierr)
C     ****************************************************     
C     
C     The subroutine for
C     processes in color-flow-sampling option of MG. 
C     Input:
C     pp    4 momentum of external particles
C     Output:
C     Amplitude squared and summed
C
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     
C     CONSTANTS
C     
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     
C     ARGUMENTS 
C     
      real*8 rnd(8),P1(0:3),mass(4),smin1,smax1,smin2,smax2
      real*8 P(0:3,4),jacob,wgt
      integer ierr
C     
C     GLOBAL VARIABLES
C     
C     
C     LOCAL VARIABLES 
C     
      real*8 sqrts,costh,phi,K(0:3,2),KK(0:3,2),smass(2)
      real*8 jacob1,jacob2,jacob3,wgt1,wgt2,wgt3
C     
C     EXTERNAL FUNCTIONS
C     
      real*8 esbeta,es4dot,es4sq
      external esbeta,es4dot,es4sq
C     ----------
C     BEGIN CODE
C     ----------
      ierr = 0

      if (es4sq(P1).ge.0d0) then
         sqrts = dsqrt(es4sq(P1))
      else
         ierr = 1
         return
      endif

      smass(1) = dsqrt( dexp( dlog(smin1**2) +( dlog(smax1**2) 
     &     -dlog(smin1**2) )*rnd(1) ) )
      smass(2) = dsqrt( dexp( dlog(smin2**2) +( dlog(smax2**2) 
     &     -dlog(smin2**2) )*rnd(2) ) )
      wgt = ( dlog(smax1**2) -dlog(smin1**2) )/(2*pi)
     &     *( dlog(smax2**2) -dlog(smin2**2) )/(2*pi)
      jacob = smass(1)**2*smass(2)**2

c      smass(1) = dsqrt( 1d0/( 1d0/smax1**2 +( 1d0/smin1**2 
c     &     -1d0/smax1**2 )*rnd(1) ) )
c      smass(2) = dsqrt( 1d0/( 1d0/(smax2**2) +( 1d0/(smin2**2) 
c     &     -1d0/(smax2**2) )*rnd(2) ) )
c      wgt = ( 1d0/(smin1**2) -1d0/(smax1**2) )/(2*pi)
c     &     *( 1d0/(smin2**2) -1d0/(smax2**2) )/(2*pi)
c      jacob = smass(1)**4*smass(2)**4

c      smass(1) = dsqrt( smin1**2 +( smax1**2 
c     &     -smin1**2 )*rnd(1) )
c      smass(2) = dsqrt( smin2**2 +( smax2**2
c     &     -smin2**2 )*rnd(2) )
c      wgt = ( smax1**2 -smin1**2 )/(2*pi)
c     &     *( smax2**2 -smin2**2 )/(2*pi)
c      jacob = 1d0

      call ps2bdt_2(rnd(3),P1,smass,K,jacob1,wgt1,ierr)
c      call ps2bd(rnd(3),P1,smass,K,jacob1,wgt1,ierr)
      call ps2bd(rnd(5),K(0,1),mass(1),P(0,1),jacob2,wgt2,ierr)
      call ps2bd(rnd(7),K(0,2),mass(3),P(0,3),jacob3,wgt3,ierr)

      wgt = wgt*wgt1*wgt2*wgt3
      jacob = jacob*jacob1*jacob2*jacob3


      return
      end
