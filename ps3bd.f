      subroutine ps3bd2(rnd,P1,mass,smin1,smax1,sgncos,P,jacob,wgt,ierr)
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
C     ***************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(5),P1(0:3),mass(3),smin1,smax1,smin2,smax2
      real*8 P(0:3,3),jacob,wgt
      integer ierr,sgncos
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,K(0:3,2),KK(0:3,2),smass(2)
      real*8 jacob0,jacob1,jacob2,jacob3,wgt0,wgt1,wgt2,wgt3
      real*8 s1,s1l,s1lmin,s1lmax
      integer ipw
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4dot,es4sq
      external esbeta,es4dot,es4sq
      real*8 ampfact
      common /ampfact/ ampfact
C     ----------
C     BEGIN CODE
C     ----------
      ierr = 0
      ipw = 2

      if (es4sq(P1).ge.0d0) then
         sqrts = dsqrt(es4sq(P1))
      else
         ierr = 1
         return
      endif

      if (ipw.eq.0) then
         s1lmin = smin1
         s1lmax = smax1
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s1 = s1l
         wgt0 = (s1lmax -s1lmin)/(2*pi)
         jacob0 = 1d0
      elseif(ipw.eq.1) then
         s1lmin = dsqrt(smin1)
         s1lmax = dsqrt(smax1)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s1 = s1l**2
         wgt0 = (s1lmax -s1lmin)/(2*pi)
         jacob0 = 2*dsqrt(s1)
      elseif(ipw.eq.2) then
         s1lmin = dlog(smin1)
         s1lmax = dlog(smax1)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s1 = dexp(s1l)
         wgt0 = (s1lmax -s1lmin)/(2*pi)
         jacob0 = s1
      elseif(ipw.eq.3) then
         s1lmin = 1d0/dsqrt(smin1)
         s1lmax = 1d0/dsqrt(smax1)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s1 = 1d0/s1l**2
         wgt0 = (s1lmax -s1lmin)/(2*pi)
         jacob0 = -s1**1.5d0
      elseif(ipw.eq.4) then
         s1lmin = 1d0/smin1
         s1lmax = 1d0/smax1
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s1 = 1d0/s1l
         wgt0 = (s1lmax -s1lmin)/(2*pi)
         jacob0 = -s1**2
      endif
      
      smass(1) = dsqrt(s1)
      smass(2) = mass(3)

      call ps2bdt(rnd(2),P1,smass,sgncos,K,jacob1,wgt1,ierr)
      call ps2bd(rnd(4),K,mass,P,jacob2,wgt2,ierr)
      P(0,3) = K(0,2) 
      P(1,3) = K(1,2) 
      P(2,3) = K(2,2) 
      P(3,3) = K(3,2) 

      wgt = wgt0*wgt1*wgt2
      jacob = jacob0*jacob1*jacob2
c      ampfact = -s1**2 +smax1**2
c      ampfact = -(s1-smax1)**2 +smax1**2 
c      ampfact = dsqrt(14000**2-s1)

      return
      end
