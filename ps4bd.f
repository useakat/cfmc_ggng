      subroutine ps4bdt(rnd,P1,mass,smin1,smax1,smin2,smax2,sgncos,P
     &     ,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(8),P1(0:3),mass(4),smin1,smax1,smin2,smax2
      real*8 P(0:3,4),jacob,wgt
      integer ierr,sgncos
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,K(0:3,2),KK(0:3,2),smass(2)
      real*8 jacob0,jacob1,jacob2,jacob3,wgt0,wgt1,wgt2,wgt3
      real*8 s1,s1l,s1lmin,s1lmax
      real*8 s2,s2l,s2lmin,s2lmax
      integer ipw
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4dot,es4sq
      external esbeta,es4dot,es4sq
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
         s2lmin = smin2
         s2lmax = smax2
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = s1l
         s2 = s2l
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = 1d0
      elseif(ipw.eq.1) then
         s1lmin = dsqrt(smin1)
         s1lmax = dsqrt(smax1)
         s2lmin = dsqrt(smin2)
         s2lmax = dsqrt(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = s1l**2
         s2 = s2l**2
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = 2*dsqrt(s1)*2*dsqrt(s2)
      elseif(ipw.eq.2) then
         s1lmin = dlog(smin1)
         s1lmax = dlog(smax1)
         s2lmin = dlog(smin2)
         s2lmax = dlog(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = dexp(s1l)
         s2 = dexp(s2l)
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = s1*s2
      elseif(ipw.eq.3) then
         s1lmin = 1d0/dsqrt(smin1)
         s1lmax = 1d0/dsqrt(smax1)
         s2lmin = 1d0/dsqrt(smin2)
         s2lmax = 1d0/dsqrt(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = 1d0/s1l**2
         s2 = 1d0/s2l**2
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = (-2*s1**1.5d0)*(-2*s2**1.5d0)
      elseif(ipw.eq.4) then
         s1lmin = 1d0/smin1
         s1lmax = 1d0/smax1
         s2lmin = 1d0/smin2
         s2lmax = 1d0/smax2
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = 1d0/s1l
         s2 = 1d0/s2l
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = (-s1**2)*(-s2**2)
      endif
      
      smass(1) = dsqrt(s1)
      smass(2) = dsqrt(s2)

      call ps2bdt(rnd(3),P1,smass,sgncos,K,jacob1,wgt1,ierr)
      call ps2bd(rnd(5),K(0,1),mass(1),P(0,1),jacob2,wgt2,ierr)
      call ps2bd(rnd(7),K(0,2),mass(3),P(0,3),jacob3,wgt3,ierr)

      wgt = wgt0*wgt1*wgt2*wgt3
      jacob = jacob0*jacob1*jacob2*jacob3

      return
      end


      subroutine ps4bd22(rnd,P1,mass,smin1,smax1,smin2,smax2,sgncos,P
     &     ,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(8),P1(0:3),mass(4),smin1,smax1,smin2,smax2
      real*8 P(0:3,4),jacob,wgt
      integer ierr,sgncos
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,K(0:3,2),KK(0:3,2),smass(2)
      real*8 jacob0,jacob1,jacob2,jacob3,wgt0,wgt1,wgt2,wgt3
      real*8 s1,s1l,s1lmin,s1lmax
      real*8 s2,s2l,s2lmin,s2lmax
      integer ipw
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4dot,es4sq
      external esbeta,es4dot,es4sq
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
         s2lmin = smin2
         s2lmax = smax2
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = s1l
         s2 = s2l
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = 1d0
      elseif(ipw.eq.1) then
         s1lmin = dsqrt(smin1)
         s1lmax = dsqrt(smax1)
         s2lmin = dsqrt(smin2)
         s2lmax = dsqrt(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = s1l**2
         s2 = s2l**2
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = 2*dsqrt(s1)*2*dsqrt(s2)
      elseif(ipw.eq.2) then
         s1lmin = dlog(smin1)
         s1lmax = dlog(smax1)
         s2lmin = dlog(smin2)
         s2lmax = dlog(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = dexp(s1l)
         s2 = dexp(s2l)
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = s1*s2
      elseif(ipw.eq.3) then
         s1lmin = 1d0/dsqrt(smin1)
         s1lmax = 1d0/dsqrt(smax1)
         s2lmin = 1d0/dsqrt(smin2)
         s2lmax = 1d0/dsqrt(smax2)
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = 1d0/s1l**2
         s2 = 1d0/s2l**2
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = (-2*s1**1.5d0)*(-2*s2**1.5d0)
      elseif(ipw.eq.4) then
         s1lmin = 1d0/smin1
         s1lmax = 1d0/smax1
         s2lmin = 1d0/smin2
         s2lmax = 1d0/smax2
         s1l = s1lmin +(s1lmax -s1lmin)*rnd(1)
         s2l = s2lmin +(s2lmax -s2lmin)*rnd(2)
         s1 = 1d0/s1l
         s2 = 1d0/s2l
         wgt0 = (s1lmax -s1lmin)/(2*pi)*(s2lmax -s2lmin)/(2*pi)
         jacob0 = (-s1**2)*(-s2**2)
      endif
      
      smass(1) = dsqrt(s1)
      smass(2) = dsqrt(s2)

      call ps2bdt(rnd(3),P1,smass,sgncos,K,jacob1,wgt1,ierr)
      call ps2bd(rnd(5),K(0,1),mass(1),P(0,1),jacob2,wgt2,ierr)
      call ps2bd(rnd(7),K(0,2),mass(3),P(0,3),jacob3,wgt3,ierr)

      wgt = wgt0*wgt1*wgt2*wgt3
      jacob = jacob0*jacob1*jacob2*jacob3

      return
      end
