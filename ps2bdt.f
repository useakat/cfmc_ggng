      subroutine ps2bdt(rnd,P1,mass,sgncos,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      integer sgncos
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer ipw
      real*8 sqrts,costh,phi,s1,s2,sshat,adenom,a,ll,lmin,lmax,logl
      real*8 bbeta,ZP1(0:3),wgt0,wgt1,wgt2,jacob0,jacob1,jacob2
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
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

      s1 = mass(1)**2
      s2 = mass(2)**2
      sshat = sqrts**2
      bbeta = esbeta(s1/sshat,s2/sshat)
      adenom = 1d0 -(s1+s2)/sshat +bbeta
      a = 4*s1*s2/sshat**2/adenom

      if (sgncos.eq.1) then 
         if (ipw.eq.0) then
            costh = -1d0 +2*rnd(1)
            jacob1 = 1d0
            wgt1 = 2d0
         elseif (ipw.eq.1) then
            lmin = dsqrt(a/bbeta +2d0)
            lmax = dsqrt(a/bbeta) 
            logl = lmin +(lmax -lmin)*rnd(1)
            ll = logl**2
            jacob1 = -2*dsqrt(ll)
            wgt1 = lmax -lmin
            costh = 1d0 +a/bbeta -ll
         elseif (ipw.eq.2) then
            lmin = dlog(a/bbeta +2d0)
            lmax = dlog(a/bbeta)
            logl = lmin +(lmax -lmin)*rnd(1)
            ll = dexp(logl)
c            jacob1 = -ll
c            wgt1 = lmax -lmin
            jacob1 = ll
            wgt1 = -(lmax -lmin)
            costh = 1d0 +a/bbeta -ll
         endif
      elseif (sgncos.eq.-1) then
         if (ipw.eq.0) then
            costh = -1d0 +2*rnd(1)
            jacob1 = 1d0
            wgt1 = 2d0
         elseif (ipw.eq.1) then
            lmin = dsqrt(a/bbeta)
            lmax = dsqrt(a/bbeta +2d0) 
            logl = lmin +(lmax -lmin)*rnd(1)
            ll = logl**2
            jacob1 = 2*dsqrt(ll)
            wgt1 = lmax -lmin
            costh = -1d0 -a/bbeta +ll
         elseif (ipw.eq.2) then
            lmin = dlog(a/bbeta)
            lmax = dlog(a/bbeta +2d0)
            logl = lmin +(lmax -lmin)*rnd(1)
            ll = dexp(logl)
            jacob1 = ll
            wgt1 = lmax -lmin
            costh = -1d0 +ll -a/bbeta
         endif
      else
         ierr = 1
         return
      endif

c      lmin = 1d0/dsqrt(a+2*bbeta)
c      lmax = 1d0/dsqrt(a)
c      logl = lmin +(lmax -lmin)*rnd(1)
c      ll = 1d0/logl**2
c      jacob = 2*ll**2d0/bbeta
c      wgt = lmax -lmin
c      costh = 1d0 -( ll -a )/bbeta      

c      lmin = 1d0/(a+2*bbeta) 
c      lmax = 1d0/a
c      logl = lmin +(lmax -lmin)*rnd(1)
c      ll = 1d0/logl
c      jacob = 1d0/bbeta*ll**2
c      wgt = lmax -lmin
c      costh = 1d0 -( ll -a )/bbeta

      phi = 2*pi*rnd(2)
      jacob2 = 1d0
      wgt2 = 2*pi

      call mom2cx(sqrts,mass(1),mass(2),costh,phi, PPP(0,1),PPP(0,2))
c      call mom2cx(sqrts,mass(1),mass(2),costh,phi, P(0,1),P(0,2))
c      ZP1(0) = P1(0)
c      ZP1(1) = 0d0
c      ZP1(2) = 0d0
c      ZP1(3) = dsqrt(P1(1)**2+P1(2)**2+P1(3)**2)
      call boostx(PPP(0,1),P1,P(0,1))
c      call rotxxx(PP(0,1),P1,P(0,1))
      call boostx(PPP(0,2),P1,P(0,2))
c      call rotxxx(PP(0,2),P1,P(0,2))

      jacob0 = 1d0
      wgt0 = esbeta((mass(1)/sqrts)**2,(mass(2)/sqrts)**2)/(8*pi)
     &     /2d0/(2*pi)

      jacob = jacob0*jacob1*jacob2
      wgt = wgt0*wgt1*wgt2


      return
      end

