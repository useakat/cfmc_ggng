      subroutine ps2bd(rnd,P1,mass,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2),ZP1(0:3)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,eps,maxcosth,minx,maxx,x
      real*8 jacob1,jacob2,wgt1,wgt2 
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
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
      costh = -1d0 +2d0*rnd(1)
      jacob1 = 1d0
      wgt1 = 2d0

      phi = 2*pi*rnd(2)
      jacob2 = 1d0
      wgt2 = 2*pi

      call mom2cx(sqrts,mass(1),mass(2),costh,phi, PPP(0,1),PPP(0,2))
      call boostx(PPP(0,1),P1,P(0,1))
      call boostx(PPP(0,2),P1,P(0,2))

      wgt = esbeta((mass(1)/sqrts)**2,(mass(2)/sqrts)**2)/(8*pi)/(4*pi)
     &     *wgt1*wgt2
      jacob = jacob1*jacob2

      return
      end


      subroutine ps2bd1(rnd,P1,mass,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2),ZP1(0:3)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,eps,maxcosth,minx,maxx,x
      real*8 jacob1,jacob2,wgt1,wgt2 
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
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
c      costh = -1d0 +2d0*rnd(1)
c      jacob1 = 1d0
c      wgt1 = 2d0
      eps = 1d-2
      maxcosth = 1d0 -eps
      minx = 0.5*dlog((1d0-maxcosth)/(1d0+maxcosth))
      maxx = 0.5*dlog((1d0+maxcosth)/(1d0-maxcosth))
      x = minx +(maxx -minx)*rnd(1)
      costh = (dexp(2*x) -1d0)/(dexp(2*x) +1d0) 
      jacob1 = (1d0 -costh**2)
      wgt1 = maxx -minx

      phi = 2*pi*rnd(2)
      jacob2 = 1d0
      wgt2 = 2*pi

      call mom2cx(sqrts,mass(1),mass(2),costh,phi, PPP(0,1),PPP(0,2))
      call boostx(PPP(0,1),P1,P(0,1))
      call boostx(PPP(0,2),P1,P(0,2))

      wgt = esbeta((mass(1)/sqrts)**2,(mass(2)/sqrts)**2)/(8*pi)/(4*pi)
     &     *wgt1*wgt2
      jacob = jacob1*jacob2

      return
      end

      subroutine ps2bd2(rnd,P1,mass,ptcut,rcut,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2),ZP1(0:3)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 m12,costh,pphi,eps,minx,maxx,x
      real*8 jacob1,jacob2,wgt1,wgt2,Et12,gamma,gambeta,sinth 
      real*8 beta,pt12,mincosth,maxcosth,ptcut,rcut,a1,a2
      real*8 sinphi_max,maxphi,sinphi_min,minphi
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq,pt,y,phi
      external esbeta,es4sq,pt,y,phi
C     ----------
C     BEGIN CODE
C     ----------
      ierr = 0

      if (es4sq(P1).ge.0d0) then
         m12 = dsqrt(es4sq(P1))
      else
         ierr = 1
         return
      endif
      pt12 = pt(P1)
c      y12 = y(P1)
c      phi12 = phi(P1)

      ET12 = dsqrt(m12**2 +pt12**2)
      gamma = ET12/m12
      gambeta = pt12/m12     

c     generation of costh*
c
       if (m12.lt.2*ptcut) then
          maxcosth = 1d0/pt12*(ET12 -2*ptcut)
          if (rcut.le.2*m12/pt12) then
             mincosth = 0d0
          else
             mincosth = pt12/ET12*dsqrt(1d0 -(2*m12/(pt12*rcut))**2)
          endif
          
          if (rnd(1).lt.0.5) then
             costh = -maxcosth +(-mincosth +maxcosth)*2*rnd(1)
          else
             costh = mincosth +(maxcosth -mincosth)*2*(rnd(1)-0.5d0)
          endif
          wgt1 = (maxcosth -mincosth)*2
          jacob1 = 1d0
       else
          maxcosth = 1d0
          mincosth = -1d0
          costh = mincosth +(maxcosth -mincosth)*rnd(1)
          wgt1 = (maxcosth -mincosth)
          jacob1 = 1d0
       endif
       if ((costh.ge.-1d0).and.(costh.le.1d0)) then
          continue
       else
c          write(*,*) "ierr=1: costh is out of range"
c          write(*,*) "costh =",costh
          ierr = 1
          return
       endif
       sinth = dsqrt(1d0 -costh**2)

c     generation of phi*
c
       if (m12.lt.2*ptcut) then
c      checking if the argument of dsqrt is positive
          a1 = 2*dsqrt(4*m12**4 +costh**2*pt12**2*ET12**2*rcut**4)
          a2 = rcut**2*(pt12**2 +costh**2*ET12**2)
          if (a1.lt.a2) then
c             write(*,*) "ierr=1: a1 < a2"
c             write(*,*) a1,a2
             ierr = 1
             return
          endif
          
          sinphi_max = 1d0/(sinth*m12*rcut)*dsqrt(a1-a2)
          if (sinphi_max.le.1d0) then
             maxphi = dasin(sinphi_max)
          else
             maxphi = pi/2d0
          endif
 
          if (costh.le.0d0) then
             if (ptcut.ge.0.5*dabs(pt12 +ET12*costh)) then
                sinphi_min = 1d0/sinth*dsqrt( (2*ptcut/m12)**2 
     &               -gamma**2*(gambeta/gamma +costh)**2 )
                minphi = dasin(sinphi_min)
             else
                minphi = 0d0
             endif
          else
             if (ptcut.ge.0.5d0*dabs(pt12 -ET12*costh)) then
                sinphi_min = 1d0/sinth*dsqrt( (2*ptcut/m12)**2 
     &               -gamma**2*(gambeta/gamma -costh)**2 )
                minphi = dasin(sinphi_min)
             else
                minphi = 0d0
             endif
          endif

          if (rnd(2).lt.0.25d0) then
             pphi = minphi +4*rnd(2)*(maxphi-minphi)
          elseif (rnd(2).lt.0.5d0) then
             pphi = pi -maxphi +4*(rnd(2)-0.25d0)*(maxphi-minphi)
          elseif (rnd(2).lt.0.75d0) then
             pphi = pi +minphi +4*(rnd(2)-0.5d0)*(maxphi-minphi)
          elseif (rnd(2).le.1d0) then
             pphi = 2*pi -maxphi +4*(rnd(2)-0.75d0)*(maxphi-minphi)
          endif
          wgt2 = 4*(maxphi-minphi)
          jacob2 = 1d0
       else
          maxphi = 2*pi
          minphi = 0d0
          pphi = minphi +(maxphi -minphi)*rnd(2)
          wgt2 = (maxphi-minphi)
          jacob2 = 1d0
       endif
       if ((pphi.ge.0d0).and.(pphi.le.2*pi)) then
          continue
       else
c          write(*,*) "IFLG=1: phi is out of range"
c          write(*,*) "phi =",pphi
          ierr = 1
          return
       endif

c     construct 4-momentum of the particle 1 and 2 in the center of mass frame of the particle 12
c
       beta = esbeta(mass(1)**2/m12**2,mass(2)**2/m12**2)
       PPP(0,1) = m12/2d0*(1d0 +(mass(1)**2 -mass(2)**2)/m12**2)
       PPP(1,1) = m12/2d0*beta*costh
       PPP(2,1) = m12/2d0*beta*sinth*dsin(pphi)
       PPP(3,1) = -m12/2d0*beta*sinth*dcos(pphi)
       PPP(0,2) = m12/2d0*(1d0 +(mass(2)**2 -mass(1)**2)/m12**2)
       PPP(1,2) = -PPP(1,1)
       PPP(2,2) = -PPP(2,1)
       PPP(3,2) = -PPP(3,1)

       call boost2(PPP(0,1),P1,P(0,1))
       call boost2(PPP(0,2),P1,P(0,2))
       
      wgt = beta/(8*pi)/(4*pi)*wgt1*wgt2
      jacob = jacob1*jacob2

      return
      end


      subroutine ps2bd2_costh(rnd,P1,mass,ptcut,rcut,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2),ZP1(0:3)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,eps,minx,maxx,x
      real*8 jacob1,jacob2,wgt1,wgt2,Et12,gamma,gambeta,sinth 
      real*8 beta,pt12,mincosth,maxcosth,ptcut,rcut
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
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

      pt12 = dsqrt(P1(1)**2 +P1(2)**2)

      ET12 = dsqrt(sqrts**2 +pt12**2)
      gamma = ET12/sqrts
      gambeta = pt12/sqrts    

c     generation of costh*
c
       if (sqrts.lt.2*ptcut) then
          maxcosth = 1d0/pt12*(ET12 -2*ptcut) 
c          maxcosth = 1d0
          if (rcut.le.2*sqrts/pt12) then
             mincosth = 0d0
          else
             mincosth = pt12/ET12*dsqrt(1d0 -(2*sqrts/(pt12*rcut))**2)
c             mincosth = 0d0
          endif
          
          if (rnd(1).lt.0.5) then
             costh = -maxcosth +(-mincosth +maxcosth)*2*(rnd(1)-0d0)
          else
             costh = mincosth +(maxcosth -mincosth)*2*(rnd(1)-0.5d0)
          endif
          wgt1 = (maxcosth -mincosth)*2
          jacob1 = 1d0
       else
          maxcosth = 1d0
          mincosth = -1d0
          costh = mincosth +(maxcosth -mincosth)*rnd(1)
          wgt1 = (maxcosth -mincosth)
          jacob1 = 1d0
       endif
       if ((costh.ge.-1d0).and.(costh.le.1d0)) then
          continue
       else
c          write(*,*) "IFLG=1: costh is out of range"
c          write(*,*) "costh =",costh
          ierr = 1
          return
       endif
       sinth = dsqrt(1d0 -costh**2)
c       write(*,*) costh,sinth


      phi = 2*pi*rnd(2)
      wgt2 = 2*pi
      jacob2 = 1d0
       if ((phi.ge.0d0).and.(phi.le.2*pi)) then
          continue
       else
c          write(*,*) "IFLG=1: phi is out of range"
c          write(*,*) "phi =",phi
          ierr = 1
          return
       endif
c       write(*,*) phi

c     construct 4-momentum of the particle 1 and 2 in the center of mass frame of the particle 12
c
       beta = esbeta(mass(1)**2/sqrts**2,mass(2)**2/sqrts**2)
       PPP(0,1) = sqrts/2d0*(1d0 +(mass(1)**2 -mass(2)**2)/sqrts**2)
       PPP(1,1) = sqrts/2d0*beta*costh
       PPP(2,1) = sqrts/2d0*beta*sinth*dsin(phi)
       PPP(3,1) = -sqrts/2d0*beta*sinth*dcos(phi)
       PPP(0,2) = sqrts/2d0*(1d0 +(mass(2)**2 -mass(1)**2)/sqrts**2)
       PPP(1,2) = -PPP(1,1)
       PPP(2,2) = -PPP(2,1)
       PPP(3,2) = -PPP(3,1)

       P1(0) = ET12
       P1(1) = pt12
       P1(2) = 0d0
       P1(3) = 0d0

       call boostx(PPP(0,1),P1,P(0,1))
       call boostx(PPP(0,2),P1,P(0,2))

      wgt = beta/(8*pi)/(4*pi)*wgt1*wgt2
      jacob = jacob1*jacob2

      return
      end


      subroutine ps2bd2_phi(rnd,P1,mass,ptcut,rcut,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2),ZP1(0:3)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,eps,minx,maxx,x
      real*8 jacob1,jacob2,wgt1,wgt2,Et12,gamma,gambeta,sinth 
      real*8 beta,pt12,mincosth,maxcosth,ptcut,rcut,a1,a2
      real*8 sinphi_max,maxphi,sinphi_min,minphi
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
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
      ptcut = 20d0
      rcut = 0.4d0

      pt12 = dsqrt(P1(1)**2 +P1(2)**2)

      ET12 = dsqrt(sqrts**2 +pt12**2)
      gamma = ET12/sqrts
      gambeta = pt12/sqrts
     
      costh = -1d0 +2*rnd(1)
      wgt1 = 2d0
      jacob1 = 1d0       
      if ((costh.ge.-1d0).and.(costh.le.1d0)) then
         continue
      else
c         write(*,*) "ierr=1: costh is out of range"
c         write(*,*) "costh =",costh
         ierr = 1
         return
      endif
      sinth = dsqrt(1d0 -costh**2)

c     generation of phi*
c
       if (sqrts.lt.2*ptcut) then
c      checking if the argument of dsqrt is positive
          a1 = 2*dsqrt(4*sqrts**4 +costh**2*pt12**2*ET12**2*rcut**4)
          a2 = rcut**2*(pt12**2 +costh**2*ET12**2)
          if (a1.lt.a2) then
c             write(*,*) "ierr=1: a1 < a2"
c             write(*,*) a1,a2
             ierr = 1
             return
          endif
          
          sinphi_max = 1d0/(sinth*sqrts*rcut)*dsqrt(a1-a2)
          if (sinphi_max.le.1d0) then
             maxphi = dasin(sinphi_max)
          else
             maxphi = pi/2d0
          endif
 
          if (costh.le.0d0) then
             if (ptcut.ge.0.5*dabs(pt12 +ET12*costh)) then
                sinphi_min = 1d0/sinth*dsqrt( (2*ptcut/sqrts)**2 
     &               -gamma**2*(gambeta/gamma +costh)**2 )
                minphi = dasin(sinphi_min)
             else
                minphi = 0d0
             endif
          else
             if (ptcut.ge.0.5d0*dabs(pt12 -ET12*costh)) then
                sinphi_min = 1d0/sinth*dsqrt( (2*ptcut/sqrts)**2 
     &               -gamma**2*(gambeta/gamma -costh)**2 )
                minphi = dasin(sinphi_min)
             else
                minphi = 0d0
             endif
          endif

c          minphi = minphi +(0d0-minphi)*1d0
c          maxphi = maxphi +(pi/2d0 -maxphi)*1d0
c          minphi = 0d0
c          maxphi = pi/2d0
          if (rnd(2).lt.0.25d0) then
             phi = minphi +4*rnd(2)*(maxphi-minphi)
          elseif (rnd(2).lt.0.5d0) then
             phi = pi -maxphi +4*(rnd(2)-0.25d0)*(maxphi-minphi)
          elseif (rnd(2).lt.0.75d0) then
             phi = pi +minphi +4*(rnd(2)-0.5d0)*(maxphi-minphi)
          elseif (rnd(2).le.1d0) then
             phi = 2*pi -maxphi +4*(rnd(2)-0.75d0)*(maxphi-minphi)
          endif
          wgt2 = 4*(maxphi-minphi)
          jacob2 = 1d0
       else
          maxphi = 2*pi
          minphi = 0d0
          phi = minphi +(maxphi -minphi)*rnd(2)
          wgt2 = (maxphi-minphi)
          jacob2 = 1d0
       endif
       if ((phi.ge.0d0).and.(phi.le.2*pi)) then
          continue
       else
c          write(*,*) "IFLG=1: phi is out of range"
c          write(*,*) "phi =",phi
          ierr = 1
          return
       endif

c     construct 4-momentum of the particle 1 and 2 in the center of mass frame of the particle 12
c
       beta = esbeta(mass(1)**2/sqrts**2,mass(2)**2/sqrts**2)
       PPP(0,1) = sqrts/2d0*(1d0 +(mass(1)**2 -mass(2)**2)/sqrts**2)
       PPP(1,1) = sqrts/2d0*beta*costh
       PPP(2,1) = sqrts/2d0*beta*sinth*dsin(phi)
       PPP(3,1) = -sqrts/2d0*beta*sinth*dcos(phi)
       PPP(0,2) = sqrts/2d0*(1d0 +(mass(2)**2 -mass(1)**2)/sqrts**2)
       PPP(1,2) = -PPP(1,1)
       PPP(2,2) = -PPP(2,1)
       PPP(3,2) = -PPP(3,1)

       call boostx(PPP(0,1),P1,P(0,1))
       call boostx(PPP(0,2),P1,P(0,2))

      wgt = beta/(8*pi)/(4*pi)*wgt1*wgt2
      jacob = jacob1*jacob2

      return
      end


      subroutine ps2bd2_0(rnd,P1,mass,ptcut,rcut,P,jacob,wgt,ierr)
C     ****************************************************     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     ****************************************************
      implicit none
C     CONSTANTS
      real*8 pi
      parameter(pi=3.14159265358979d0)
C     ARGUMENTS 
      real*8 rnd(2),P1(0:3),mass(2),P(0:3,2),jacob,wgt
      real*8 PP(0:3,2),PPP(0:3,2),ZP1(0:3)
      integer ierr
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      real*8 sqrts,costh,phi,eps,minx,maxx,x
      real*8 jacob1,jacob2,wgt1,wgt2,Et12,gamma,gambeta,sinth 
      real*8 beta,pt12,mincosth,maxcosth,ptcut,rcut,a1,a2
      real*8 sinphi_max,maxphi,sinphi_min,minphi
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
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

c      pt12 = dsqrt(P1(1)**2 +P1(2)**2)

c      ET12 = dsqrt(sqrts**2 +pt12**2)
c      gamma = ET12/sqrts
c      gambeta = pt12/sqrts
     
      costh = -1d0 +2*rnd(1)
      wgt1 = 2d0
      jacob1 = 1d0       
      if ((costh.ge.-1d0).and.(costh.le.1d0)) then
         continue
      else
c         write(*,*) "ierr=1: costh is out of range"
c         write(*,*) "costh =",costh
         ierr = 1
         return
      endif
      sinth = dsqrt(1d0 -costh**2)

      phi = 2*pi*rnd(2)
      wgt2 = 2*pi
      jacob2 = 1d0
      if ((phi.ge.0d0).and.(phi.le.2*pi)) then
         continue
      else
c         write(*,*) "IFLG=1: phi is out of range"
c         write(*,*) "phi =",phi
         ierr = 1
         return
      endif

c     construct 4-momentum of the particle 1 and 2 in the center of mass frame of the particle 12
c
       beta = esbeta(mass(1)**2/sqrts**2,mass(2)**2/sqrts**2)
       PPP(0,1) = sqrts/2d0*(1d0 +(mass(1)**2 -mass(2)**2)/sqrts**2)
       PPP(1,1) = sqrts/2d0*beta*costh
       PPP(2,1) = sqrts/2d0*beta*sinth*dsin(phi)
       PPP(3,1) = -sqrts/2d0*beta*sinth*dcos(phi)
       PPP(0,2) = sqrts/2d0*(1d0 +(mass(2)**2 -mass(1)**2)/sqrts**2)
       PPP(1,2) = -PPP(1,1)
       PPP(2,2) = -PPP(2,1)
       PPP(3,2) = -PPP(3,1)

       call boostx(PPP(0,1),P1,P(0,1))
       call boostx(PPP(0,2),P1,P(0,2))

      wgt = beta/(8*pi)/(4*pi)*wgt1*wgt2
      jacob = jacob1*jacob2

      return
      end


      subroutine ps2bd_in_pshk(rnd,P1,mass,P,wgt,jacob,iflg)
      implicitnone
      integer iflg
      real*8 rnd(2),P1(5),mass(2),P(5,2),wgt,jacob
      real*8 K(0:3),KK(0:3,2)

      K(0) = P1(4)
      K(1) = P1(1)
      K(2) = P1(2)
      K(3) = P1(3)
      call ps2bd(rnd,K,mass,KK,Jacob,wgt,iflg)
      P(1,1) = KK(1,1)
      P(2,1) = KK(2,1)
      P(3,1) = KK(3,1)
      P(4,1) = KK(0,1)
      P(5,1) = dsqrt( P(1,1)**2 +P(2,1)**2 +P(3,1)**2 ) 
      P(1,2) = KK(1,2)
      P(2,2) = KK(2,2)
      P(3,2) = KK(3,2)
      P(4,2) = KK(0,2)
      P(5,2) = dsqrt( P(1,2)**2 +P(2,2)**2 +P(3,2)**2 ) 

      return
      end


      subroutine ps2bd2_in_pshk(rnd,P1,mass,ptcut,rcut,P,wgt,jacob,iflg)
      implicitnone
      integer iflg
      real*8 rnd(2),P1(5),mass(2),P(5,2),wgt,jacob
      real*8 K(0:3),KK(0:3,2),ptcut,rcut

      K(0) = P1(4)
      K(1) = P1(1)
      K(2) = P1(2)
      K(3) = P1(3)
c      call ps2bd(rnd,K,mass,KK,Jacob,wgt,iflg)
      call ps2bd2(rnd,K,mass,ptcut,rcut,KK,Jacob,wgt,iflg)
c      call ps2bd2_costh(rnd,K,mass,ptcut,rcut,KK,Jacob,wgt,iflg)
c      call ps2bd2_phi(rnd,K,mass,ptcut,rcut,KK,Jacob,wgt,iflg)
c      call ps2bd2_0(rnd,K,mass,ptcut,rcut,KK,Jacob,wgt,iflg)
      P(1,1) = KK(1,1)
      P(2,1) = KK(2,1)
      P(3,1) = KK(3,1)
      P(4,1) = KK(0,1)
      P(5,1) = dsqrt( P(1,1)**2 +P(2,1)**2 +P(3,1)**2 ) 
      P(1,2) = KK(1,2)
      P(2,2) = KK(2,2)
      P(3,2) = KK(3,2)
      P(4,2) = KK(0,2)
      P(5,2) = dsqrt( P(1,2)**2 +P(2,2)**2 +P(3,2)**2 ) 

      return
      end
