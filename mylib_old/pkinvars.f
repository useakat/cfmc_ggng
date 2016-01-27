      real*8 function pmass(p)
      implicit none
      real*8 p(0:3)

      pmass = dsqrt( p(0)**2 -p(1)**2 -p(2)**2 -p(3)**2 )

      return
      end


      real*8 function ppt(p)
      implicit none
      real*8 p(0:3)

      ppt = dsqrt( p(1)**2 +p(2)**2 )

      return
      end


      real*8 function py(p)
      implicit none
      real*8 p(0:3)

      py = 0.5d0*dlog((p(0) +p(3))/(p(0) -p(3)))

      return
      end


       DOUBLE PRECISION  FUNCTION peta(p)
c************************************************************************
c     Returns pseudo-rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision pcosth

      if (abs(pcosth(p)).lt.0.9999999999d0)then
         peta = -dlog(dtan(0.5d0*dacos(pcosth(p))))
      elseif (pcosth(p).ge.0.9999999999d0)then
         peta =  10d0
      else
         peta = -10d0
      endif

      end


       DOUBLE PRECISION  FUNCTION pcosth(p)
c************************************************************************
c     Returns the cosine of the angle between the particle and the +z axis
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision ppt,pm

      pm = dsqrt(p(3)**2+ppt(p)**2)
      if (abs(pm-p(3)).ge.abs(1d-10*p(3)))then
         pcosth = p(3)/pm
      else
         pcosth = 1d0
      endif

      end


      double precision function pdelphi(phi1, phi2)
      implicit none

      double precision phi1, phi2

      double precision x1,y1
      double precision x2,y2

      double precision sinphi,cosphi
      double precision phis,phic
      
      x1 = dcos(phi1)
      y1 = dsin(phi1)

      x2 = dcos(phi2)
      y2 = dsin(phi2)

      sinphi = x1*y2-y1*x2
      cosphi = x1*x2+y1*y2

      phis = dasin(sinphi)
      phic = dacos(cosphi)

      if (cosphi .ge. 0.d0) then
         pdelphi = phis
      else if (sinphi .ge. 0.d0) then
         pdelphi = phic
      else 
         pdelphi =  -phic
      endif

      end


      double precision function pdeltar(p1,p2)
      implicit none

      double precision p1(0:3),p2(0:3)

      double precision eta1,eta2
      double precision phi1,phi2

      double precision peta,pdelphi
      external peta,pdelphi

      eta1 = peta(p1)
      eta2 = peta(p2)

      phi1 = datan2(p1(2),p1(1))
      phi2 = datan2(p2(2),p2(1))

      pdeltar = dsqrt((eta1-eta2)**2 +pdelphi(phi1,phi2)**2)

      end


      double precision function pR(p1,p2)
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision eta1,eta2
      double precision phi1,phi2
      double precision peta,pdelphi
      external peta,pdelphi

      eta1 = peta(p1)
      eta2 = peta(p2)

      phi1 = datan2(p1(2),p1(1))
      phi2 = datan2(p2(2),p2(1))

      pR = dsqrt( 2*dcosh(eta1-eta2) -dcos(pdelphi(phi1,phi2)) )

      return
      end


      real*8 function pmij(p1,p2)
c by Yoshitaro Takaesu -Aug/03/2011 @KEK
      implicit none

      integer i

      real*8 p1(0:3),p2(0:3)

      real*8 q(0:3)

      real*8 fvdot
      external fvdot

      do i = 0,3
         q(i) = p1(i) +p2(i)
      enddo
      pmij = dsqrt(fvdot(q))

      end

      double precision function pktij(p1,p2)
c************************************************************************
c     Returns kT distance of two particles (jet algorithm)
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision ppt,pdeltar

      pktij = max( ppt(p1) , ppt(p2) )*pdeltar(p1,p2)

      end



      double precision function pphi(p)
c************************************************************************
c     Returns azimuthal angle phi of a particle
c************************************************************************
      implicit none
      double precision p(0:3)
      double precision denom,temp
      double precision pi
      parameter ( pi=3.1415926535897932385 )

      denom = dsqrt(p(1)**2 + p(2)**2)
      temp = max(-0.99999999d0, p(1) / denom)
      temp = min( 0.99999999d0, temp)
      if (p(2).ge.0d0)then
         pphi =  dacos(temp)
      else
         pphi = -dacos(temp) + 2d0*pi
      endif

      end


      subroutine boostx(p,q , pboost)
c
c This subroutine performs the Lorentz boost of a four-momentum.  The
c momentum p is assumed to be given in the rest frame of q.  pboost is
c the momentum p boosted to the frame in which q is given.  q must be a
c timelike momentum.
c
c input:
c       real    p(0:3)         : four-momentum p in the q rest  frame
c       real    q(0:3)         : four-momentum q in the boosted frame
c
c output:
c       real    pboost(0:3)    : four-momentum p in the boosted frame
c
      implicit none
      double precision p(0:3),q(0:3),pboost(0:3),pq,qq,m,lf

      double precision rZero
      parameter( rZero = 0.0d0 )

      qq = q(1)**2+q(2)**2+q(3)**2

      if ( qq.ne.rZero ) then
         pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
         m = sqrt(q(0)**2-qq)
         lf = ((q(0)-m)*pq/qq+p(0))/m
         pboost(0) = (p(0)*q(0)+pq)/m
         pboost(1) =  p(1)+q(1)*lf
         pboost(2) =  p(2)+q(2)*lf
         pboost(3) =  p(3)+q(3)*lf
      else
         pboost(0) = p(0)
         pboost(1) = p(1)
         pboost(2) = p(2)
         pboost(3) = p(3)
      endif
c
      return
      end


      subroutine boost2(p,q,pboost)
c
c This subroutine performs the spacetime coordinate transformation of a four-momentum.
c This routine performs 
c     1. boost along x-axis by gam = ET(q)/m(q)
c     2. boost along z-axis by the rapidity of q
c     3. x-y rotation by azmuthal angle of q
C The momentum p is assumed to be given in the rest frame of q.
c pboost is the momentum p boosted to the frame in which q is given.
c q must be atimelike momentum.
c
c input:
c       real    p(0:3)         : four-momentum p in the q rest  frame
c       real    q(0:3)         : four-momentum q in the boosted frame
c
c output:
c       real    pboost(0:3)    : four-momentum p in the boosted frame
c
      implicit none
      double precision p(0:3),q(0:3),pboost(0:3)
      real*8 beta,gam,ET,pt,y,phi,m
      double precision rZero
      parameter( rZero = 0.0d0 )
      real*8 ppt,py,pphi,pmass
      external ppt,py,pphi,pmass

      m = pmass(q) 
      pt = ppt(q)
      y = py(q)
      phi = pphi(q)

      ET = dsqrt(m**2 +pt**2)
      gam = ET/m
      beta = dsqrt(gam**2 -1d0)/gam

      pboost(0) = gam*p(0)*dcosh(y) +gam*beta*p(1)*dcosh(y) 
     &           +p(3)*dsinh(y)
      pboost(1) = gam*beta*p(0)*dcos(phi) +gam*p(1)*dcos(phi)
     &           -p(2)*dsin(phi)
      pboost(2) = p(2)*dcos(phi) +gam*beta*p(0)*dsin(phi) 
     &           +gam*p(1)*dsin(phi)
      pboost(3) = p(3)*dcosh(y) +gam*p(0)*dsinh(y) 
     &           +gam*beta*p(1)*dsinh(y)
      
      return
      end
