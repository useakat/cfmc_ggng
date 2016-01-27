cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c KINEMATIC FUNCTIONS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function pt(p)
c************************************************************************
c     Returns transverse momentum of particle
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)

      pt = dsqrt(p(1)**2+p(2)**2)

      end

      double precision function et(p)
c************************************************************************
c     Returns transverse energy of particle
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)
      real*8 pt

      pt = dsqrt(p(1)**2+p(2)**2)
      if (pt .gt. 0d0) then
         et = p(0)*pt/dsqrt(pt**2+p(3)**2)
      else
         et = 0d0
      endif
      end

      double precision function e(p)
c************************************************************************
c     Returns energy of particle in lab frame
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)

      e = p(0)

      end

      double precision function mom(p)
c************************************************************************
c     Returns absolute value of momentum of particle in lab frame
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)

      mom = dsqrt(p(1)**2+p(2)**2+p(3)**2)

      end

      DOUBLE PRECISION  FUNCTION y(p)
c************************************************************************
c     Returns rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision pm

      pm = p(0)
      y = .5d0*dlog((pm+p(3))/(pm-p(3)))

      end

       DOUBLE PRECISION  FUNCTION eta(p)
c************************************************************************
c     Returns pseudo-rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision costh

      if (abs(costh(p)).lt.0.9999999999d0)then
         eta = -dlog(dtan(0.5d0*dacos(costh(p))))
      elseif (costh(p).ge.0.9999999999d0)then
         eta =  10d0
      else
         eta = -10d0
      endif

      end

       DOUBLE PRECISION  FUNCTION costh(p)
c************************************************************************
c     Returns the cosine of the angle between the particle and the +z axis
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision pt,pm

      pm = dsqrt(p(3)**2+pt(p)**2)
      if (abs(pm-p(3)).ge.abs(1d-10*p(3)))then
         costh = p(3)/pm
      else
         costh = 1d0
      endif

      end

      double precision function phi(p)
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
         phi =  dacos(temp)
      else
         phi = -dacos(temp) + 2d0*pi
      endif

      end

      double precision function dRij(P1,P2)
c************************************************************************
c     Distance in eta,phi between two particles.
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision ETA,DELTA_PHI,R2
      external eta,delta_phi

      R2 = (DELTA_PHI(P1,P2))**2+(ETA(p1)-ETA(p2))**2
      dRij = dsqrt(R2)

      RETURN
      END

      double precision function dRyij(P1,P2)
c************************************************************************
c     Distance in y,phi between two particles.
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision Y,DELTA_PHI,R2
      external y,delta_phi

      R2 = (DELTA_PHI(P1,P2))**2+(Y(p1)-Y(p2))**2
      dRyij = dsqrt(R2)

      RETURN
      END

      double precision function Ryij(P1,P2)
c************************************************************************
c     Distance in y,phi between two particles.
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision Y,DELTA_PHI,R2
      external y,delta_phi

      R2 = 2*( dcosh(Y(p1)-Y(p2)) -dcos(DELTA_PHI(P1,P2)) )
      Ryij = dsqrt(R2)

      RETURN
      END

      double precision function mij(P1,P2)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision sumdot
      external sumdot

      mij = dsqrt(Sumdot(p1,p2,1d0))

      RETURN
      END

      double precision function mijk(P1,P2,P3)
c************************************************************************
c     Invarient mass of 3 particles
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3),p3(0:3)
      double precision sumdot3
      external sumdot3

      mijk = dsqrt(Sumdot3(p1,p2,p3))

      RETURN
      END


      double precision function cosij(P1,P2)
c************************************************************************
c     Angle between direction of P1 in the restframe of P1+P2 and the 
c     direction of P1+P2 in the labframe.
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision ptot(0:3),pb1(0:3),l1,l2
      integer i

      do i=0,3
         ptot(i)=-p1(i)-p2(i)
      enddo
      ptot(0)=-ptot(0)

      call boostx(p1,ptot,pb1)

      l1 = dSqrt( pb1(1)**2+ pb1(2)**2+ pb1(3)**2)
      l2 = dSqrt(ptot(1)**2+ptot(2)**2+ptot(3)**2)

      if (l1.gt.0d0.and.l2.gt.0d0)then
         cosij =-(pb1(1)*ptot(1)+pb1(2)*ptot(2)+pb1(3)*ptot(3))/l1/l2
      else
         cosij = 0d0
      endif

      RETURN
      END

      double precision function delta_phi(p1,p2)
c************************************************************************
c     Returns separation in phi of two particles
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      REAL*8 DENOM, TEMP

      DENOM = SQRT(P1(1)**2 + P1(2)**2) * SQRT(P2(1)**2 + P2(2)**2)
      TEMP = MAX(-0.99999999D0, (P1(1)*P2(1) + P1(2)*P2(2)) / DENOM)
      TEMP = MIN( 0.99999999D0, TEMP)
      DELTA_PHI = ACOS(TEMP)

      END

      double precision function delta_eta(p1,p2)
c************************************************************************
c     Returns separation in eta of two particles
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision eta

      delta_eta=abs(eta(p1)-eta(p2))

      end

      double precision function kTij(p1,p2)
c************************************************************************
c     Returns kT distance of two particles (jet algorithm)
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision pt,dRij

      kTij = max( pt(p1) , pt(p2) )*dRij(p1,p2)

      end

      double precision function pTij(p1,p2)
c************************************************************************
c     Returns pT distance of two particles (jet algorithm)
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision pt,dRyij

      pTij = min( pt(p1) , pt(p2) )*dRyij(p1,p2)

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c AUXILIARY FUNCTIONS & SUBROUTINES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function dot(p1,p2)
C****************************************************************************
C     4-Vector Dot product
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      dot = p1(0)*p2(0) -p1(1)*p2(1) -p1(2)*p2(2) -p1(3)*p2(3)
      end

      double precision function SumDot(P1,P2,dsign)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3),dsign
      integer i
      double precision ptot(0:3)
      double precision dot
      external dot

      do i=0,3
         ptot(i)=p1(i)+dsign*p2(i)
      enddo
      SumDot = dot(ptot,ptot)

      RETURN
      END

      double precision function SumDot3(P1,P2,P3)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3),P3(0:3)
      integer i
      double precision ptot(0:3)
      double precision dot
      external dot

      do i=0,3
         ptot(i)=p1(i)+p2(i)+p3(i)
      enddo
      SumDot3 = dot(ptot,ptot)

      RETURN
      END
