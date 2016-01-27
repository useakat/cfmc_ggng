CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This subroutine generate 2 body phase space in a CM frame, 
C beta/(8Pi)*dcos/2*dphi/(2Pi).
C inputs 
C   s: totall energy of the system
C   z(2): random numbers
C   mass(4): masses of external particles
C   PTcut: pT cut value
C   ETAcut: rapidity cut value
C
C outputs  
C   Q(0:3,4): momenta of external particles
C   wspn: weight of the phase space
C   ierr: error flag. ierr = 0: OK, ierr=1: momenta does not 
C   satisfy user defined cuts
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ps2b(s,z,mass,PTcut,ETAcut, Q,wpsn,ierr)

      implicit none

      integer i
      real*8 pi
      integer ierr,n
      real*8 s,z(2),wpsn,mass(4)
      DOUBLE PRECISION Q(0:3,50)
      real*8 PTcut,ETAcut
      real*8 a1,a2,beta,costh,phi,sqrts
      real*8 PT(4),eta(4)

      pi = dacos(-1d0)
      ierr = 0
      sqrts = dsqrt(s)

      a1 = mass(3)**2/s
      a2 = mass(4)**2/s
      beta = dsqrt(1d0 +a1**2 +a2**2 -2*(a1+a2+a1*a2))
      beta = 1d0
      wpsn = beta/(8*pi)

      costh = -1d0 +2*z(1)
      phi = 2*pi*z(2)

      call mom2cx(sqrts,0d0,0d0,1d0,0d0, Q(0,1),Q(0,2))
      call mom2cx(sqrts,0d0,0d0,costh,phi, Q(0,3),Q(0,4))

      PT(3) = dsqrt(Q(1,3)**2 +Q(2,3)**2)
      PT(4) = dsqrt(Q(1,4)**2 +Q(2,4)**2)
      eta(3) = 0.5d0*dlog((Q(0,3)+Q(3,3))/(Q(0,3)-Q(3,3)))
      eta(4) = 0.5d0*dlog((Q(0,4)+Q(3,4))/(Q(0,4)-Q(3,4)))

      do i = 3,4
         if (PT(i).lt.PTcut) ierr = 1
         if (eta(i).gt.ETAcut) ierr = 1
      enddo

      return
      end
