      program test
C     ****************************************************
C     
C     By Yoshitaro Takaesu @KEK Nov.19 2011
C     
C     The driver for making include and dat files for QCD
C     processes in color-flow-sampling option of MG. 
C
C     ****************************************************
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
C     CONSTANTS
C     ARGUMENTS 
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i,j,iflg
      real*8 P1(0:3,10),P2(0:3,10),zz(30),pini(0:3),pfin(0:3)
      real rand
      real*8 wgt,P12(5),P13(0:3),m12,pt12,ET12,pt1,pt2,eta12,phi12
      integer nnext,nnfin
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq,pt,dRyij,Ryij
      external esbeta,es4sq,pt,dRyij,Ryij
C     ----------
C     BEGIN CODE
C     ----------
      nnext = 8
      nnfin = nnext-2
      ss = (14000d0)**2
      ndim = 3*(nnext-2) -2
      do i = 1,53
         zz(1) = rand()
      enddo

      nptcut_fails = 0
      nrcut_fails = 0

      do j = 1,10

      do i = 1,ndim
         zz(i) = rand()
      enddo

      ptcut = 20d0
      ayecut = 5d0
      rcut = 0.4d0
      ptijcut = 0d0
      lptcut = 0d0
      ptcut2 = ptcut
      ayecut2 = ayecut
      ayecut3 = ayecut
      ayecut4 = ayecut
c      Qcut = dsqrt(2*ptcut**2*(1d0-dcos(rcut)))
      Qcut = ptcut*rcut
      Qmax = dsqrt(ss)
      pi = 3.1415926535d0
      do i = 1,nnfin
         ptmin_cut(i) = ptcut
         ptmax_cut(i) = Qmax/2d0
      enddo
c      mass(1) = 0d0
c      mass(2) = 0d0
c      do i = 1,next
      do i = 1,nnext
         M(1,i) = 0d0
         M(2,i) = M(1,i)**2
         M(3,i) = 0
      enddo
c$$$      m12 = 15d0
c$$$      pt12 = 75d0
c$$$      ET12 = dsqrt(m12**2 +pt12**2)
c$$$      eta12 = 1.2d0
c$$$      phi12 = 0.54d0;
c$$$
c$$$      P12(4) = ET12*dcosh(eta12)
c$$$      P12(1) = pt12*dcos(phi12)
c$$$      P12(2) = pt12*dsin(phi12)
c$$$      P12(3) = ET12*dsinh(eta12)
c$$$      P12(5) = dsqrt(P12(1)**2 +P12(2)**2 +P12(3)**2)
      call pshkpos2_3(ss,nnfin,zz,x1,x2,wpsn,jacob,iflg)
c      call ps2bd2_in_pshk(zz,P12,mass,ptcut,rcut,P,wgt,jacob,iflg)
      do i = 1,nnext
         P1(0,i) = P(4,i)
         P1(1,i) = P(1,i)
         P1(2,i) = P(2,i)
         P1(3,i) = P(3,i)
      enddo
      if (iflg.eq.0) then
         call mom_check(nnext,P1,M)
      endif
c$$$
c$$$      call ptcut_chk(2,P1,ptcut,iflg)
c$$$      if (iflg.ne.0) then
c$$$         nptcut_fails = nptcut_fails +1
c$$$         write(*,*) pt(P1(0,1)),pt(P1(0,2))
c$$$      endif
c$$$      call rcut_chk(2,P1,rcut,iflg)
c$$$c      call drcut_chk(2,P1,rcut,iflg)
c$$$      if (iflg.ne.0) then
c$$$         nrcut_fails = nrcut_fails +1
c$$$         write(*,*) Ryij(P1(0,1),P1(0,2))
c$$$      endif
c$$$
c$$$      Enddo
c$$$
c$$$      write(*,*) "pt_fails",nptcut_fails
c$$$      write(*,*) "dr_fails",nrcut_fails
c$$$
c$$$      P2(0,1) = m12
c$$$      P2(1,1) = 0d0
c$$$      P2(2,1) = 0d0
c$$$      P2(3,1) = 0d0
c$$$      P13(0) = ET12*dcosh(eta12)
c$$$      P13(1) = pt12*dcos(phi12)
c$$$      P13(2) = pt12*dsin(phi12)
c$$$      P13(3) = ET12*dsinh(eta12)
c$$$c      call boostx(P2(0,1),P13,P1(0,1))
c$$$      call boost2(P2(0,1),P13,P1(0,1))
c$$$      write(*,*) P13
c$$$      write(*,*) P1(0,1),P1(1,1),P1(2,1),P1(3,1)
      enddo

      end
