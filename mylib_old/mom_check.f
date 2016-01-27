      subroutine mom_check(n,P,mass)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG.30 2013
C     
C     The subroutine to check the validity of momenta of 
C     scattering particles
C     ****************************************************
      implicitnone
C     CONSTANTS
C     ARGUMENTS 
      integer n
      real*8 P(0:3,n),mass(n)
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i
      real*8 pini(0:3),pfin(0:3)
C     EXTERNAL FUNCTIONS
      real*8 esbeta,es4sq
      external esbeta,es4sq
C     ----------
C     BEGIN CODE
C     ----------
      call write_mom(n,P,6)
      do i = 0,3
         pini(i) = P(i,1) +P(i,2)
         pfin(i) = P(i,3) +P(i,4) +P(i,5)
      enddo
      write(*,*) "Sum of initial momenta"
      call write_mom(1,pini,6)
      write(*,*) "Sum of final momenta"
      call write_mom(1,pfin,6)
      write(*,*) "4-momentum squared and mass^2"
      do i =1,n
         write(*,*) i,es4sq(P(0,i)),mass(i)**2
      enddo
      write(*,*)

      return
      end
