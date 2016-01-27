      subroutine randperm_finmom(P)
C     ****************************************************     
C     
C     The subroutine for permutating final-state momenta
C     randomely.
C
C     inout:
C        P: momenta to be permutated
C
C      output:
C         P: randomely permutated momenta
C
C     By Yoshitaro Takaesu @KEK Apr.30 2012
C     ****************************************************
      implicitnone
C     
C     GLOBAL VARIABLES
C     
      include 'cparam.inc'
C     
C     ARGUMENTS
C     
      real*8 P(0:3,next)
C     
C     LOCAL VARIABLES 
C     
      integer i
      integer nperm,ri,a(next),ipflag
C     
C     EXTERNAL FUNCTIONS
C     
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
      nperm = fact(nfin)
      ri = int(nperm*rand())+1
      do i = 1,next
         a(i) = i
      enddo
c      call iperm(ri,nfin,a(3))
      do i = 1,ri
         if(i.ne.1) then
            call ipnext(a(3),nfin,ipflag)
         endif
      enddo
      call porder(next,P,a)

      return
      end
