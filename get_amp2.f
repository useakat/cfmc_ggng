      subroutine get_amp2(Q,amp2)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'cparam.inc'
      include 'opt_info.inc'
C     CONSTANTS
C     ARGUMENTS 
      real rand
      real*8 Q(0:3,next),amp2
C     LOCAL VARIABLES 
      integer i,icflow,ichopt
      real*8 matele
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (icf.eq.0) then
c         call smatrix_mg_gg4g(Q,amp2)
         call smatrix_mg_per(Q,amp2)
      elseif (icf.eq.1) then
         amp2 = 0d0
         do i = 1,nsample_cf
cccc
cccc For amplitude obtained by Johan's recursive implementation (CF sampling)
cccc
c            icflow = int(ncf*rand()) +1 
c            call smatrix_co1_gg4g(icflow,Q,matele)
cccc
cccc  For amplitude obtined by CFQCD model (CF sum/sample)
cccc
c            call GetCF_i(i)          ! For Leading Color Sum
            call GetCF(dble(rand()))  ! For Leading Color Sample
            call smatrix_cfqcd(Q,matele)
            amp2 = amp2 +matele
         enddo
         amp2 = amp2/dble(nsample_cf)
      elseif (icf.eq.2) then
cccc
cccc     For the Optimized LCF sampling
cccc
c         call GetCF_opt
c         call smatrix_cfqcd(Q,amp2)
c
cccc
cccc     For the exact Optimized LCF sum
cccc
         ichopt = per_opt(ich)  !  ch number > opt-ch number
         cffact = opt_ncf(ichopt) ! the number of CFs for the opt-ch
         amp2 = 0d0
         do i = 1,cffact
            ichcf = opt_cf(i,ichopt) ! get the i_th CF for the opt-ch
            call GetCF_basic(ichcf)  ! get the CF ordering
            call smatrix_cfqcd(Q,matele) 
            call GetChWgt_opt(Q,ichcf)
            amp2 = amp2 +matele/AAA
c            amp2 = amp2 +matele
         enddo
c         call GetCF_opt_co1(dble(rand()),icflow)
c         call smatrix_co1(icflow,Q,amp2)
      elseif (icf.eq.3) then
         amp2 = 0d0
         do i = 1,next-1
            call GetCF_basic(i)
            call smatrix_cfqcd(Q,matele)
            amp2 = amp2 +matele
         enddo
      endif

      return
      end
