      subroutine analyse_optch(cfused)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 25 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'chinfo.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer cfused(ngluons-1)
C     LOCAL VARIABLES 
      integer i,nch_bcf
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nch_bcf = 0
      do i = 1,ngluons-1
         if (cfused(i).eq.1) then
            nch_bcf = nch_bcf +1
            opt_cf(nch_bcf,nch_opt) = i
            cf_nopt(i) = cf_nopt(i) +1
            cf_opt(cf_nopt(i),i) = nch_opt
         endif
      enddo
      opt_ncf(nch_opt) = nch_bcf

      return
      end
