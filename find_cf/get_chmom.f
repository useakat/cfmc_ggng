      subroutine get_chmom(natm,atm,nsch)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 25 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'chinfo.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer natm,nsch
      integer atm(2,natm)
C     LOCAL VARIABLES 
      integer i,ipos1,ipos2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      ipos1 = 2
      ipos2 = ngluons -nsch*2
      do i = 1,natm
         if (i.le.2) then
            if (atm(2,i).ne.0) then
               ipos1 = ipos1 +1
               chmom(atm(2,i)-2,nch_per) = ipos1
            endif
         else
            if (atm(1,i).ne.0) then
               if (atm(2,i).eq.0) then
                  ipos1 = ipos1 +1
                  chmom(atm(1,i)-2,nch_per) = ipos1
               else
                  ipos2 = ipos2 +1
                  chmom(atm(2,i)-2,nch_per) = ipos2
                  ipos2 = ipos2 +1
                  chmom(atm(1,i)-2,nch_per) = ipos2
               endif
            endif
         endif
      enddo       

      return
      end
