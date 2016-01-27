      subroutine perm_atm(natm,atm)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer natm
      integer atm(2,natm)
C     LOCAL VARIABLES 
      integer nperm
C     EXTERNAL FUNCTIONS
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
      do i = 3,natm
         aatm(1,3) = atm(1,3)
         aatm(2,3) = atm(2,3)
      enddo
      if ((atm(2,1)*atm(2,2)).ne.0) then
         do i = 3,natm
            a(i) = i
         enddo
         nperm = fact(natm-2)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a,natm-2,iflag)
            endif
            do j = 3,natm
               atm(1,j) = aatm(1,a(j))
               atm(2,j) = aatm(2,a(j))
            enddo
         enddo

         
      return
      end
