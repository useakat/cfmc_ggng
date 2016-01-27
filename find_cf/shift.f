      subroutine shift(ia,in,ipos)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C
C     
C     ****************************************************
      implicit none

C     GLOBAL VARIABLES

C     CONSTANTS

C     ARGUMENTS 
      integer in
      integer ia(in),ipos
C     LOCAL VARIABLES 
      integer i,ib(in)
C     EXTERNAL FUNCTIONS

C     ----------
C     BEGIN CODE
C     ----------
      do i = 1,in
         ib(i) = ia(i)
      enddo
      
      ia(1) = ib(ipos)
      do i = 2,in
         if (i.le.ipos) then
            ia(i) = ib(i-1)
         elseif (i.gt.ipos) then
            ia(i) = ib(i)
         endif
      enddo

      return
      end
