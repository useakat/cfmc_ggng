      subroutine check_atm(maxnatm,atm,max_nsch,natm,iflag)
C     ****************************************************
C     Check the atm corresponds to a peripheral channel
C     
C     By Yoshitaro Takaesu @KIAS JAN 14 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer maxnatm,natm,iflag,max_nsch
      integer atm(2,maxnatm)
C     LOCAL VARIABLES 
      integer i,j
      integer zero_flag,nzero,zero_atm1(natm),nnatm,nsch
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 1

C The first element always should be larger than the 2nd one.    
      do i = 3,maxnatm
         if (atm(1,i).lt.atm(2,i)) return
      enddo
         
C (0,0) atms should be placed at the end of the atm chain. 
      do i = 3,maxnatm-1
         if ((atm(1,i)+atm(2,i)).eq.0) then
            do j = i+1,maxnatm
               if ((atm(1,j)+atm(2,j)).ne.0) return
            enddo
         endif
      enddo

C count the number of non-zero atms
      natm = 0
      do i = 1,maxnatm
         if (atm(1,i).ne.0) natm = natm +1
      enddo

C nsch cut
      nsch = 0
      do i = 3,natm
         if (atm(1,i)*atm(2,i).ne.0) nsch = nsch +1
      enddo
      if (nsch.gt.max_nsch) return

C The specific positioning of s-type atms attatched to 1st or 2nd leg. 
      if (atm(2,1).eq.0) then
         if (atm(2,natm).eq.0) return
      endif
      if (atm(2,2).eq.0) then
         if (atm(2,3).eq.0) return
      endif
C atm attached to 2nd leg should be less than atm attached to 1st leg.
      if (atm(2,1)+atm(2,2).eq.0) then
         if (atm(1,3).gt.atm(1,natm)) return
      endif

C fix the order of intermediate legs.
      if (atm(2,1)*atm(2,2).ne.0) then
         do i = 3,natm-1
            do j = i+1,natm
               if (atm(1,i).gt.atm(1,j)) return
            enddo
         enddo
      elseif ((atm(2,1).eq.0).and.(atm(2,2).ne.0)) then
         do i = 3,natm-2
            do j = i+1,natm-1
               if (atm(1,i).gt.atm(1,j)) return
            enddo
         enddo
      elseif ((atm(2,1).ne.0).and.(atm(2,2).eq.0)) then
         do i = 4,natm-1
            do j = i+1,natm
               if (atm(1,i).gt.atm(1,j)) return
            enddo
         enddo
      elseif ((atm(2,1).eq.0).and.(atm(2,2).eq.0)) then
         do i = 4,natm-2
            do j = i+1,natm-1
               if (atm(1,i).gt.atm(1,j)) return
            enddo
         enddo
      endif

c$$$      nzero = 0
c$$$      do i = 3,maxnatm
c$$$         if (atm(2,i).eq.0) then
c$$$            nzero = nzero +1
c$$$            zero_atm1(nzero) = atm(1,i) 
c$$$         endif
c$$$      enddo
c$$$      do i = 1,nzero-1
c$$$         do j = i+1,nzero
c$$$            if (zero_atm1(i).lt.zero_atm1(j)) return
c$$$         enddo
c$$$      enddo

      iflag = 0
      
      return
      end
