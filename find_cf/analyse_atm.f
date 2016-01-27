      subroutine analyse_atm(ngluons,cf,maxnatm,atm,nchannels
     &     ,nonzero_channels,ierr)
C     ****************************************************
C     Make peripheral channel vs. basic cf table from atm
C     generated naively.
C
C     By Yoshitaro Takaesu @KIAS JAN 15 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer ngluons,maxnatm,cf(ngluons,ngluons-1),ierr
      integer atm(2,maxnatm)
C     LOCAL VARIABLES 
      integer i,j
      integer nchannels,cflow(ngluons,100000),ncf,nonzero_cf,max_nsch
      integer iflag_all,cfused(ngluons-1),iflag(100000),nonzero_channels
      integer natm
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      max_nsch = 10
      call check_atm(maxnatm,atm,max_nsch,natm,ierr)
      if (ierr.eq.1) return

      nchannels = nchannels +1
      call find_perchcf(ngluons,natm,atm,cflow,ncf)
      iflag_all = 0
      do i = 1,ngluons-1
         cfused(i) = 0
      enddo
      do i = 1,ncf
         do j = 1,ngluons-1
            call check_cf(ngluons,cf(1,j),cflow(1,i),iflag(i))
            if (iflag(i).eq.1) then
               if (cfused(j).eq.0) then
                  iflag_all = 1
                  cfused(j) = 1
               else
                  iflag(i) = 0
               endif
               exit
            endif
         enddo
      enddo
      if (iflag_all.eq.1) then
         nonzero_channels = nonzero_channels +1
         write(6,*)
         write(6,*) nonzero_channels,":",atm(1,1),atm(2,1)
     &        ,atm(1,2),atm(2,2),atm(1,3),atm(2,3),atm(1,4),atm(2,4)
     &        ,atm(1,5),atm(2,5),atm(1,6),atm(2,6),atm(1,7),atm(2,7)
     &        ,atm(1,8),atm(2,8)
         nonzero_cf = 0
         do i = 1,ncf
            if (iflag(i).eq.1) then
               nonzero_cf = nonzero_cf +1
               write(6,*) "flow",nonzero_cf,": ",cflow(1,i),cflow(2,i),cflow(3,i)
     &              ,cflow(4,i),cflow(5,i),cflow(6,i),cflow(7,i),cflow(8,i)
     &              ,cflow(9,i),cflow(10,i)
            endif
         enddo      
      endif

      call analyse_perch(natm,atm,iflag_all)
      if (iflag_all.eq.1) then
         call analyse_optch(cfused)
      endif

      return
      end
