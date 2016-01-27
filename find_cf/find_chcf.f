      subroutine find_chcf(next,natm,atm,cflow,ncf)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Nov.19 2012
C     ****************************************************
      implicitnone
C     
C     CONSTANTS
C     
C     
C     ARGUMENTS 
C     
      integer next,natm
      integer atm(2,natm),cflow(next,100000),ncf
C     
C     GLOBAL VARIABLES
C     
C     
C     LOCAL VARIABLES 
C     
      integer i,j,k
      integer iflag,nperm,pos2,ierr,a(natm),nncf,maxk(natm)
      integer cflowtmp
C     
C     EXTERNAL FUNCTIONS
C     
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 0
      do i = 1,natm
         a(i) = i
      enddo
      do i = 1,natm
         maxk(i) = 1
         if (atm(2,i).ne.0) maxk(i) = 2
      enddo
       
      ncf = 0
      nperm = fact(natm-1)
      do i = 1,nperm
         if (i.ne.1) then
            call ipnext(a(2),natm-1,iflag)
         endif
         do j = 2,natm
            if (a(j).eq.2) pos2 = j
         enddo
         ierr = 0
         do j = 2,pos2-2
            do k = j+1,pos2-1
               if (a(j).lt.a(k)) ierr = 1
            enddo
         enddo
         if (ierr.ne.0) cycle
         ierr = 0
         do j = pos2+1,natm-1
            do k = j+1,natm
               if (a(j).gt.a(k)) ierr = 1
            enddo
         enddo
         if (ierr.ne.0) cycle

         call get_cflow(next,natm,atm,a,maxk,nncf,cflow(1,ncf+1))
         ncf = ncf +nncf
      enddo

      do i = 1,ncf
         if (cflow(1,i).ne.1) then
            cflowtmp = cflow(1,i)
            do j = 1,next-1
               cflow(j,i) = cflow(j+1,i)
            enddo
            cflow(next,i) = cflowtmp
         endif
      enddo

      return
      end


      subroutine find_perchcf(next,natm,atm,cflow,ncf)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN.14 2013
C     ****************************************************
      implicitnone
C     CONSTANTS
C     ARGUMENTS 
      integer next,natm
      integer atm(2,natm),cflow(next,100000),ncf
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i,j,k
      integer iflag,nperm,a(natm),nonzero,zero,nonzero_atm(2,natm)
      integer zero_atm(2,natm),incf,iatm,aatm(2,natm),aaatm(2,natm)
C     EXTERNAL FUNCTIONS
      integer fact
      external fact
C     ----------
C     BEGIN CODE
C     ----------
      iflag = 0
      ncf = 0
      do i = 3,natm
         a(i) = i
      enddo
      do i = 1,natm
         aatm(1,i) = atm(1,i)
         aatm(2,i) = atm(2,i)
      enddo

      if ((atm(2,1)*atm(2,2)).ne.0) then
         nperm = fact(natm-2)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(3),natm-2,iflag)
            endif
            do j = 3,natm
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      elseif ((atm(2,1).eq.0).and.(atm(2,2).ne.0)) then
         nperm = fact(natm-3)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(3),natm-3,iflag)
            endif
            do j = 3,natm-1
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      elseif ((atm(2,1).ne.0).and.(atm(2,2).eq.0)) then
         nperm = fact(natm-3)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(4),natm-3,iflag)
            endif
            do j = 4,natm
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      elseif ((atm(2,1).eq.0).and.(atm(2,2).eq.0)) then
         nperm = fact(natm-4)
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(4),natm-4,iflag)
            endif
            do j = 4,natm-1
               aatm(1,j) = atm(1,a(j))
               aatm(2,j) = atm(2,a(j))
            enddo
            call find_chcf(next,natm,aatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
         do i = 1,natm
            aatm(1,i) = atm(1,i)
            aatm(2,i) = atm(2,i)
         enddo
         aatm(1,3) = atm(1,natm)
         aatm(2,3) = atm(2,natm)
         aatm(1,natm) = atm(1,3)
         aatm(2,natm) = atm(2,3)

         do i = 1,natm
            aaatm(1,i) = aatm(1,i)
            aaatm(2,i) = aatm(2,i)
         enddo
         iflag = 0
         do i = 1,natm
            a(i) = i
         enddo
         do i = 1,nperm
            if (i.ne.1) then
               call ipnext(a(4),natm-4,iflag)
            endif
            do j = 4,natm-1
               aaatm(1,j) = aatm(1,a(j))
               aaatm(2,j) = aatm(2,a(j))
            enddo
            call find_chcf(next,natm,aaatm,cflow(1,ncf+1),incf)
            ncf = ncf +incf
         enddo
      endif

      return
      end
