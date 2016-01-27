      program find_perch
      implicitnone
C     CONSTANTS
C     GLOBAL VARIABLES
      include 'chinfo.inc'
C     ARGUMENTS
C     LOCAL VARIABLES 
      integer i,j
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14
      integer p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14
      integer natm,used(ngluons+1,4),nchannels,cfflag
      integer atm(2,max_natm),cf(ngluons,ngluons-1),ierr,zero_flag
      integer nonzero_channels
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nchannels = 0
      nonzero_channels = 0
      nch_per = 0
      nch_opt = 0
      ntch_per = 0
      ntch_opt = 0
      nsch_per = 0
      nsch_opt = 0
      n2sch_per = 0
      n2sch_opt = 0
      n3sch_per = 0
      n3sch_opt = 0
      n4sch_per = 0
      n4sch_opt = 0
      do i = 1,maxnch
         opt_ncf(i) = 0
      enddo
      do i = 1,nbcf
         cf_nopt(i) = 0
      enddo

      do i =1,ngluons
         call gen_basiccf(i,ngluons,cf(1,i))
      enddo

      do i1 = 3,ngluons+1
         do i = 3,ngluons+1
            used(i,1) = 0
         enddo
         if (i1.eq.ngluons+1) then
            p1 = 0
            used(i1,1) = used(i1,1) +1
            if (used(i1,1).gt.(ngluons-4)) cycle
         else
            if (used(i1,1).ne.0) then
               cycle
            else
               p1 = i1
               used(i1,1) = 1
            endif
         endif
      do i2 = 3,ngluons+1
         do i = 3,ngluons+1
            used(i,2) = used(i,1)
         enddo
         if (i2.eq.ngluons+1) then
            p2 = 0
            used(i2,2) = used(i2,2) +1
            if (used(i2,2).gt.(ngluons-4)) cycle
         else
            if (used(i2,2).ne.0) then
               cycle
            else
               p2 = i2
               used(i2,2) = 1
            endif
         endif
      do i3 = 3,ngluons+1
         do i = 3,ngluons+1
            used(i,3) = used(i,2)
         enddo
         if (i3.eq.ngluons+1) then
            p3 = 0
            used(i3,3) = used(i3,3) +1
            if (used(i3,3).gt.(ngluons-4)) cycle
         else
            if (used(i3,3).ne.0) then
               cycle
            else
               p3 = i3
               used(i3,3) = 1
            endif
         endif
      do i4 = 3,ngluons+1
         do i = 3,ngluons+1
            used(i,4) = used(i,3)
         enddo
         if (i4.eq.ngluons+1) then
            p4 = 0
            used(i4,4) = used(i4,4) +1
            if (used(i4,4).gt.(ngluons-4)) cycle
         else
            if (used(i4,4).ne.0) then
               cycle
            else
               p4 = i4
               used(i4,4) = 1
            endif
         endif
         do i = 1,max_natm
            do j = 1,2
               atm(j,i) = 0
            enddo
         enddo
         atm(1,1) = 1
         atm(2,1) = p1
         atm(1,2) = 2
         atm(2,2) = p2
         atm(1,3) = p3
         atm(2,3) = p4
         call analyse_atm(ngluons,cf,max_natm,atm
     &        ,nchannels,nonzero_channels,ierr)
         if (ierr.eq.1) cycle
      enddo
      enddo
      enddo
      enddo
      write(6,*) 
      write(6,*) "nchannels =", nchannels
      write(6,*) "non-zero channels =", nonzero_channels

      call write_chinfo

      end
