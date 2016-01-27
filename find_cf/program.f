      program program
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Nov.19 2012
C     ****************************************************
      implicitnone
C     CONSTANTS
C     GLOBAL VARIABLES
C     ARGUMENTS
      character*2 cngluons
      integer ngluons
C     LOCAL VARIABLES 
      integer i,j,i1,i2,i3,i4,i5,i6,p1,p2,p3,p4,p5,p6
      integer natm,used(10,6),nchannels,max_natm,cfflag,iflag(10000)
      integer atm(2,10),cflow(ngluons,10000),ncf,cf(ngluons,ngluons-1)
      integer iflag_all
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      call getarg(1,cngluons)
      read(cnext,*) ngluons

      max_natm = ngluons -2
      nchannels = 0

      if (ngluons.eq.5) then
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
              do j = 3,ngluons+1
                 used(j,2) = used(j,1)
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
                 do j = 3,ngluons+1
                    used(j,3) = used(j,2)
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
                    do j = 3,ngluons+1
                       used(j,4) = used(j,3)
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
                    if (p4.gt.p3) cycle
                    nchannels = nchannels +1
                    natm = 2
                    do i = 1,max_natm
                       do j = 1,2
                          atm(j,i) = 0
                       enddo
                    enddo
                    atm(1,1) = 1
                    atm(2,1) = p1
                    atm(1,2) = 2
                    atm(2,2) = p2
                    if ((p3+p4).ne.0) then
                       natm = natm +1
                       atm(1,natm) = p3
                       atm(2,natm) = p4
                    endif
                    call find_chcf(ngluons,natm,atm,cflow,ncf)
                    iflag_all = 0
                    do i = 1,ncf
                       do j = 1,ngluons-1
                          call check_cf(ngluons,cf(1,j),cflow(1,i),iflag(i))
                          if (iflag_all.eq.0) then
                             if (iflag(i).eq.1) iflag_all = 1
                          endif
                          if (iflag(i).eq.1) exit
                       enddo
                    enddo
                    if (iflag_all.eq.1) then
                       write(6,*)
                       write(6,*) nchannels,":",atm(1,1),atm(2,1)
     &                      ,atm(1,2),atm(2,2),atm(1,3),atm(2,3)
                       do i = 1,ncf
                          if (iflag(i).eq.1) then
                             write(6,*) "flow",i,": ",cflow(1,i),cflow(2,i),cflow(3,i)
     &                            ,cflow(4,i),cflow(5,i)
                          endif
                       enddo      
                    endif
                 enddo
              enddo
           enddo
        enddo
      elseif (ngluons.eq.6) then
        cf(1,1) = 1	
        cf(2,1) = 2	
        cf(3,1) = 3
        cf(4,1) = 4
        cf(5,1) = 5
        cf(6,1) = 6
        cf(1,2) = 1
        cf(2,2) = 3
        cf(3,2) = 2
        cf(4,2) = 4
        cf(5,2) = 5
        cf(6,2) = 6
        cf(1,3) = 1
        cf(2,3) = 3
        cf(3,3) = 4
        cf(4,3) = 2
        cf(5,3) = 5
        cf(6,3) = 6
        cf(1,4) = 1
        cf(2,4) = 3
        cf(3,4) = 4
        cf(4,4) = 5
        cf(5,4) = 2
        cf(6,4) = 6
        cf(1,5) = 1
        cf(2,5) = 3
        cf(3,5) = 4
        cf(4,5) = 5
        cf(5,5) = 6
        cf(6,5) = 2

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
              do j = 3,ngluons+1
                 used(j,2) = used(j,1)
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
                 do j = 3,ngluons+1
                    used(j,3) = used(j,2)
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
                    do j = 3,ngluons+1
                       used(j,4) = used(j,3)
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
                    do i5 = 3,ngluons+1
                       do j = 3,ngluons+1
                          used(j,5) = used(j,4)
                       enddo
                       if (i5.eq.ngluons+1) then
                          p5 = 0
                          used(i5,5) = used(i5,5) +1
                          if (used(i5,5).gt.(ngluons-4)) cycle
                       else
                          if (used(i5,5).ne.0) then
                             cycle
                          else
                             p5 = i5
                             used(i5,5) = 1
                          endif
                       endif
                       do i6 = 3,ngluons+1
                          do j = 3,ngluons+1
                             used(j,6) = used(j,5)
                          enddo
                          if (i6.eq.ngluons+1) then
                             p6 = 0
                             used(i6,6) = used(i6,6) +1
                             if (used(i6,6).gt.(ngluons-4)) cycle
                          else
                             if (used(i6,6).ne.0) then
                                cycle
                             else
                                p6 = i6
                                used(i6,6) = 1
                             endif
                          endif
                          if (p4.lt.p3) cycle
                          if (p5.gt.p3) cycle
                          if (p6.gt.p4) cycle
                          nchannels = nchannels +1
                          natm = 2
                          do i = 1,max_natm
                             do j = 1,2
                                atm(j,i) = 0
                             enddo
                          enddo
                          atm(1,1) = 1
                          atm(2,1) = p1
                          atm(1,2) = 2
                          atm(2,2) = p2
                          if ((p3+p5).ne.0) then
                             natm = natm +1
                             atm(1,natm) = p3
                             atm(2,natm) = p5
                          endif
                          if ((p4+p6).ne.0) then
                             natm = natm +1
                             atm(1,natm) = p4
                             atm(2,natm) = p6
                          endif
                          call find_chcf(ngluons,natm,atm,cflow,ncf)
                          iflag_all = 0
                          do i = 1,ncf
                             do j = 1,ngluons-1
                                call check_cf(ngluons,cf(1,j),cflow(1,i),iflag(i))
                                if (iflag_all.eq.0) then
                                   if (iflag(i).eq.1) iflag_all = 1
                                endif
                                if (iflag(i).eq.1) exit
                             enddo
                          enddo
                          if (iflag_all.eq.1) then
                             write(6,*)
                             write(6,*) nchannels,":",atm(1,1),atm(2,1)
     &                            ,atm(1,2),atm(2,2),atm(1,3),atm(2,3),atm(1,4),atm(2,4)
                             do i = 1,ncf
                                if (iflag(i).eq.1) then
                                   write(6,*) "flow",i,": ",cflow(1,i),cflow(2,i),cflow(3,i)
     &                                  ,cflow(4,i),cflow(5,i),cflow(6,i)
                                endif
                             enddo      
                          endif
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo     

      endif

      end
      
