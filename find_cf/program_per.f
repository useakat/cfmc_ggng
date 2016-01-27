      program program
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Nov.19 2012
C     ****************************************************
      implicitnone
C     CONSTANTS
C     GLOBAL VARIABLES
C     LOCAL VARIABLES 
      integer i,j,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12
      integer p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
      integer natm,ngluons,used(20,20),nchannels,max_natm,cfflag
      parameter (ngluons=9)
      parameter (max_natm=ngluons-2)
      integer atm(2,10),cf(ngluons,ngluons-1),ierr,zero_flag
      integer nonzero_channels
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nchannels = 0
      nonzero_channels = 0

      if (ngluons.eq.5) then
         call gen_basiccf(ngluons,cf)
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
                     natm = 2
                     do i = 1,max_natm
                        do j = 1,2
                           atm(j,i) = 0
                        enddo
                     enddo
                     zero_flag = 0
                     atm(1,1) = 1
                     atm(2,1) = p1
                     atm(1,2) = 2
                     atm(2,2) = p2
                     natm = natm +1
                     atm(1,natm) = p3
                     atm(2,natm) = p4
                     
                     call analyse_atm(ngluons,cf,natm,atm
     &                    ,nchannels,nonzero_channels,ierr)
                     if (ierr.eq.1) cycle
                     
                  enddo
               enddo
            enddo
         enddo
         write(6,*) 
         write(6,*) "nchannels =", nchannels
         write(6,*) "non-zero channels =", nonzero_channels
      elseif (ngluons.eq.6) then
         call gen_basiccf(ngluons,cf)

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
                           natm = natm +1
                           atm(1,natm) = p3
                           atm(2,natm) = p5
                           natm = natm +1
                           atm(1,natm) = p4
                           atm(2,natm) = p6
                           
                           call analyse_atm(ngluons,cf,natm,atm
     &                          ,nchannels,nonzero_channels,ierr)
                           if (ierr.eq.1) cycle

                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo     
         write(6,*) 
         write(6,*) "nchannels =", nchannels
         write(6,*) "non-zero channels =", nonzero_channels
      elseif (ngluons.eq.7) then
c      elseif (ngluons.eq.0) then
         call gen_basiccf(ngluons,cf)

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
                           do i7 = 3,ngluons+1
                              do j = 3,ngluons+1
                                 used(j,7) = used(j,6)
                              enddo
                              if (i7.eq.ngluons+1) then
                                 p7 = 0
                                 used(i7,7) = used(i7,7) +1
                                 if (used(i7,7).gt.(ngluons-4)) cycle
                              else
                                 if (used(i7,7).ne.0) then
                                    cycle
                                 else
                                    p7 = i7
                                    used(i7,7) = 1
                                 endif
                              endif                             
                              do i8 = 3,ngluons+1
                                 do j = 3,ngluons+1
                                    used(j,8) = used(j,7)
                                 enddo
                                 if (i8.eq.ngluons+1) then
                                    p8 = 0
                                    used(i8,8) = used(i8,8) +1
                                    if (used(i8,8).gt.(ngluons-4)) cycle
                                 else
                                    if (used(i8,8).ne.0) then
                                       cycle
                                    else
                                       p8 = i8
                                       used(i8,8) = 1
                                    endif
                                 endif                             
                                 
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
                                 natm = natm +1
                                 atm(1,natm) = p3
                                 atm(2,natm) = p6
                                 natm = natm +1
                                 atm(1,natm) = p4
                                 atm(2,natm) = p7
                                 natm = natm +1
                                 atm(1,natm) = p5
                                 atm(2,natm) = p8
                                 
                                 call analyse_atm(ngluons,cf,natm,atm
     &                                ,nchannels,nonzero_channels,ierr)
                                 if (ierr.eq.1) cycle
                                 
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo     
            enddo
         enddo
         write(6,*) 
         write(6,*) "nchannels =", nchannels
         write(6,*) "non-zero channels =", nonzero_channels
      elseif (ngluons.eq.8) then
         call gen_basiccf(ngluons,cf)

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
                           do i7 = 3,ngluons+1
                              do j = 3,ngluons+1
                                 used(j,7) = used(j,6)
                              enddo
                              if (i7.eq.ngluons+1) then
                                 p7 = 0
                                 used(i7,7) = used(i7,7) +1
                                 if (used(i7,7).gt.(ngluons-4)) cycle
                              else
                                 if (used(i7,7).ne.0) then
                                    cycle
                                 else
                                    p7 = i7
                                    used(i7,7) = 1
                                 endif
                              endif                             
                              do i8 = 3,ngluons+1
                                 do j = 3,ngluons+1
                                    used(j,8) = used(j,7)
                                 enddo
                                 if (i8.eq.ngluons+1) then
                                    p8 = 0
                                    used(i8,8) = used(i8,8) +1
                                    if (used(i8,8).gt.(ngluons-4)) cycle
                                 else
                                    if (used(i8,8).ne.0) then
                                       cycle
                                    else
                                       p8 = i8
                                       used(i8,8) = 1
                                    endif
                                 endif                             
                                 do i9 = 3,ngluons+1
                                    do j = 3,ngluons+1
                                       used(j,9) = used(j,8)
                                    enddo
                                    if (i9.eq.ngluons+1) then
                                       p9 = 0
                                       used(i9,9) = used(i9,9) +1
                                       if (used(i9,9).gt.(ngluons-4)) cycle
                                    else
                                       if (used(i9,9).ne.0) then
                                          cycle
                                       else
                                          p9 = i9
                                          used(i9,9) = 1
                                       endif
                                    endif                             
                                    do i10 = 3,ngluons+1
                                       do j = 3,ngluons+1
                                          used(j,10) = used(j,9)
                                       enddo
                                       if (i10.eq.ngluons+1) then
                                          p10 = 0
                                          used(i10,10) = used(i10,10) +1
                                          if (used(i10,10).gt.(ngluons-4)) cycle
                                       else
                                          if (used(i10,10).ne.0) then
                                             cycle
                                          else
                                             p10 = i10
                                             used(i10,10) = 1
                                          endif
                                       endif                             
                                       
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
                                       natm = natm +1
                                       atm(1,natm) = p3
                                       atm(2,natm) = p7
                                       natm = natm +1
                                       atm(1,natm) = p4
                                       atm(2,natm) = p8
                                       natm = natm +1
                                       atm(1,natm) = p5
                                       atm(2,natm) = p9
                                       natm = natm +1
                                       atm(1,natm) = p6
                                       atm(2,natm) = p10
                                       
                                       call analyse_atm(ngluons,cf,natm,atm
     &                                      ,nchannels,nonzero_channels,ierr)
                                       if (ierr.eq.1) cycle
                                       
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo     
                  enddo
               enddo
            enddo
         enddo
         write(6,*) 
         write(6,*) "nchannels =", nchannels
         write(6,*) "non-zero channels =", nonzero_channels
      elseif (ngluons.eq.9) then
c      elseif (ngluons.eq.0) then
         call gen_basiccf(ngluons,cf)

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
                           do i7 = 3,ngluons+1
                              do j = 3,ngluons+1
                                 used(j,7) = used(j,6)
                              enddo
                              if (i7.eq.ngluons+1) then
                                 p7 = 0
                                 used(i7,7) = used(i7,7) +1
                                 if (used(i7,7).gt.(ngluons-4)) cycle
                              else
                                 if (used(i7,7).ne.0) then
                                    cycle
                                 else
                                    p7 = i7
                                    used(i7,7) = 1
                                 endif
                              endif                             
                              do i8 = 3,ngluons+1
                                 do j = 3,ngluons+1
                                    used(j,8) = used(j,7)
                                 enddo
                                 if (i8.eq.ngluons+1) then
                                    p8 = 0
                                    used(i8,8) = used(i8,8) +1
                                    if (used(i8,8).gt.(ngluons-4)) cycle
                                 else
                                    if (used(i8,8).ne.0) then
                                       cycle
                                    else
                                       p8 = i8
                                       used(i8,8) = 1
                                    endif
                                 endif                             
                                 do i9 = 3,ngluons+1
                                    do j = 3,ngluons+1
                                       used(j,9) = used(j,8)
                                    enddo
                                    if (i9.eq.ngluons+1) then
                                       p9 = 0
                                       used(i9,9) = used(i9,9) +1
                                       if (used(i9,9).gt.(ngluons-4)) cycle
                                    else
                                       if (used(i9,9).ne.0) then
                                          cycle
                                       else
                                          p9 = i9
                                          used(i9,9) = 1
                                       endif
                                    endif                             
                                    do i10 = 3,ngluons+1
                                       do j = 3,ngluons+1
                                          used(j,10) = used(j,9)
                                       enddo
                                       if (i10.eq.ngluons+1) then
                                          p10 = 0
                                          used(i10,10) = used(i10,10) +1
                                          if (used(i10,10).gt.(ngluons-4)) cycle
                                       else
                                          if (used(i10,10).ne.0) then
                                             cycle
                                          else
                                             p10 = i10
                                             used(i10,10) = 1
                                          endif
                                       endif                             
                                 do i11 = 3,ngluons+1
                                    do j = 3,ngluons+1
                                       used(j,11) = used(j,10)
                                    enddo
                                    if (i11.eq.ngluons+1) then
                                       p11 = 0
                                       used(i11,11) = used(i11,11) +1
                                       if (used(i11,11).gt.(ngluons-4)) cycle
                                    else
                                       if (used(i11,11).ne.0) then
                                          cycle
                                       else
                                          p11 = i11
                                          used(i11,11) = 1
                                       endif
                                    endif                             
                                    do i12 = 3,ngluons+1
                                       do j = 3,ngluons+1
                                          used(j,12) = used(j,11)
                                       enddo
                                       if (i12.eq.ngluons+1) then
                                          p12 = 0
                                          used(i12,12) = used(i12,12) +1
                                          if (used(i12,12).gt.(ngluons-4)) cycle
                                       else
                                          if (used(i12,12).ne.0) then
                                             cycle
                                          else
                                             p12 = i12
                                             used(i12,12) = 1
                                          endif
                                       endif                             
                                       
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
                                       natm = natm +1
                                       atm(1,natm) = p3
                                       atm(2,natm) = p8
                                       natm = natm +1
                                       atm(1,natm) = p4
                                       atm(2,natm) = p9
                                       natm = natm +1
                                       atm(1,natm) = p5
                                       atm(2,natm) = p10
                                       natm = natm +1
                                       atm(1,natm) = p6
                                       atm(2,natm) = p11
                                       natm = natm +1
                                       atm(1,natm) = p7
                                       atm(2,natm) = p12
                                       
                                       call analyse_atm(ngluons,cf,natm,atm
     &                                      ,nchannels,nonzero_channels,ierr)
                                       if (ierr.eq.1) cycle
                                       
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo     
                  enddo
               enddo
            enddo
         enddo
         enddo
         enddo
         write(6,*) 
         write(6,*) "nchannels =", nchannels
         write(6,*) "non-zero channels =", nonzero_channels

c      elseif (ngluons.eq.9) then
      elseif (ngluons.eq.0) then
         call gen_basiccf(ngluons,cf)
         
c$$$         write(6,*) cf(1,1),cf(2,1),cf(3,1),cf(4,1),cf(5,1),cf(6,1),cf(7,1),cf(8,1),cf(9,1)
c$$$         write(6,*) cf(1,2),cf(2,2),cf(3,2),cf(4,2),cf(5,2),cf(6,2),cf(7,2),cf(8,2),cf(9,2)
c$$$         write(6,*) cf(1,3),cf(2,3),cf(3,3),cf(4,3),cf(5,3),cf(6,3),cf(7,3),cf(8,3),cf(9,3)
c$$$         write(6,*) cf(1,4),cf(2,4),cf(3,4),cf(4,4),cf(5,4),cf(6,4),cf(7,4),cf(8,4),cf(9,4)
c$$$         write(6,*) cf(1,5),cf(2,5),cf(3,5),cf(4,5),cf(5,5),cf(6,5),cf(7,5),cf(8,5),cf(9,5)
c$$$         write(6,*) cf(1,6),cf(2,6),cf(3,6),cf(4,6),cf(5,6),cf(6,6),cf(7,6),cf(8,6),cf(9,6)
c$$$         write(6,*) cf(1,7),cf(2,7),cf(3,7),cf(4,7),cf(5,7),cf(6,7),cf(7,7),cf(8,7),cf(9,7)
c$$$         write(6,*) cf(1,8),cf(2,8),cf(3,8),cf(4,8),cf(5,8),cf(6,8),cf(7,8),cf(8,8),cf(9,8)

         atm(1,1) = 1
         atm(2,1) = 3
         atm(1,2) = 2
         atm(2,2) = 4
         atm(1,3) = 9
         atm(2,3) = 0
         atm(1,4) = 8
         atm(2,4) = 0
         atm(1,5) = 7
         atm(2,5) = 0
         atm(1,6) = 6
         atm(2,6) = 0
         atm(1,7) = 5
         atm(2,7) = 0
         natm = 7

         call analyse_atm(ngluons,cf,natm,atm
     &        ,nchannels,nonzero_channels,ierr)
         
         write(6,*) 
         write(6,*) "nchannels =", nchannels
         write(6,*) "non-zero channels =", nonzero_channels         
      endif
      
      end
      
