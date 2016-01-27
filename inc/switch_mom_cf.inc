C -*- fortran -*-
         do i = 3,next
            do j = 0,3
               Qtemp(j,i) = Q(j,i)
            enddo
         enddo
         if (ich.eq.1) then
            do i = 0,3
               Q(i,3) = Qtemp(i,3)
               Q(i,4) = Qtemp(i,5)
               Q(i,5) = Qtemp(i,4)
            enddo
         elseif (ich.eq.2) then
            do i = 0,3
               Q(i,3) = Qtemp(i,3)
               Q(i,4) = Qtemp(i,4)
               Q(i,5) = Qtemp(i,5)
            enddo
         elseif (ich.eq.3) then
            do i = 0,3
               Q(i,3) = Qtemp(i,5)
               Q(i,4) = Qtemp(i,3)
               Q(i,5) = Qtemp(i,4)
            enddo
c            call switch_mom(Q(0,3),Q(0,4))
c            call switch_mom(Q(0,3),Q(0,5))
         elseif (ich.eq.4) then
            do i = 0,3
               Q(i,3) = Qtemp(i,4)
               Q(i,4) = Qtemp(i,3)
               Q(i,5) = Qtemp(i,5)
            enddo
c            call switch_mom(Q(0,3),Q(0,4))
         elseif (ich.eq.5) then
            do i = 0,3
               Q(i,3) = Qtemp(i,5)
               Q(i,4) = Qtemp(i,4)
               Q(i,5) = Qtemp(i,3)
            enddo
c            call switch_mom(Q(0,3),Q(0,5))
         elseif (ich.eq.6) then
            do i = 0,3
               Q(i,3) = Qtemp(i,4)
               Q(i,4) = Qtemp(i,5)
               Q(i,5) = Qtemp(i,3)
            enddo
c            call switch_mom(Q(0,3),Q(0,5))
c            call switch_mom(Q(0,3),Q(0,4))
         elseif (ich.eq.7) then
            do i = 0,3
               Q(i,3) = Qtemp(i,5)
               Q(i,4) = Qtemp(i,3)
               Q(i,5) = Qtemp(i,4)
            enddo
c            call switch_mom(Q(0,3),Q(0,5))
         elseif (ich.eq.8) then
            do i = 0,3
               Q(i,3) = Qtemp(i,3)
               Q(i,4) = Qtemp(i,5)
               Q(i,5) = Qtemp(i,4)
            enddo
c            call switch_mom(Q(0,4),Q(0,5))
         elseif (ich.eq.9) then
            do i = 0,3
               Q(i,3) = Qtemp(i,3)
               Q(i,4) = Qtemp(i,4)
               Q(i,5) = Qtemp(i,5)
            enddo
c            continue
         elseif (ich.eq.10) then
            do i = 0,3
               Q(i,3) = Qtemp(i,5)
               Q(i,4) = Qtemp(i,3)
               Q(i,5) = Qtemp(i,4)
            enddo
c            call switch_mom(Q(0,3),Q(0,5))
         elseif (ich.eq.11) then
            do i = 0,3
               Q(i,3) = Qtemp(i,3)
               Q(i,4) = Qtemp(i,5)
               Q(i,5) = Qtemp(i,4)
            enddo
c            call switch_mom(Q(0,4),Q(0,5))
         elseif (ich.eq.12) then
            do i = 0,3
               Q(i,3) = Qtemp(i,3)
               Q(i,4) = Qtemp(i,4)
               Q(i,5) = Qtemp(i,5)
            enddo
c            continue
         endif
