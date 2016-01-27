C -*- fortran -*-
         if (ich.eq.1) then
            iperm = int(6*rand())+1
            if (iperm.eq.2) then
               call switch_mom(Q(0,4),Q(0,5))
            elseif (iperm.eq.3) then
               call switch_mom(Q(0,3),Q(0,4))
            elseif (iperm.eq.4) then
               call switch_mom(Q(0,3),Q(0,4))
               call switch_mom(Q(0,3),Q(0,5))
            elseif (iperm.eq.5) then
               call switch_mom(Q(0,4),Q(0,5))
               call switch_mom(Q(0,3),Q(0,5))
            elseif (iperm.eq.6) then
               call switch_mom(Q(0,3),Q(0,5))
            endif
         elseif (ich.eq.2) then
            iperm = int(3*rand())+1
            if (iperm.eq.2) then
               call switch_mom(Q(0,3),Q(0,5))
               call switch_mom(Q(0,4),Q(0,5))
            elseif (iperm.eq.3) then
               call switch_mom(Q(0,4),Q(0,5))
            endif
         elseif (ich.eq.3) then
            iperm = int(3*rand())+1
            if (iperm.eq.2) then
               call switch_mom(Q(0,3),Q(0,5))
               call switch_mom(Q(0,4),Q(0,5))
            elseif (iperm.eq.3) then
               call switch_mom(Q(0,4),Q(0,5))
            endif
         endif