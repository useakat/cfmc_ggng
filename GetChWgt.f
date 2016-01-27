      subroutine GetChWgt_SDE
      implicitnone
      include 'cparam.inc'
      include 'opt_info.inc'
      include 'icf_ch.inc'
      integer i,iperch,ioptch

      AAA = 0d0
      do i = 1,nch
         if (ich_scheme.eq.0) then
            if (ich_type.eq.0) then 
               iperch = i
            elseif (ich_type.eq.1) then 
               iperch = icf_perch1(i)
            elseif (ich_type.eq.2) then 
               iperch = icf_perch2(i)
            elseif (ich_type.eq.3) then 
               iperch = icf_perch3(i)
            elseif (ich_type.eq.4) then 
               iperch = icf_perch4(i)
            endif
         elseif (ich_scheme.eq.1) then
            if (ich_type.eq.0) then 
               iperch = opt_per(i)
            elseif (ich_type.eq.1) then 
               iperch = opt_per(icf_optch1(i))
            elseif (ich_type.eq.2) then 
               iperch = opt_per(icf_optch2(i))
            elseif (ich_type.eq.3) then 
               iperch = opt_per(icf_optch3(i))
            elseif (ich_type.eq.4) then 
               iperch = opt_per(icf_optch4(i))
            endif
         endif
         AAA = AAA +AA(iperch)
      enddo
      CH(ich) = AA(ich)

      return
      end


      subroutine GetChWgt_per(P)
      implicitnone
      include 'cparam.inc'
      include 'icf_ch.inc'
      integer i,iperch
      real*8 P(0:3,next)

      call GenWgt_per(P)  ! Generate weights of all the channel
      call GetChWgt_SDE
c$$$      AAA = 0d0
c$$$      do i = 1,nch
c$$$         if (ich_type.eq.0) then
c$$$            AAA = AAA +AA(i)
c$$$         elseif (ich_type.eq.1) then
c$$$            iperch = icf_perch1(i)
c$$$            AAA = AAA +AA(iperch)
c$$$         endif
c$$$      enddo
c$$$      CH(ich) = AA(ich)

      return
      end


      subroutine GetChWgt_opt(P,icf_selected)
c
c     Calculate the sum of the singularity weights related to the CF(icf_selected)
c     in the optimized Leading color Sum scheme, and the channel weight
c
      implicitnone
      include 'cparam.inc'
      include 'opt_info.inc'
      real*8 P(0:3,next)
      integer i,ch_topo,ioptch,iperch,icf_selected

      call GenWgt_per(P) ! generate weights for all channels
      AAA = 0d0
      do i = 1,cf_nopt(icf_selected) ! do over the opt channels which are related to the CF
         ioptch = cf_opt(i,icf_selected) ! get the opt-ch number of the i_th opt-ch related to the CF 
         iperch = opt_per(ioptch) ! opt-ch number > ch number
         ch_topo = optch_type(ioptch) ! get the channel topology of the opt-ch
         if (ich_type.lt.4) then
            if (ich_type.eq.0) then
               AAA = AAA +AA(iperch)
            endif
            if (ich_type.gt.0) then
               if (ch_topo.eq.0) AAA = AAA +AA(iperch)
            endif
            if (ich_type.gt.1) then
               if (ch_topo.eq.1) AAA = AAA +AA(iperch)
            endif
            if (ich_type.gt.2) then
               if (ch_topo.eq.2) AAA = AAA +AA(iperch)
               if (ch_topo.eq.20) AAA = AAA +AA(iperch)
            endif
            if (ich_type.gt.3) then
               if (ch_topo.eq.3) AAA = AAA +AA(iperch)
            endif
         else
            if (ich_type.eq.4) then
               if (ch_topo.eq.1) AAA = AAA +AA(iperch)
            endif
         endif
      enddo
      CH(ich) = AA(ich)

      return
      end
