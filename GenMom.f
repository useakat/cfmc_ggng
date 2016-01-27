      subroutine GenMom_per(rnd,Q,ierr)
      implicitnone
      include 'cparam.inc'
      include 'pshk.inc'
      include 'opt_info.inc'
      integer i,j
      integer ierr,sgncos,iperm,mom(next),ipos,ch_topo
      real*8 rnd(3*nfin-2),Q(0:3,next)
      real*8 q12(0:3),m12
      real*8 pmass,ppt,pR,pdeltar
      external pmass,ppt,pR,pdeltar

      ierr = 0
      do i = 1,next
         M(1,i) = mass(i) 
         M(2,i) = M(1,i)**2
         M(3,i) = 0
      enddo
      ch_topo = perch_type(ich) ! channel topology of this job
      sgncos = 1      
      ipos = 0

      if (posflag.eq.1) then
         if (ch_topo.eq.0) ipos = ipos_tch(ich)
      endif

      if (ch_topo.eq.0) then ! t-channel
         call wrap_pshkpos(ipos,rnd,sgncos,Q,ierr)
      elseif (ch_topo.eq.1) then ! 1s-channel
         call wrap_pshkpos2(rnd,sgncos,Q,ierr)
      elseif (ch_topo.eq.2) then ! 2s-channel
         call wrap_pshkpos2_2(rnd,sgncos,Q,ierr)
      elseif (ch_topo.eq.20) then ! ss-channel
c         call wrap_pshkpos2_2_0(rnd,sgncos,Q,ierr)
         call wrap_pshkpos2_2(rnd,sgncos,Q,ierr)
      elseif (ch_topo.eq.3) then ! 3s-channel
         call wrap_pshkpos2_3(rnd,sgncos,Q,ierr)
      elseif (ch_topo.eq.4) then ! 4s-channel
         call wrap_pshkpos2_4(rnd,sgncos,Q,ierr)
      endif
      if (ierr.ne.0) return

c$$$      q12(0) = Q(0,4)+Q(0,5)
c$$$      q12(1) = Q(1,4)+Q(1,5)
c$$$      q12(2) = Q(2,4)+Q(2,5)
c$$$      q12(3) = Q(3,4)+Q(3,5)
c$$$      m12 = pmass(q12)
c$$$
c$$$      if (m12.lt.2*ptcut) then
c$$$         if ((ppt(Q(0,4)).lt.ptcut).or.(ppt(Q(0,5)).lt.ptcut)) then
c$$$            nptcut_fails = nptcut_fails +1
c$$$         endif
c$$$         if (pdeltar(Q(0,4),Q(0,5)).lt.rcut) then
c$$$            nrcut_fails = nrcut_fails +1
c$$$         endif
c$$$      endif

      if (posflag.eq.0) then
         mom(1) = 1 
         mom(2) = 2
         do j = 1,nfin
            mom(2+j) = chmom(j,ich)
         enddo
         call mom_perm(next,mom,Q)
      endif

      return
      end
