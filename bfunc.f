      double precision function bfunc(z)
      implicit none
      include 'cparam.inc'
      include 'pshk.inc'
      include 'bfunc.inc'
      include 'lhe.inc'
      include 'opt_info.inc'
      real rand
      real*8 amp2

      bfunc = 0.d0
      flowtype = 1
      if (icf.ne.0) then
         ncf = fact( next -1 )
      endif
      call bfunc_ini

      posflag = 0
      delyflag = 0
      call GenMom_per(z,Q,ierr)
      if (ierr.ne.0) return

      npoints = npoints +1
      call impose_cuts(ierr)
      if (ierr.ne.0) return

      if (ipdf.eq.0) then
         do i = -5,5
            pdf1(i) = 1d0
            pdf2(i) = 1d0
         enddo
      elseif (ipdf.eq.1) then
         call pftopdg(1,x1,QQ,pdf1)
         call pftopdg(1,x2,QQ,pdf2)
      endif

c$$$      Q(0,1) = 5.94344757d0
c$$$      Q(1,1) = 0d0
c$$$      Q(2,1) = 0d0
c$$$      Q(3,1) = 5.94344757d0
c$$$      Q(0,2) = 270.42542d0
c$$$      Q(1,2) = 0d0
c$$$      Q(2,2) = 0d0
c$$$      Q(3,2) = -270.42542d0
c$$$      Q(0,3) = 100.434119d0
c$$$      Q(1,3) = 6.68537915d0
c$$$      Q(2,3) = 24.0872099d0
c$$$      Q(3,3) = -97.2734506d0
c$$$      Q(0,4) = 65.8901559d0
c$$$      Q(1,4) = 20.7754249d0
c$$$      Q(2,4) = -5.41193258d0
c$$$      Q(3,4) = -62.294505d0
c$$$      Q(0,5) = 110.044593d0
c$$$      Q(1,5) = -27.4607998d0
c$$$      Q(2,5) = -18.6752785d0
c$$$      Q(3,5) = -104.914017d0
      call get_amp2_cfz(Q,amp2,z(3*nfin-1))
c      call get_amp2(Q,amp2)

      if (icf.eq.0) then
         if (ichwgt.eq.0) call GetChWgt_SDE
         if (ichwgt.eq.1) call GetChWgt_per(Q)
      elseif (icf.eq.1) then
         call GetChWgt_per(Q)
c         call GetChWgt_opt(Q)
      elseif (icf.eq.2) then
c         call GetChWgt_per(Q)
         continue
      elseif (icf.eq.3) then
         if (ichwgt.eq.0) call GetChWgt_SDE
         if (ichwgt.eq.1) call GetChWgt_per(Q)
      endif

      if (icf.eq.0) then
         bfunc = pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
      elseif (icf.eq.1) then
         bfunc = ncf*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
      elseif (icf.eq.2) then
c         bfunc = cffact*fact(nfin)*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
c         bfunc = fact(nfin)*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
         bfunc = fact(nfin)*pdf1(0)*pdf2(0)*amp2*CH(ich)*jacob*wpsn
      elseif (icf.eq.3) then
         bfunc = fact(nfin)*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
c         call randperm_finmom(Q)
      endif

c      write(*,*) bfunc

c      if ((bfunc.lt.0).or.(bfunc.gt.1d9)) then
c         write(*,*) bfunc,amp2,jacob,wpsn,AAA,CH(ich)
c         call write_mom(next,Q,6)
c      endif

      if (1.eq.1) then
         if (abs(1d0-CH(ich)*jacob).gt.1d-3) then
            write(luerr,*) "CH*jacob .ne. 1"
            return
         else
            continue
         endif
      endif

      call fill_histo(z,Q,amp2,bfunc)

      return
      end
