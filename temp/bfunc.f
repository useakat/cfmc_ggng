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

      call get_amp2(Q,amp2)

      if (icf.le.1) then
         if (ichwgt.eq.0) call GetChWgt_SDE
         if (ichwgt.eq.1) call GetChWgt_per(Q)
      elseif (icf.eq.2) then
         call GetChWgt_opt(Q)
      elseif (icf.eq.3) then
         if (ichwgt.eq.0) call GetChWgt_SDE
         if (ichwgt.eq.1) call GetChWgt_per(Q)
      endif

      if (icf.eq.0) then
         bfunc = pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
      elseif (icf.eq.1) then
         bfunc = ncf*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
      elseif (icf.eq.2) then
         bfunc = cffact*fact(nfin)*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
         call randperm_finmom(Q)
      elseif (icf.eq.3) then
         bfunc = fact(nfin)*pdf1(0)*pdf2(0)*amp2/AAA*CH(ich)*jacob*wpsn
         call randperm_finmom(Q)
      endif

      if (0.eq.1) then
         if (abs(1d0-CH(ich)*jacob).gt.1d-3) then
            write(99,*) "CH*jacob .ne. 1"
            return
         endif
      endif
c      write(6,*) bfunc,cffact,amp2,AAA,CH(ich)

      call fill_histo(z,Q,amp2,bfunc)

      return
      end
