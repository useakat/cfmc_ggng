      subroutine userin
      implicit none
      integer maxdim
      parameter (maxdim=50)
      real*8 xl(maxdim),xu(maxdim)
      integer ig(maxdim)
      include 'cparam.inc'
      include 'coupl.inc'
      include 'pshk.inc'
      include 'lhe.inc'
      include 'opt_info.inc'
      include 'alfas.inc'
      integer i,idim,ipdfset,ch_type,ialfas
      real*8 ncall_def,ncall_f(0:chtype_max),ncall_sch,ncall_2sch
      real*8 Qcut_def,alphas,ren_scale
      integer ncall,acc1,acc2,itmx1,itmx2
      external bfunc,alphas

      pi = dacos(-1d0)      
CCCCCCC parameters CCCCCCCCCC
      EBMUP(1) = 7000d0
      EBMUP(2) = 7000d0

      PTcut   = 10d0
      AYEcut  = 3.5d0
      rcut    = 0.3d0
      ptijcut = 0d0
      lptcut = 0d0
      Qcut_def = 1d0

      ipdfset = 4 ! 3:CTEQ6L, 4:CTEQ6L1
      QQ = 20d0
      ren_scale = 20d0
      ialfas = 1   ! 0:alfas by hand, 1:alfas @ ren_scale 
      alfas = 0.1712021d0

      nsample_cf = 1
      ncall_def = 10000
c      ncall_def = 4000
c      ncall_def = 1024*2
      ncall_def = 128*2
      acc1 = 0.1d0
      acc2 = 0.1d0
      itmx1 = 3
      itmx2 = 1

      ncall_f(0) = 1d0  ! tch
      ncall_f(1) = 2d0  ! sch
      ncall_f(2) = 2d0  ! 2sch
      ncall_f(20) = 2d0  ! ssch
      ncall_f(3) = 2d0  ! 3sch
      jfact(0) = 1d0
      jfact(1) = 1d0
      jfact(2) = 1d0
      jfact(20) = 1d0
      jfact(3) = 1d0
c      mgluon = 91.188d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ch_type = perch_type(ich)
      NCALL = ncall_f(ch_type)*ncall_def

      if (rcut.gt.0d0) then
         Qcut = dsqrt( 2*ptcut**2*( 1d0 -dcos(rcut) ) )
      else
         Qcut = Qcut_def
      endif
      Qmax = EBMUP(1) +EBMUP(2)

      do i = 1,next
         mass(i) = 0d0
      enddo
      if (ialfas.eq.1) then
c         asmz=0.118d0 
         asmz=0.129783d0 ! compatible with Sherpa CTEQ6L1
c         asmz=0.13d0 ! compatible with MG5 CTEQ6L1 
         nloop=1 ! for CTEQ6L1
         alfas = alphas(ren_scale)
      endif
      call setCtq6(ipdfset)
      G = dsqrt(4*pi*alfas)
      GC_2 = (0d0,0.205302511496246d0)
      GC_3 = (0d0,-0.307953767244369d0)
      GC_4 = -G
      GC_5 = ci*G
      GC_6 = ci*G**2
      PTcut2  = PTcut
      do i = 1,nfin
         ptmin_cut(i) = ptcut
         ptmax_cut(i) = Qmax/2d0
      enddo
      AYEcut2 = AYEcut
      AYEcut3 = AYEcut
      sqrts = EBMUP(1) +EBMUP(2)
      AQCDUP = alfas
      AQEDUP = 1d0/(4*pi)
      SCALUP = QQ
      PDFGUP(1) = 0
      PDFSUP(1) = 10042
      PDFGUP(2) = 0
      PDFSUP(2) = 10042
      ss = sqrts**2
      if (icf.eq.1) then
         ndim = 3*nfin -4 +ipdf*2 +1
c         ndim = 3*nfin -4 +ipdf*2
      else
         ndim = 3*nfin -4 +ipdf*2
      endif
      nwild = min(ndim,16)
      if ( nfin.eq.6) then
         do idim=1,nwild
            xl(idim) = 0d0
            xu(idim) = 1d0
            ig(idim) = 1
         enddo
         do idim=nwild+1,ndim
            xl(idim) = 0d0
            xu(idim) = 1d0
            ig(idim) = 1
         enddo
      else
         do idim=1,ndim
            xl(idim) = 0d0
            xu(idim) = 1d0
            ig(idim) = 1
         enddo
      endif

      call bssetp(ncall,itmx1,itmx2,acc1,acc2)
      call bssetd(ndim,nwild,xl,xu,ig)
      call set_histo(xl,xu)

      open(1,file="params.dat",status="unknown")
      write(1,*) 'Process:',' g g >',nfin,'g'
      write(1,*) ''
      write(1,*) 'mass:',mass
      write(1,*) ''
      write(1,*) 'EBMUP(1):',EBMUP(1)
      write(1,*) 'EBMUP(2):',EBMUP(2)
      write(1,*) ''
      if (ptcut.gt.0) write(1,*) 'pT Cut:',ptcut,'GeV'
      if (ayecut.gt.0) write(1,*) 'eta Cut:',AYEcut
      if (rcut.gt.0) write(1,*) 'dR Cut:',rcut
      if (ptijcut.gt.0) write(1,*) 'pTij Cut:',ptijcut
      if (lptcut.gt.0) write(1,*) 'Leading pT Cut:',lptcut,'GeV'
      write(1,*) ''
      if (ipdf.eq.0) write(1,*) 'No PDF'
      if (ipdf.eq.1) then
         if (ipdfset.eq.3) write(1,*)'PDF: CTEQ6L'
         if (ipdfset.eq.4) write(1,*)'PDF: CTEQ6L1'
      endif
      write(1,*) 'alpha_s:',alfas
      write(1,*) 'Renormalization Scale:',ren_scale,'GeV'
      write(1,*) 'Factorization Scale:',QQ,'GeV'
      write(1,*) 'inv. mass minCut(gen):',Qcut,'GeV'
      write(1,*) ''
      write(1,*) 'GC_2:',GC_2     
      write(1,*) 'GC_3:',GC_3
      write(1,*) 'GC_4:',GC_4          
      write(1,*) 'GC_5:',GC_5     
      write(1,*) 'GC_6:',GC_6     
      write(1,*) ''
      write(1,*) 'Phase Space Variables:',ndim
      write(1,*) 'Wild Variables:',nwild  
      write(1,*) ''
      if (icf.eq.0) write(1,*) 'Exact CF Sum'
      if (icf.eq.1) write(1,*) 'Simple CF Sampling: ',nsample_cf,
     &     ' CF/PS'
      if (icf.eq.2) write(1,*) 'Optimized CF Sampling: 1 CF/PS'
      if (icf.eq.3) write(1,*) 'Basic CF Sum'
      write(1,*) ''      
      if (ich_scheme.eq.0) write(1,*) 'Semi-peripheral channeling'
      if (ich_scheme.eq.1) write(1,*) 'Optimized channeling'
      write(1,*)
      if (ich_type.eq.0) write(1,*) 'All Peripheral Channels'
      if (ich_type.eq.1) write(1,*) 'Peripheral tch'
      if (ich_type.eq.2) write(1,*) 'Peripheral tch + sch'
      if (ich_type.eq.3) write(1,*) 'Peripheral tch + sch + 2sch'
      if ((ich_type.eq.0).or.(ich_type.ge.2)) write(1,*) 
     &     'ncall factor (tch):',ncall_f(0)
      if ((ich_type.eq.0).or.(ich_type.ge.2)) write(1,*) 
     &     'ncall factor (sch):',ncall_f(1)
      if ((ich_type.eq.0).or.(ich_type.ge.3)) write(1,*) 
     &     'ncall factor (2sch):',ncall_f(2)
      write(1,*) 'Channels:',nch  
      if (ichwgt.eq.0) write(1,*) 'Channel Weights:',nwgt,
     &' (Diagram wgt)'  
      if (ichwgt.eq.1) write(1,*) 'Channel Weights:',nwgt,
     &' (Jacob wgt)'  
      write(1,*) ''
      write(1,*) 'Phase Space points / iteration:', ncall_def
      write(1,*) ''
      close(1)

      return
      end
