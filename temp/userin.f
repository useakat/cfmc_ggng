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
      include 'bsffcm.inc'
      include 'opt_info.inc'
      integer i,idim,ipdfset,ch_type
      real*8 ncall_def,ncall_f(0:chtype_max),ncall_sch,ncall_2sch
      external bfunc

      pi = dacos(-1d0)      
      call hminit
CCCCCCC parameters CCCCCCCCCC
      EBMUP(1) = 7000d0
      EBMUP(2) = 7000d0
      PTcut   = 20d0
      AYEcut  = 5d0
      rcut    = 0.4d0
      lptcut = 0d0
      ipdfset = 4
      QQ = 20d0
c      alfas = 0.13d0
      alfas = 0.171202109d0
c      Qmax = 8d0
      nsample_cf = 5
      ncall_def = 10000
      ncall_f(0) = 1d0  ! tch
      ncall_f(1) = 2d0  ! sch
      ncall_f(20) = 4d0  ! ssch
      jfact(0) = 1d0
      jfact(1) = 1d0
      jfact(20) = 1d0
c      mgluon = 91.188d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ch_type = perch_type(ich)
      NCALL = ncall_f(ch_type)*ncall_def

      Qcut = dsqrt( 2*ptcut**2*( 1d0 -dcos(rcut) ) )
      Qmax = EBMUP(1) +EBMUP(2)
      do i = 1,next
         mass(i) = 0d0
      enddo
      G = dsqrt(4*pi*alfas)
      GC_2 = (0d0,0.205302511496246d0)
      GC_3 = (0d0,-0.307953767244369d0)
      GC_4 = -G
      GC_5 = ci*G
      GC_6 = ci*G**2
      PTcut2  = PTcut
      AYEcut2 = AYEcut
      AYEcut3 = AYEcut
      sqrts = EBMUP(1) +EBMUP(2)
      AQCDUP = alfas
      AQEDUP = 1d0/(4*pi)
      SCALUP = QQ
      call setCtq6(ipdfset)
      PDFGUP(1) = 0
      PDFSUP(1) = 10042
      PDFGUP(2) = 0
      PDFSUP(2) = 10042
      ss = sqrts**2
      ndim = 3*nfin -4 +ipdf*2 
      nwild = ndim
      do idim=1,ndim
         xl(idim) = 0.d0
         xu(idim) = 1.d0
         ig(idim) = 1
      enddo

      call bssetp(ncall,itmx1,itmx2,acc1,acc2)

      open(1,file="params.dat",status="unknown")
      write(1,*) 'Process:',' g g >',nfin,'g'
      write(1,*) ''
      write(1,*) 'mass:',mass
      write(1,*) ''
      write(1,*) 'EBMUP(1):',EBMUP(1)
      write(1,*) 'EBMUP(2):',EBMUP(2)
      write(1,*) ''
      write(1,*) 'pT Cut:',ptcut,'GeV'
      write(1,*) 'eta Cut:',AYEcut
      write(1,*) 'dR Cut:',rcut
      write(1,*) 'Leading pT Cut:',lptcut,'GeV'
      write(1,*) ''
      if (ipdf.eq.0) write(1,*) 'No PDF'
      if (ipdf.eq.1) then
         if (ipdfset.eq.4) write(1,*)'PDF: CTEQ6L1'
      endif
      write(1,*) 'alpha_s:',alfas
      write(1,*) 'Fact Scale:',QQ,'GeV'
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
      if (ich_type.eq.0) write(1,*) 'All Peripheral Channels'
      if (ich_type.eq.1) write(1,*) 'Peripheral tch'
      if (ich_type.eq.2) write(1,*) 'Peripheral tch + sch'
      if ((ich_type.eq.0).or.(ich_type.ge.2)) write(1,*) 
     &     'ncall factor (tch): 1.'
      if ((ich_type.eq.0).or.(ich_type.ge.2)) write(1,*) 
     &     'ncall factor (sch):',ncall_sch
      if ((ich_type.eq.0).or.(ich_type.ge.3)) write(1,*) 
     &     'ncall factor (sch):',ncall_2sch
      write(1,*) 'Channels:',nch  
      if (ichwgt.eq.0) write(1,*) 'Channel Weights:',nwgt,
     &' (Diagram wgt)'  
      if (ichwgt.eq.1) write(1,*) 'Channel Weights:',nwgt,
     &' (Jacob wgt)'  
      write(1,*) ''
      write(1,*) 'Phase Space points / iteration:', ncall_def
      write(1,*) ''
      close(1)

      call bssetd(ndim,nwild,xl,xu,ig)
      call set_histo(xl,xu)

      return
      end
