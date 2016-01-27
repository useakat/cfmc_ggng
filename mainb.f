      program main
      implicit none

      double precision bfunc
      external bfunc
      integer i,seed
      character*4 iseed,iii
      character*1 genmode
      character*10 filename1
      character*11 filename2
      character*12 filename3
      integer it1,it2
      double precision estim,sigma,ctime

      include 'lhe.inc'
      include 'cparam.inc'
      include 'coupl.inc'
      include 'opt_info.inc'
      include 'icf_ch.inc'

***********************
* initialization
***********************
      call getarg(1,ii)
      call getarg(2,iseed)
      read(ii,*) ich
      read(iseed,*) seed

c      seed = 12345
      ich_tmp = ich ! the job number
      call bsinit(seed)
      ifprog = 1

      if (icf.eq.1) then
         if (ich_scheme.eq.0) then
            if (ich_type.eq.0) ich = ich_tmp
            if (ich_type.eq.1) ich = icf_perch1(ich_tmp) ! job number > ch number
            if (ich_type.eq.2) ich = icf_perch2(ich_tmp) 
            if (ich_type.eq.3) ich = icf_perch3(ich_tmp) 
            if (ich_type.eq.4) ich = icf_perch4(ich_tmp) 
         elseif (ich_scheme.eq.1) then
            if (ich_type.eq.0) ich = opt_per(ich_tmp) ! job number > opt-ch number > ch number
            if (ich_type.eq.1) ich = opt_per(icf_optch1(ich_tmp))            
            if (ich_type.eq.2) ich = opt_per(icf_optch2(ich_tmp))
            if (ich_type.eq.3) ich = opt_per(icf_optch3(ich_tmp))
            if (ich_type.eq.4) ich = opt_per(icf_optch4(ich_tmp))
         endif
      elseif (icf.eq.2) then
         if (ich_scheme.ne.1) then ! force to Optimized channel scheme
            write(*,*) "ERROR: ich_scheme should be 1 for icf = 2."
            write(*,*) "Stopping the program..."
            stop
         endif
         if (ich_type.eq.0) ich = opt_per(ich_tmp) ! job number > opt-ch number > ch number
         if (ich_type.eq.1) ich = opt_per(icf_optch1(ich_tmp))            
         if (ich_type.eq.2) ich = opt_per(icf_optch2(ich_tmp))
         if (ich_type.eq.3) ich = opt_per(icf_optch3(ich_tmp))
         if (ich_type.eq.4) ich = opt_per(icf_optch4(ich_tmp))
      endif
      call hmffrd
      call bsffrd
      call userin
****************************
* prepare output files
****************************
      if (ich_tmp.lt.10) then
         filename1 = "xsec_"//ii(1:1)//".dat"
         open(luxsec,file=filename1,status="replace")
         open(luerr,file="bases_err"//ii(1:1)//".txt",status="replace")
         open(lugrid,file="grid"//ii(1:1)//".dat",status="replace")
c         open(lumaxwgts,file="maxwgts"//ii(1:1)//".dat",status="replace")
c        open(lumaxwgts2,file="maxwgts_ref"//ii(1:1)//".dat",status="old")
      elseif (ich_tmp.lt.100) then
         filename2 = "xsec_"//ii(1:2)//".dat"
         open(luxsec,file=filename2,status="replace")
         open(luerr,file="bases_err"//ii(1:2)//".txt",status="replace")
         open(lugrid,file="grid"//ii(1:2)//".dat",status="replace")
c         open(lumaxwgts,file="maxwgts"//ii(1:2)//".dat",status="replace")
c        open(lumaxwgts2,file="maxwgts_ref"//ii(1:2)//".dat",status="old")
      elseif (ich_tmp.lt.1000) then
         filename3 = "xsec_"//ii(1:3)//".dat"
         open(luxsec,file=filename3,status="replace")
         open(luerr,file="bases_err"//ii(1:3)//".txt",status="replace")
         open(lugrid,file="grid"//ii(1:3)//".dat",status="replace")
c         open(lumaxwgts,file="maxwgts"//ii(1:3)//".dat",status="replace")
c        open(lumaxwgts2,file="maxwgts_ref"//ii(1:3)//".dat",status="old")
      endif

*************************
* integration
*************************
      npoints = 0
      nptcut_fails = 0
      nrcut_fails = 0
      call bases(bfunc, estim,sigma,ctime,it1,it2)
c      write(*,*) npoints,nptcut_fails,nrcut_fails
      ich = ich_tmp
      call bsinfo(lustdo)
c      call bhplot(lustdo)
***************
* save grids
***************
      call bswrit(lubsf,lugrid)

**********************************
* write bases information
**********************************
      call save_histo(ii)
c      call usrout

      write(luxsec,*) estim,sigma
      close(luxsec)
      close(luerr)
      close(lugrid)
c      close(lumaxwgts)
c      close(lumaxwgts2)

      end
