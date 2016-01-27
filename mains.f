
      program main
      implicit none

      include 'lhe.inc'
      include 'cparam.inc'
      include 'opt_info.inc'
      include 'icf_ch.inc'
      include 'inc/base4.h'
      integer i,seed
      double precision bfunc
      external bfunc
      character*4 iseed
      character*1 genmode
      integer iret,mxtry
      parameter (mxtry=100000)
c      parameter (mxtry=122)

***********************
* initialization
***********************
      call getarg(1,ii)
      call getarg(2,iseed)
      call getarg(3,genmode)
      read(genmode,*) IDWTUP
      read(ii,*) ich
      read(iseed,*) seed

      seed = 12345
      ich_tmp = ich

****************************
* prepare output files
****************************
      if (ich_tmp.lt.10) then
         open(luerr,file="spring_err"//ii(1:1)//".txt",status="replace")
      elseif (ich_tmp.lt.100) then
         open(luerr,file="spring_err"//ii(1:2)//".txt",status="replace")
      elseif (ich_tmp.lt.1000) then
         open(luerr,file="spring_err"//ii(1:3)//".txt",status="replace")
      endif
*****************************

      call bsinit(seed)
      ifprog = 2

      if (icf.eq.1) then
         if (ich_type.eq.0) ich = ich_tmp
         if (ich_type.eq.1) ich = icf_perch1(ich_tmp)
         if (ich_type.eq.2) ich = icf_perch2(ich_tmp)
         if (ich_type.eq.3) ich = icf_perch3(ich_tmp)
         if (ich_type.eq.4) ich = icf_optch4(ich_tmp)
      elseif (icf.eq.2) then
         if (ich_type.eq.0) ich = opt_per(ich_tmp)
         if (ich_type.eq.1) ich = opt_per(icf_optch1(ich_tmp))            
         if (ich_type.eq.2) ich = opt_per(icf_optch2(ich_tmp))
         if (ich_type.eq.3) ich = opt_per(icf_optch3(ich_tmp))
         if (ich_type.eq.4) ich = opt_per(icf_optch4(ich_tmp))
      endif
      call hmffrd
      call bsffrd
      call userin
      call bsread(lubsf)
      close(lubsf)

      call drnset(seed)
*******************************
* read process information
*******************************
      open(1,file="Cards/process.dat",status="unknown")
      read(1,*) IDBMUP(1),IDBMUP(2),NPRUP,IDPRUP
      do i = 1,next
         read(1,*) IDUP(i),ISTUP(i),MOTHUP(1,i),MOTHUP(2,i)
      enddo
      NUP = next
      XSECUP(1) = 1
      XERRUP(1) = 0.1d0
      XMAXUP(1) = 1

      write(luevtf,*) '<LesHouchesEvents version="1.0">'
      write(luevtf,*) '<!--'
      write(luevtf,*) '-->'
      write(luevtf,*) '<header>'
      write(luevtf,*) '</header>'
      write(luevtf,*) '<init>'
      write(luevtf,200) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     &     PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
      do i = 1,NPRUP
         write(luevtf,*) XSECUP(1),XERRUP(1),XMAXUP(1),NPRUP
      enddo
      write(luevtf,*) '</init>'
      
      do i = 1, nospev
         call spring(bfunc,mxtry)
         call spevnt(iret)
      enddo

      write(99,*) "After generation"
      do i = 1,10
         write(99,*) "HeperCube",I
         write(99,*) "DXP",DXP(I)
         write(99,*) "DXD",DXD(I)
      enddo

      ich = ich_tmp
      write(luevtf,*) '</LesHouchesEvents>'
      call spinfo(lustdo)
      call shplot(lustdo)

      close(luerr)

 200  format(2(i9),2(e19.11),2(i2),2(i6),i2,i3)

      stop
      end
