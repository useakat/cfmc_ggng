      subroutine write_chinfo
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 24 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'chinfo.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      character*3 cformat,cii,cnopt,cncf
      integer i,j,k,nrem,iformat
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      open(1,file="opt_info.inc",status='replace')

      write(1,*) "     integer oi,oj"
      write(1,*) "     integer nnext"
      write(1,*) "     parameter (nnext =",ngluons,")"
      write(1,*) "     integer nbcf,nnfin"
      write(1,*) "     parameter (nbcf = nnext-1, nnfin = nnext-2)"
      write(1,*) "     integer nch_per,ntch_per,nch_opt,ntch_opt"
      write(1,*) "     integer nsch_per,nsch_opt"
      write(1,*) "     integer n2sch_per,n2sch_opt"
      write(1,*) "     integer n3sch_per,n3sch_opt"
      write(1,*) "     integer n4sch_per,n4sch_opt"
      write(1,*) "     parameter (nch_per =",nch_per,", ntch_per ="
     &     ,ntch_per,")"
      write(1,*) "     parameter (nch_opt =",nch_opt,", ntch_opt ="
     &     ,ntch_opt,")"
      if (ngluons.ge.5) then
         write(1,*) "     parameter (nsch_per =",
     &        nsch_per,", nsch_opt =",nsch_opt,")"
      endif
      if (ngluons.ge.6) then
         write(1,*) "     parameter (n2sch_per =",
     &        n2sch_per,", n2sch_opt =",n2sch_opt,")"
      endif
      if (ngluons.ge.8) then
         write(1,*) "     parameter (n3sch_per =",
     &        n3sch_per,", n3sch_opt =",n3sch_opt,")"
      endif
      if (ngluons.ge.10) then
         write(1,*) "     parameter (n4sch_per =",
     &        n4sch_per,", n4sch_opt =",n4sch_opt,")"
      endif
      write(1,*) "     integer opt_ncf(nch_opt),cf_nopt(nbcf)"
      call data_array(opt_ncf,nch_opt,"opt_ncf",7,1)
      call data_array(cf_nopt,nbcf,"cf_nopt",7,1)
      write(1,*) ""
      write(1,*) "     integer opt_per(nch_opt),per_opt(nch_per)"
      write(1,*) "     integer opt_cf(nnext-1,nch_opt)"
     &     ,", cf_opt(nch_opt,nnext-1)"
      call data_array(opt_per,nch_opt,"opt_per",7,1)
      call data_array(per_opt,nch_per,"per_opt",7,1)
      write(1,*) ""
      if (nch_opt.lt.10) then
         do i = 1,nch_opt
            write(cii,'(i1)') i
            write(cncf,'(i1)') opt_ncf(i)
            call data_array(opt_cf(1,i),opt_ncf(i),
     &          "(opt_cf(oi,"//cii(1:1)//"),oi=1,"//cncf(1:1)//")",21,1)
         enddo
      elseif (nch_opt.lt.100) then
         do i = 1,9
            write(cii,'(i1)') i
            write(cncf,'(i1)') opt_ncf(i)
            call data_array(opt_cf(1,i),opt_ncf(i),
     &          "(opt_cf(oi,"//cii(1:1)//"),oi=1,"//cncf(1:1)//")",21,1)
         enddo
         do i = 10,nch_opt
            write(cii,'(i2)') i
            write(cncf,'(i1)') opt_ncf(i)
            call data_array(opt_cf(1,i),opt_ncf(i),
     &          "(opt_cf(oi,"//cii(1:2)//"),oi=1,"//cncf(1:1)//")",22,1)
         enddo
      elseif (nch_opt.lt.1000) then
         do i = 1,9
            write(cii,'(i1)') i
            write(cncf,'(i1)') opt_ncf(i)
            call data_array(opt_cf(1,i),opt_ncf(i),
     &          "(opt_cf(oi,"//cii(1:1)//"),oi=1,"//cncf(1:1)//")",21,1)
         enddo
         do i = 10,99
            write(cii,'(i2)') i
            write(cncf,'(i1)') opt_ncf(i)
            call data_array(opt_cf(1,i),opt_ncf(i),
     &          "(opt_cf(oi,"//cii(1:2)//"),oi=1,"//cncf(1:1)//")",22,1)
         enddo
         do i = 100,999
            write(cii,'(i3)') i
            write(cncf,'(i1)') opt_ncf(i)
            call data_array(opt_cf(1,i),opt_ncf(i),
     &          "(opt_cf(oi,"//cii(1:3)//"),oi=1,"//cncf(1:1)//")",23,1)
         enddo
      else
         write(99,*) "Error: nch_opt >= 1000 is not supported."
      endif

      write(1,*) ""
      do i = 1,nbcf
         write(cii,'(i1)') i
         if (cf_nopt(i).lt.10) then
            write(cnopt,'(i1)') cf_nopt(i)
            call data_array(cf_opt(1,i),cf_nopt(i),
     &           "(cf_opt(oi,"//cii(1:1)//"),oi=1,"//cnopt(1:1)//")"
     &           ,21,1)
         elseif (cf_nopt(i).lt.100) then
            write(cnopt,'(i2)') cf_nopt(i)
            call data_array(cf_opt(1,i),cf_nopt(i),
     &           "(cf_opt(oi,"//cii(1:1)//"),oi=1,"//cnopt(1:2)//")"
     &           ,22,1)
         elseif (cf_nopt(i).lt.1000) then
            write(cnopt,'(i3)') cf_nopt(i)
            call data_array(cf_opt(1,i),cf_nopt(i),
     &           "(cf_opt(oi,"//cii(1:1)//"1),oi=1,"//cnopt(1:3)//")"
     &           ,23,1)
         endif
      enddo
      write(1,*) ""
      write(1,*) "     integer ipos_tch(ntch_per)"
      call data_array(ipos_tch,ntch_per,"ipos_tch",8,1)

      write(1,*) ""
      write(1,*) "     integer chmom(nnfin,nch_per)"
      do i = 1,nch_per
         if (i.lt.10) then
            write(cii,'(i1)') i
            call data_array(chmom(1,i),nnfin,
     &           "(chmom(oi,"//cii(1:1)//"),oi=1,nnfin)",24,1)
         elseif (i.lt.100) then
            write(cii,'(i2)') i
            call data_array(chmom(1,i),nnfin,
     &           "(chmom(oi,"//cii(1:2)//"),oi=1,nnfin)",25,1)
         elseif (i.lt.1000) then
            write(cii,'(i3)') i
            call data_array(chmom(1,i),nnfin,
     &           "(chmom(oi,"//cii(1:3)//"),oi=1,nnfin)",26,1)
         endif
      enddo

      write(1,*) ""
      write(1,*) "     integer max_natm"
      write(1,*) "     parameter (max_natm = nnext-2)"
      write(1,*) "     integer perch(max_natm,2,nch_per)"
      do i = 1,nch_per
         if (i.lt.10) then
            write(cii,'(i1)') i
            call data_array(perch(1,1,i),max_natm,
     &           "(perch(oi,1,"//cii(1:1)//"),oi=1,max_natm)",29,1)
            call data_array(perch(1,2,i),max_natm,
     &           "(perch(oi,2,"//cii(1:1)//"),oi=1,max_natm)",29,1)
         elseif (i.lt.100) then
            write(cii,'(i2)') i
            call data_array(perch(1,1,i),max_natm,
     &           "(perch(oi,1,"//cii(1:2)//"),oi=1,max_natm)",30,1)
            call data_array(perch(1,2,i),max_natm,
     &           "(perch(oi,2,"//cii(1:2)//"),oi=1,max_natm)",30,1)
         elseif (i.lt.1000) then
            write(cii,'(i3)') i
            call data_array(perch(1,1,i),max_natm,
     &           "(perch(oi,1,"//cii(1:3)//"),oi=1,max_natm)",31,1)
            call data_array(perch(1,2,i),max_natm,
     &           "(perch(oi,2,"//cii(1:3)//"),oi=1,max_natm)",31,1)
         endif
      enddo

      write(1,*) ""
      write(1,*) "     integer perch_type(nch_per),optch_type(nch_opt)"
      call data_array(perch_type,nch_per,"perch_type",10,1)
      call data_array(optch_type,nch_opt,"optch_type",10,1)

      return
      end
