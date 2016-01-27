      subroutine make_cparam(nini,nfin,maxnch)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2013
C     ****************************************************
      implicit none
      include 'opt_info.inc'
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nini,nfin,maxnch,ndist
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      open(2,file='dist_param.inc',status='old')
      read(2,*) ndist
      close(2)

      open(1,file='cparam.inc',status='replace')
      write(1,*) "     complex*16 ci"
      write(1,*) "     parameter (ci = (0d0,1d0))"
      write(1,*) ""
      write(1,*) "     integer nini,nfin,next,ndist,nch,nwgt,maxnch"
      write(1,*) "     integer ipdf,ich_type,icf,ichwgt,isingle_ch,ich_scheme"
      write(1,*) "     parameter (nini =",nini,", nfin =", nfin,
     &     ", next=nini+nfin)"
      write(1,*) "     parameter (maxnch =",maxnch,")"
      write(1,'(a1,a24,a47,a16)') "c","     icf: 0:Exact CF Sum",
     &     ", 1:Simple CF Sampling, 2:Optimized CF Sampling",
     &     ", 3:Basic CF Sum"
      write(1,'(a1,a45,a24)') "c",
     &     "     ich_scheme: 0:Semi-peripheral channeling",
     &     ", 1:Optimized channeling"
      write(1,'(a1,a29,a34)') "c","     ich_type: 0:allch, 1:tch",
     &     ", 2:tch+sch, 3:tch+sch+2sch, 4:sch"
      write(1,'(a40,a21,i3,a8,i3,a17)') 
     &     "      parameter (icf = 2, ich_scheme = 1",
     &     ", ich_type = 0, nch =",nch_opt,
     &     ", nwgt =",nch_opt,", isingle_ch = 0)"
      write(1,'(a1,a25,a16)') "c","     ichwgt: 0:Diagram wgt",
     &     ", 1:Jacobian wgt"
      write(1,*) "     parameter (ichwgt = 1)"
      write(1,*) "     parameter (ipdf = 1)"
      write(1,*) "     parameter (ndist =",ndist,")"
      write(1,*) ""
      write(1,*) "     integer ndim"
      write(1,*) "     integer nwild"
      write(1,*) "     common /ndim/ ndim,nwild"
      write(1,*) ""
      write(1,*) "     integer ich,nsample_cf,ich_tmp,ichcf,cffact"
      write(1,*) "     common /ich/ ich,ich_tmp,ichcf,cffact"
      write(1,*) "     common /sample/ nsample_cf"
      write(1,*) ""
      write(1,*) "     integer chtype_max"
      write(1,*) "     parameter (chtype_max = 20)"
      write(1,*) "     real*8 AA(maxnch),AAA,CH(maxnch)",
     &     ",jfact(0:chtype_max)"
      write(1,*) "     common /AS/ AA,AAA,CH,jfact"
      write(1,*) ""
      write(1,*) "     integer ncf,floworder(next,2),flowtype"
      write(1,*) "     common /colorflow/ floworder,flowtype,ncf"
      write(1,*) "     character*5 ii"
      write(1,*) "     character*5 jj"
      write(1,*) "     character*5 kk"
      write(1,*) ""
      write(1,*) "     real*8 mgluon"
      write(1,*) "     common /mgluon/ mgluon" 
      write(1,*) ""
      write(1,*) "     real*8 mass(next)"
      write(1,*) "     common /masses/ mass"
      write(1,*) ""
      write(1,*) "     integer posflag,delyflag"
      write(1,*) "     common /pshk_flag/ posflag,delyflag"
      write(1,*) "     integer luxsec,luerr,lustdi,lustdo,lubsf"
     &     ,",luhist,luevtf,lugrid"
      write(1,*) "     parameter (lustdi=5,lustdo=6,lubsf=23"
     &     ,"    &     ,luhist=25,luevtf=26,luerr=99,luxsec=98"
     &     ,",lugrid=27)"
      write(1,*) ""
      write(1,*) "     integer nospev,nevprt"
      write(1,*) "     common /bsfile/ nospev,nevprt"
      write(1,*) "     character*128 bsoutf,bshstf,spevtf"
      write(1,*) "     common /bsfil0/ bsoutf,bshstf,spevtf"
      write(1,*) "     save /bsfile/,/bsfil0/"
      write(1,*) ""
      write(1,*) "     integer ifprog"
      write(1,*) "     common /hmflag/ ifprog"
      write(1,*) "     save /hmflag/"
      write(1,*) 
      write(1,*) "     integer npoints,nptcut_fails,nrcut_fails"
      write(1,*) "     common /npoints/ npoints,nptcut_fails"
     &     ,",nrcut_fails"
      close(1)

      return
      end
