      complex*16 ci
      parameter (ci = (0d0,1d0))
 
      integer nini,nfin,next,ndist,nch,nwgt,maxnch
      integer ipdf,ich_type,icf,ichwgt,isingle_ch,ich_scheme
      parameter (nini = 2, nfin = 3, next=nini+nfin)
      parameter (maxnch = 200)
c     icf: 0:Exact CF Sum, 1:Simple CF Sampling, 2:Optimized CF Sampling, 3:Basic CF Sum
c     ich_scheme: 0:Semi-peripheral channeling, 1:Optimized channeling
c     ich_type: 0:allch, 1:tch, 2:tch+sch, 3:tch+sch+2sch, 4:sch
      parameter (icf = 2, ich_scheme = 1, ich_type = 0, nch =  8, nwgt =  8, isingle_ch = 0)
c     ichwgt: 0:Diagram wg, 1:Jacobian wgt
      parameter (ichwgt = 1)
      parameter (ipdf = 1)
      parameter (ndist = 12)
 
      integer ndim
      integer nwild
      common /ndim/ ndim,nwild
 
      integer ich,nsample_cf,ich_tmp,ichcf,cffact
      common /ich/ ich,ich_tmp,ichcf,cffact
      common /sample/ nsample_cf
 
      integer chtype_max
      parameter (chtype_max = 20)
      real*8 AA(maxnch),AAA,CH(maxnch),jfact(0:chtype_max)
      common /AS/ AA,AAA,CH,jfact
 
      integer ncf,floworder(next,2),flowtype
      common /colorflow/ floworder,flowtype,ncf
      character*5 ii
      character*5 jj
      character*5 kk
 
      real*8 mgluon
      common /mgluon/ mgluon
 
      real*8 mass(next)
      common /masses/ mass
 
      integer posflag,delyflag
      common /pshk_flag/ posflag,delyflag
      integer luxsec,luerr,lustdi,lustdo,lubsf,luhist,luevtf,lugrid
      parameter (lustdi=5,lustdo=6,lubsf=23
     &     ,luhist=25,luevtf=26,luerr=99,luxsec=98,lugrid=27)
 
      integer nospev,nevprt
      common /bsfile/ nospev,nevprt
      character*128 bsoutf,bshstf,spevtf
      common /bsfil0/ bsoutf,bshstf,spevtf
      save /bsfile/,/bsfil0/
 
      integer ifprog
      common /hmflag/ ifprog
      save /hmflag/

      integer npoints,nptcut_fails,nrcut_fails
      common /npoints/ npoints,nptcut_fails,nrcut_fails
