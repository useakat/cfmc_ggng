      integer ngluons
      parameter (ngluons = 5)
      integer max_natm,nnfin
      parameter (max_natm = ngluons-2, nnfin = ngluons-2)      
      integer nbcf
      parameter (nbcf = ngluons-1)
      integer maxnch
      parameter (maxnch = 1000)
      integer nch_per,ntch_per,nsch_per,n2sch_per,n3sch_per,n4sch_per
      integer nch_opt,ntch_opt,nsch_opt,n2sch_opt,n3sch_opt,n4sch_opt
      integer opt_ncf(maxnch),cf_nopt(nbcf)
      integer opt_per(maxnch),per_opt(maxnch)
      integer opt_cf(nbcf,maxnch),cf_opt(maxnch,nbcf)
      integer ipos_tch(maxnch)
      integer chmom(nnfin,maxnch)
      integer perch(max_natm,2,maxnch)
      integer perch_type(maxnch),optch_type(maxnch)		
      common /chinfo1/ nch_per,ntch_per
      common /chinfo2/ nsch_per,n2sch_per,n3sch_per,n4sch_per
      common /chinfo3/ nch_opt,ntch_opt
      common /chinfo4/ nsch_opt,n2sch_opt,n3sch_opt,n4sch_opt
      common /chinfo5/ opt_ncf,cf_nopt,opt_per,per_opt
      common /chinfo6/ opt_cf,cf_opt,ipos_tch,chmom
      common /chinfo7/ perch,perch_type,optch_type


