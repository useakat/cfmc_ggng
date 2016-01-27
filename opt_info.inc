      integer oi,oj
      integer nnext
      parameter (nnext = 5)
      integer nbcf,nnfin
      parameter (nbcf = nnext-1, nnfin = nnext-2)
      integer nch_per,ntch_per,nch_opt,ntch_opt
      integer nsch_per,nsch_opt
      integer n2sch_per,n2sch_opt
      integer n3sch_per,n3sch_opt
      integer n4sch_per,n4sch_opt
      parameter (nch_per = 10, ntch_per = 6)
      parameter (nch_opt = 8, ntch_opt = 4)
      parameter (nsch_per = 4, nsch_opt = 4)
      integer opt_ncf(nch_opt),cf_nopt(nbcf)
      data opt_ncf /  2,  2,  2,  2,  2,  2,  2,  2 /
      data cf_nopt /  3,  5,  5,  3 /
 
      integer opt_per(nch_opt),per_opt(nch_per)
      integer opt_cf(nnext-1,nch_opt), cf_opt(nch_opt,nnext-1)
      data opt_per /  1,  2,  3,  6,  7,  8,  9, 10 /
      data per_opt /  1,  2,  3,  0,  0,  4,  5,  6,  7,  8 /
 
      data (opt_cf(oi,1),oi=1,2) /  2,  3 /
      data (opt_cf(oi,2),oi=1,2) /  3,  4 /
      data (opt_cf(oi,3),oi=1,2) /  2,  4 /
      data (opt_cf(oi,4),oi=1,2) /  1,  2 /
      data (opt_cf(oi,5),oi=1,2) /  2,  3 /
      data (opt_cf(oi,6),oi=1,2) /  1,  3 /
      data (opt_cf(oi,7),oi=1,2) /  1,  2 /
      data (opt_cf(oi,8),oi=1,2) /  3,  4 /
 
      data (cf_opt(oi,1),oi=1,3) /  4,  6,  7 /
      data (cf_opt(oi,2),oi=1,5) /  1,  3,  4,  5,  7 /
      data (cf_opt(oi,3),oi=1,5) /  1,  2,  5,  6,  8 /
      data (cf_opt(oi,4),oi=1,3) /  2,  3,  8 /
 
      integer ipos_tch(ntch_per)
      data ipos_tch /  0,  1,  0,  2,  1,  2 /
 
      integer chmom(nnfin,nch_per)
      data (chmom(oi,1),oi=1,nnfin) /  3,  4,  5 /
      data (chmom(oi,2),oi=1,nnfin) /  3,  5,  4 /
      data (chmom(oi,3),oi=1,nnfin) /  3,  4,  5 /
      data (chmom(oi,4),oi=1,nnfin) /  4,  3,  5 /
      data (chmom(oi,5),oi=1,nnfin) /  5,  3,  4 /
      data (chmom(oi,6),oi=1,nnfin) /  4,  5,  3 /
      data (chmom(oi,7),oi=1,nnfin) /  5,  4,  3 /
      data (chmom(oi,8),oi=1,nnfin) /  4,  5,  3 /
      data (chmom(oi,9),oi=1,nnfin) /  3,  4,  5 /
      data (chmom(oi,10),oi=1,nnfin) /  4,  5,  3 /
 
      integer max_natm
      parameter (max_natm = nnext-2)
      integer perch(max_natm,2,nch_per)
      data (perch(oi,1,1),oi=1,max_natm) /  1,  2,  5 /
      data (perch(oi,2,1),oi=1,max_natm) /  3,  4,  0 /
      data (perch(oi,1,2),oi=1,max_natm) /  1,  2,  4 /
      data (perch(oi,2,2),oi=1,max_natm) /  3,  5,  0 /
      data (perch(oi,1,3),oi=1,max_natm) /  1,  2,  5 /
      data (perch(oi,2,3),oi=1,max_natm) /  3,  0,  4 /
      data (perch(oi,1,4),oi=1,max_natm) /  1,  2,  5 /
      data (perch(oi,2,4),oi=1,max_natm) /  4,  3,  0 /
      data (perch(oi,1,5),oi=1,max_natm) /  1,  2,  3 /
      data (perch(oi,2,5),oi=1,max_natm) /  4,  5,  0 /
      data (perch(oi,1,6),oi=1,max_natm) /  1,  2,  4 /
      data (perch(oi,2,6),oi=1,max_natm) /  5,  3,  0 /
      data (perch(oi,1,7),oi=1,max_natm) /  1,  2,  3 /
      data (perch(oi,2,7),oi=1,max_natm) /  5,  4,  0 /
      data (perch(oi,1,8),oi=1,max_natm) /  1,  2,  4 /
      data (perch(oi,2,8),oi=1,max_natm) /  5,  0,  3 /
      data (perch(oi,1,9),oi=1,max_natm) /  1,  2,  5 /
      data (perch(oi,2,9),oi=1,max_natm) /  0,  3,  4 /
      data (perch(oi,1,10),oi=1,max_natm) /  1,  2,  4 /
      data (perch(oi,2,10),oi=1,max_natm) /  0,  5,  3 /
 
      integer perch_type(nch_per),optch_type(nch_opt)
      data perch_type /  0,  0,  1,  0,  0,  0,  0,  1,  1,  1 /
      data optch_type /  0,  0,  1,  0,  0,  1,  1,  1 /
