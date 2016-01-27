#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: fortran_temp.sh [type]"
    echo ""
    echo " [type] template type: pro or sub"
    echo ""
    exit
fi

ngluons=$1

out_file="chinfo.inc"
if [ -e $out_file ];then
    rm -rf $out_file
fi
touch $out_file

echo "      integer ngluons" >> ${out_file}
echo "      parameter (ngluons = "$ngluons")" >> ${out_file}
echo "      integer max_natm,nnfin" >> ${out_file}
echo "      parameter (max_natm = ngluons-2, nnfin = ngluons-2)      " >> ${out_file}
echo "      integer nbcf" >> ${out_file}
echo "      parameter (nbcf = ngluons-1)" >> ${out_file}
echo "      integer maxnch" >> ${out_file}
echo "      parameter (maxnch = 1000)" >> ${out_file}
echo "      integer nch_per,ntch_per,nsch_per,n2sch_per,n3sch_per,n4sch_per" >> ${out_file}
echo "      integer nch_opt,ntch_opt,nsch_opt,n2sch_opt,n3sch_opt,n4sch_opt" >> ${out_file}
echo "      integer opt_ncf(maxnch),cf_nopt(nbcf)" >> ${out_file}
echo "      integer opt_per(maxnch),per_opt(maxnch)" >> ${out_file}
echo "      integer opt_cf(nbcf,maxnch),cf_opt(maxnch,nbcf)" >> ${out_file}
echo "      integer ipos_tch(maxnch)" >> ${out_file}
echo "      integer chmom(nnfin,maxnch)" >> ${out_file}
echo "      integer perch(max_natm,2,maxnch)" >> ${out_file}
echo "      integer perch_type(maxnch),optch_type(maxnch)		" >> ${out_file}
echo "      common /chinfo1/ nch_per,ntch_per" >> ${out_file}
echo "      common /chinfo2/ nsch_per,n2sch_per,n3sch_per,n4sch_per" >> ${out_file}
echo "      common /chinfo3/ nch_opt,ntch_opt" >> ${out_file}
echo "      common /chinfo4/ nsch_opt,n2sch_opt,n3sch_opt,n4sch_opt" >> ${out_file}
echo "      common /chinfo5/ opt_ncf,cf_nopt,opt_per,per_opt" >> ${out_file}
echo "      common /chinfo6/ opt_cf,cf_opt,ipos_tch,chmom" >> ${out_file}
echo "      common /chinfo7/ perch,perch_type,optch_type" >> ${out_file}
echo "" >> ${out_file}
echo "" >> ${out_file}



    