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
if [ $ngluons -gt 9 ]; then
    echo "Sorry, nlugons should be < 10."
    exit
fi 

out_file="find_perch.f"

if [ -e $out_file ]; then
    rm -rf $out_file
fi
touch $out_file

maxi=`expr 2 \* $ngluons - 6`
max_natm=`expr $ngluons - 2`

echo "      program find_perch" >> ${out_file}
echo "      implicitnone" >> ${out_file}
echo "C     CONSTANTS" >> ${out_file}
echo "C     GLOBAL VARIABLES" >> ${out_file}
echo "      include 'chinfo.inc'" >> ${out_file}
echo "C     ARGUMENTS" >> ${out_file}
echo "C     LOCAL VARIABLES " >> ${out_file}
echo "      integer i,j" >> ${out_file}
echo "      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14" >> ${out_file}
echo "      integer p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14" >> ${out_file}
echo "      integer natm,used(ngluons+1,$maxi),nchannels,cfflag" >> ${out_file}
echo "      integer atm(2,max_natm),cf(ngluons,ngluons-1),ierr,zero_flag" >> ${out_file}
echo "      integer nonzero_channels" >> ${out_file}
echo "C     EXTERNAL FUNCTIONS" >> ${out_file}
echo "C     ----------" >> ${out_file}
echo "C     BEGIN CODE" >> ${out_file}
echo "C     ----------" >> ${out_file}
echo "      nchannels = 0" >> ${out_file}
echo "      nonzero_channels = 0" >> ${out_file}
echo "      nch_per = 0" >> ${out_file}
echo "      nch_opt = 0" >> ${out_file}
echo "      ntch_per = 0" >> ${out_file}
echo "      ntch_opt = 0" >> ${out_file}
echo "      nsch_per = 0" >> ${out_file}
echo "      nsch_opt = 0" >> ${out_file}
echo "      n2sch_per = 0" >> ${out_file}
echo "      n2sch_opt = 0" >> ${out_file}
echo "      n3sch_per = 0" >> ${out_file}
echo "      n3sch_opt = 0" >> ${out_file}
echo "      n4sch_per = 0" >> ${out_file}
echo "      n4sch_opt = 0" >> ${out_file}
echo "      do i = 1,maxnch" >> ${out_file}
echo "         opt_ncf(i) = 0" >> ${out_file}
echo "      enddo" >> ${out_file}
echo "      do i = 1,nbcf" >> ${out_file}
echo "         cf_nopt(i) = 0" >> ${out_file}
echo "      enddo" >> ${out_file}
echo "" >> ${out_file}
echo "      do i =1,ngluons" >> ${out_file}
echo "         call gen_basiccf(i,ngluons,cf(1,i))" >> ${out_file}
echo "      enddo" >> ${out_file}
echo "" >> ${out_file}
i=1
while [ $i -le $maxi ]; do
    echo "      do i$i = 3,ngluons+1" >> ${out_file}
    echo "         do i = 3,ngluons+1" >> ${out_file}
    if [ $i -eq 1 ]; then
	echo "            used(i,1) = 0" >> ${out_file}
    else
	echo "            used(i,$i) = used(i,`expr $i - 1`)" >> ${out_file}
    fi
    echo "         enddo" >> ${out_file}
    echo "         if (i$i.eq.ngluons+1) then" >> ${out_file}
    echo "            p$i = 0" >> ${out_file}
    echo "            used(i$i,$i) = used(i$i,$i) +1" >> ${out_file}
    echo "            if (used(i$i,$i).gt.(ngluons-4)) cycle" >> ${out_file}
    echo "         else" >> ${out_file}
    echo "            if (used(i$i,$i).ne.0) then" >> ${out_file}
    echo "               cycle" >> ${out_file}
    echo "            else" >> ${out_file}
    echo "               p$i = i$i" >> ${out_file}
    echo "               used(i$i,$i) = 1" >> ${out_file}
    echo "            endif" >> ${out_file}
    echo "         endif" >> ${out_file}

    i=`expr $i + 1`
done
echo "         do i = 1,max_natm" >> ${out_file}
echo "            do j = 1,2" >> ${out_file}
echo "               atm(j,i) = 0" >> ${out_file}
echo "            enddo" >> ${out_file}
echo "         enddo" >> ${out_file}
echo "         atm(1,1) = 1" >> ${out_file}
echo "         atm(2,1) = p1" >> ${out_file}
echo "         atm(1,2) = 2" >> ${out_file}
echo "         atm(2,2) = p2" >> ${out_file}
i=1
while [ $i -le `expr $max_natm - 2` ];do
    j1=`expr 2 \* $i + 1`
    j2=`expr ${j1} + 1`
    iatm=`expr $i + 2`
#    echo "         natm = natm +1" >> ${out_file}
    echo "         atm(1,"$iatm") = p${j1}" >> ${out_file}
    echo "         atm(2,"$iatm") = p${j2}" >> ${out_file}
    i=`expr $i + 1`
done
echo "         call analyse_atm(ngluons,cf,max_natm,atm" >> ${out_file}
echo "     &        ,nchannels,nonzero_channels,ierr)" >> ${out_file}
echo "         if (ierr.eq.1) cycle" >> ${out_file}
i=1
while [ $i -le $maxi ];do 
    echo "      enddo" >> ${out_file}
    i=`expr $i + 1`
done
echo "      write(6,*) " >> ${out_file}
echo '      write(6,*) "nchannels =", nchannels' >> ${out_file}
echo '      write(6,*) "non-zero channels =", nonzero_channels' >> ${out_file}
echo "" >> ${out_file}
echo "      call write_chinfo" >> ${out_file}
echo "" >> ${out_file}
echo "      end" >> ${out_file}