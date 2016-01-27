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
nfin=`expr $ngluons - 2`
nfinm1=`expr $nfin - 1`
ndim=`expr 3 \* $nfin - 2`

out_file=set_histo.f
if [ -e $out_file ];then
    rm -rf $out_file
fi
touch $out_file

id_init=30
dist_count=0

pt_flag=1
y_flag=1
phi_flag=0
dr_flag=1
if [ $ngluons -ge 8 ];then
    m_flag=0
else
    m_flag=1
fi
kt_flag=0

nbin=50
pt_min="0d0"
pt_max="100d0"
y_min="0d0"
y_max="10d0"
phi_min="0d0"
phi_max="6.5d0"
dr_min="0d0"
dr_max="5d0"
m_min="0d0"
m_max="100d0"
kt_min="0d0"
kt_max="100d0"

echo "      subroutine set_histo(xl,xu)" >> ${out_file}
echo "      implicit none" >> ${out_file}
echo "C     GLOBAL VARIABLES" >> ${out_file}
echo "      include 'cparam.inc'" >> ${out_file}
echo "C     CONSTANTS" >> ${out_file}
echo "C     ARGUMENTS " >> ${out_file}
echo "C     LOCAL VARIABLES " >> ${out_file}
echo "      integer i" >> ${out_file}
echo "      character*80 ctit" >> ${out_file}
echo "      integer maxdim" >> ${out_file}
echo "      parameter (maxdim=50)" >> ${out_file}
echo "      real*8 xl(maxdim),xu(maxdim)" >> ${out_file}
echo "C     EXTERNAL FUNCTIONS" >> ${out_file}
echo "C     ----------" >> ${out_file}
echo "C     BEGIN CODE" >> ${out_file}
echo "C     ----------" >> ${out_file}
id=$id_init
id=`expr $id + 1`
echo "      call xhinit("$id",0d0,1d0,50,'dsigma/dx1')" >> ${out_file}
id=`expr $id + 1`
echo "      call xhinit("$id",0d0,1d0,50,'dsigma/dx2')" >> ${out_file}
id=`expr $id + 1`
echo "      call xhinit("$id",0d0,2d0,50,'dsigma/M2/SumM2')" >> ${out_file}
echo "c      call dhinit(22,0d0,1d0,30,0d0,1d0,30," >> ${out_file}
echo "c     &     'd sigma/ d x1/d x2')" >> ${out_file}
echo "      do i=1,ndim" >> ${out_file}
echo "         write(ctit,1010) i,xl(i),xu(i)" >> ${out_file}
echo " 1010    format('x(',i2,') ; ',e12.6,' - ',e12.6)" >> ${out_file}
echo "         call xhinit(i+"$id",xl(i),xu(i),50,ctit)" >> ${out_file}
echo "      enddo" >> ${out_file}
id=`expr $id + $ndim`
if [ $pt_flag -eq 1 ];then
    i=1
    while [ $i -le $nfin ];do
	id=`expr $id + 1`
	dist_count=`expr ${dist_count} + 1`
	echo "      call xhinit("$id","${pt_min}","${pt_max}","${nbin}",'dsigma/pt_g"$i"')" >> ${out_file}
	i=`expr $i + 1`
    done
fi
if [ $y_flag -eq 1 ];then
    i=1
    while [ $i -le $nfin ];do
	id=`expr $id + 1`
	dist_count=`expr ${dist_count} + 1`
	echo "      call xhinit("$id","${y_min}","${y_max}","${nbin}",'dsigma/y_g"$i"')" >> ${out_file}
	i=`expr $i + 1`
    done
fi
if [ $phi_flag -eq 1 ];then
    i=1
    while [ $i -le $nfin ];do
	id=`expr $id + 1`
	dist_count=`expr ${dist_count} + 1`
	echo "      call xhinit("$id","${phi_min}","${phi_max}","${nbin}",'dsigma/phi_g"$i"')" >> ${out_file}
	i=`expr $i + 1`
    done
fi
if [ $dr_flag -eq 1 ];then
    i=1
    while [ $i -le $nfinm1 ];do
	j=`expr $i + 1`
	while [ $j -le $nfin ];do
	    id=`expr $id + 1`
	    dist_count=`expr ${dist_count} + 1`
	    echo "      call xhinit("$id","${dr_min}","${dr_max}","${nbin}",'dsigma/dr_g"$i"g"$j"')" >> ${out_file}
	    j=`expr $j + 1`
	done
	i=`expr $i + 1`
    done
fi
if [ $m_flag -eq 1 ];then
    i=1
    while [ $i -le $nfinm1 ];do
	j=`expr $i + 1`
	while [ $j -le $nfin ];do
	    id=`expr $id + 1`
	    dist_count=`expr ${dist_count} + 1`
	    echo "      call xhinit("$id","${m_min}","${m_max}","${nbin}",'dsigma/m_g"$i"g"$j"')" >> ${out_file}
	    j=`expr $j + 1`
	done
	i=`expr $i + 1`
    done
fi
if [ $kt_flag -eq 1 ];then
    i=1
    while [ $i -le $nfinm1 ];do
	j=`expr $i + 1`
	while [ $j -le $nfin ];do
	    id=`expr $id + 1`
	    dist_count=`expr ${dist_count} + 1`
	    echo "      call xhinit("$id","${kt_min}","${kt_max}","${nbin}",'dsigma/kt_g"$i"g"$j"')" >> ${out_file}
	    j=`expr $j + 1`
	done
	i=`expr $i + 1`
    done
fi

nplots=`expr 3 + $ndim + $dist_count`
if [ $nplots -gt 50 ];then
    echo "ERROR: The number of plots exceeds 50. Please reduce plots."
    rm -rf $out_file
    exit
fi

echo "" >> ${out_file}
echo "      return" >> ${out_file}
echo "      end" >> ${out_file}

echo ${dist_count} ${id_init} ${ndim} > dist_param.inc