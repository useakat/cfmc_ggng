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

out_file=fill_histo.f
if [ -e $out_file ];then
    rm -rf $out_file
fi
touch $out_file

read ndist id_init ndim < dist_param.inc

echo "      subroutine fill_histo(z,Q,amp2,bfunc)" >> ${out_file}
echo "      implicit none" >> ${out_file}
echo "C     GLOBAL VARIABLES" >> ${out_file}
echo "      include 'cparam.inc'" >> ${out_file}
echo "      include 'pshk.inc'" >> ${out_file}
echo "C     CONSTANTS" >> ${out_file}
echo "C     ARGUMENTS " >> ${out_file}
echo "      real*8 z(30),Q(0:3,next),amp2,bfunc" >> ${out_file}
echo "C     LOCAL VARIABLES" >> ${out_file}
echo "      integer i,j,k" >> ${out_file}
echo "      real*8 dist("${ndist}"),OQ(0:3,next),m2summ2 " >> ${out_file}
echo "C     EXTERNAL FUNCTIONS" >> ${out_file}
echo "      real*8 pdeltar,ppt,py,peta,pphi,pmij,pktij" >> ${out_file}
echo "      external pdeltar,ppt,py,peta,pphi,pmij,pktij" >> ${out_file}
echo "C     ----------" >> ${out_file}
echo "C     BEGIN CODE" >> ${out_file}
echo "C     ---------" >> ${out_file}
echo "      m2summ2 = amp2/AAA" >> ${out_file}
echo "      call order_pt(nfin,Q(0,3),OQ(0,3))" >> ${out_file}
echo "      i = 1" >> ${out_file}
echo "      do k = 3,next" >> ${out_file}
echo "         dist(i) = ppt(OQ(0,k))" >> ${out_file}
echo "         i = i +1" >> ${out_file}
echo "      enddo" >> ${out_file}
echo "      do k = 3,next" >> ${out_file}
echo "         dist(i) = 5 +py(OQ(0,k))" >> ${out_file}
echo "         i = i +1" >> ${out_file}
echo "      enddo" >> ${out_file}
if [ $ngluons -le -1 ]; then
    echo "      do k = 3,next" >> ${out_file}
    echo "         dist(i) = pphi(OQ(0,k))" >> ${out_file}
    echo "         i = i +1" >> ${out_file}
    echo "      enddo" >> ${out_file}
fi
echo "      do k = 3,next-1" >> ${out_file}
echo "         do j = k+1,next" >> ${out_file}
echo "            dist(i) = pdeltar(OQ(0,k),OQ(0,j))" >> ${out_file}
echo "            i = i +1" >> ${out_file}
echo "         enddo" >> ${out_file}
echo "      enddo" >> ${out_file}
if [ $ngluons -lt 8 ];then
    echo "      do k = 3,next-1" >> ${out_file}
    echo "         do j = k+1,next" >> ${out_file}
    echo "            dist(i) = pmij(OQ(0,k),OQ(0,j))" >> ${out_file}
    echo "            i = i +1" >> ${out_file}
    echo "         enddo" >> ${out_file}
    echo "      enddo" >> ${out_file}
fi
if [ $ngluons -le -1 ]; then
    echo "      do k = 3,next-1" >> ${out_file}
    echo "         do j = k+1,next" >> ${out_file}
    echo "            dist(i) = pktij(OQ(0,k),OQ(0,j))" >> ${out_file}
    echo "            i = i +1" >> ${out_file}
    echo "         enddo" >> ${out_file}
    echo "      enddo" >> ${out_file}
fi
echo "" >> ${out_file}
id=${id_init}
id=`expr $id + 1`
echo "      call xhfill("$id",x1,bfunc)" >> ${out_file}
id=`expr $id + 1`
echo "      call xhfill("$id",x2,bfunc)" >> ${out_file}
id=`expr $id + 1`
echo "      call xhfill("$id",m2summ2,bfunc)" >> ${out_file}
echo "c      call dhfill("$id",x1,x2,bfunc)" >> ${out_file}
echo "      do i=1,ndim" >> ${out_file}
echo "         call xhfill("$id"+i,z(i),bfunc)" >> ${out_file}
echo "      enddo" >> ${out_file}
id=`expr $id + $ndim`
echo "      do i = 1,ndist" >> ${out_file}
echo "         call xhfill("$id"+i,dist(i),bfunc)" >> ${out_file}
echo "      enddo" >> ${out_file}
echo "" >> ${out_file}
echo "      return" >> ${out_file}
echo "      end" >> ${out_file}



    