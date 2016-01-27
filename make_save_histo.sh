#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: fortran_temp.sh [type]"
    echo ""
    echo " [type] template type: pro or sub"
    echo ""
    exit
fi

out_file=save_histo.f
if [ -e $out_file ];then
    rm -rf $out_file
fi
touch $out_file

read ndist id_init ndim < dist_param.inc

echo "      subroutine save_histo(ii)" >> ${out_file}
echo "      implicit none" >> ${out_file}
echo "C     GLOBAL VARIABLES" >> ${out_file}
echo "      include 'cparam.inc'" >> ${out_file}
echo "C     CONSTANTS" >> ${out_file}
echo "C     ARGUMENTS " >> ${out_file}
echo "C     LOCAL VARIABLES " >> ${out_file}
echo "      integer i" >> ${out_file}
echo "      character*4 iii" >> ${out_file}
echo "C     EXTERNAL FUNCTIONS" >> ${out_file}
echo "C     ----------" >> ${out_file}
echo "C     BEGIN CODE" >> ${out_file}
echo "C     ----------" >> ${out_file}
echo "      if (ich.lt.10) then" >> ${out_file}
id=$id_init
id=`expr $id + 1`
echo "         open("$id",file='bases_data/x1_'//ii(1:1)//'.top'" >> ${out_file}
echo "     &        ,status='replace')" >> ${out_file}
id=`expr $id + 1`
echo "         open("$id",file='bases_data/x2_'//ii(1:1)//'.top'" >> ${out_file}
echo "     &        ,status='replace')" >> ${out_file}
id=`expr $id + 1`
echo "         open("$id",file='bases_data/m2summ2_'//ii(1:1)//'.top'" >> ${out_file}
echo "     &        ,status='replace')" >> ${out_file}
echo "" >> ${out_file}
echo "         do i = 1,ndim" >> ${out_file}
echo "            write(iii,*) i " >> ${out_file}
echo "            if (i.lt.10) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/z'//iii(2:2)//'_'" >> ${out_file}
echo "     &              //ii(1:1)//'.top',status='replace')" >> ${out_file}
echo "            elseif (i.lt.100) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/z'//iii(2:3)//'_'" >> ${out_file}
echo "     &              //ii(1:1)//'.top',status='replace')" >> ${out_file}
echo "            endif" >> ${out_file}
echo "         enddo" >> ${out_file}
echo "" >> ${out_file}
id=`expr $id + $ndim`
echo "         do i = 1,ndist" >> ${out_file}
echo "            write(iii,*) i" >> ${out_file}
echo "            if (i.lt.10) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/dist'//iii(2:2)//'_'" >> ${out_file}
echo "     &              //ii(1:1)//'.top',status='replace')" >> ${out_file}
echo "            elseif (i.lt.100) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/dist'//iii(2:3)//'_'" >> ${out_file}
echo "     &              //ii(1:1)//'.top',status='replace')" >> ${out_file}
echo "            endif" >> ${out_file}
echo "         enddo" >> ${out_file}
echo "" >> ${out_file}
echo "      elseif (ich.lt.100) then" >> ${out_file}
id=$id_init
id=`expr $id + 1`
echo "         open("$id",file='bases_data/x1_'//ii(1:2)//'.top'" >> ${out_file}
echo "     &        ,status='replace')" >> ${out_file}
id=`expr $id + 1`
echo "         open("$id",file='bases_data/x2_'//ii(1:2)//'.top'" >> ${out_file}
echo "     &        ,status='replace')" >> ${out_file}
id=`expr $id + 1`
id_init1=$id
echo "         open("$id",file='bases_data/m2summ2_'//ii(1:2)//'.top'" >> ${out_file}
echo "     &        ,status='replace')" >> ${out_file}
echo "" >> ${out_file}
echo "         do i = 1,ndim" >> ${out_file}
echo "            write(iii,*) i " >> ${out_file}
echo "            if (i.lt.10) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/z'//iii(2:2)//'_'" >> ${out_file}
echo "     &              //ii(1:2)//'.top',status='replace')" >> ${out_file}
echo "            elseif (i.lt.100) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/z'//iii(2:3)//'_'" >> ${out_file}
echo "     &              //ii(1:2)//'.top',status='replace')" >> ${out_file}
echo "            endif" >> ${out_file}
echo "         enddo" >> ${out_file}
echo "" >> ${out_file}
id=`expr $id + $ndim`
echo "         do i = 1,ndist" >> ${out_file}
echo "            write(iii,*) i" >> ${out_file}
echo "            if (i.lt.10) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/dist'//iii(2:2)//'_'" >> ${out_file}
echo "     &              //ii(1:2)//'.top',status='replace')" >> ${out_file}
echo "            elseif (i.lt.100) then" >> ${out_file}
echo "               open("$id"+i,file='bases_data/dist'//iii(2:3)//'_'" >> ${out_file}
echo "     &              //ii(1:2)//'.top',status='replace')" >> ${out_file}
echo "            endif" >> ${out_file}
echo "         enddo" >> ${out_file}
echo "      endif      " >> ${out_file}
echo "" >> ${out_file}
id=`expr $id_init + 1`
echo "      do i = "$id","$id_init1 >> ${out_file}
echo "         call xhsave(i,i)" >> ${out_file}
echo "         close(i)" >> ${out_file}
echo "      enddo" >> ${out_file}
id=`expr $id_init1 + 1`
echo "      do i = "$id","${id_init1}"+ndim" >> ${out_file}
echo "         call xhsave(i,i)         " >> ${out_file}
echo "         close(i)" >> ${out_file}
echo "      enddo" >> ${out_file}
id_init2=`expr $id_init1 + $ndim`
id=`expr $id_init2 + 1`
echo "      do i = "$id","${id_init2}"+ndist" >> ${out_file}
echo "         call xhsave(i,i)" >> ${out_file}
echo "         close(i)" >> ${out_file}
echo "      enddo" >> ${out_file}
echo "" >> ${out_file}
echo "      return" >> ${out_file}
echo "      end" >> ${out_file}
    