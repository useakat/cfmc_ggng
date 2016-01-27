#!/bin/bash

source plot_info.sh

nplot=12
nplotp1=`expr ${nplot} + 1`

cd Events/$1
rm -r *.dat
grep '0.0000E+00' ma_plot003.top | cut -c5-100 > ${plot_unwgt[0]}_tmp.dat
grep '0.0000E+00' ma_plot004.top | cut -c5-100 > ${plot_unwgt[1]}_tmp.dat
grep '0.0000E+00' ma_plot005.top | cut -c5-100 > ${plot_unwgt[2]}_tmp.dat
grep '0.0000E+00' ma_plot006.top | cut -c5-100 > ${plot_unwgt[3]}_tmp.dat
grep '0.0000E+00' ma_plot007.top | cut -c5-100 > ${plot_unwgt[4]}_tmp.dat
grep '0.0000E+00' ma_plot008.top | cut -c5-100 > ${plot_unwgt[5]}_tmp.dat
grep '0.0000E+00' ma_plot009.top | cut -c5-100 > ${plot_unwgt[6]}_tmp.dat
grep '0.0000E+00' ma_plot010.top | cut -c5-100 > ${plot_unwgt[7]}_tmp.dat
grep '0.0000E+00' ma_plot011.top | cut -c5-100 > ${plot_unwgt[8]}_tmp.dat
grep '0.0000E+00' ma_plot012.top | cut -c5-100 > ${plot_unwgt[9]}_tmp.dat
grep '0.0000E+00' ma_plot013.top | cut -c5-100 > ${plot_unwgt[10]}_tmp.dat
grep '0.0000E+00' ma_plot014.top | cut -c5-100 > ${plot_unwgt[11]}_tmp.dat
grep '0.0000E+00' ma_plot015.top | cut -c5-100 > ${plot_unwgt[12]}_tmp.dat
grep '0.0000E+00' ma_plot016.top | cut -c5-100 > ${plot_unwgt[13]}_tmp.dat
grep '0.0000E+00' ma_plot017.top | cut -c5-100 > ${plot_unwgt[14]}_tmp.dat
grep '0.0000E+00' ma_plot018.top | cut -c5-100 > ${plot_unwgt[15]}_tmp.dat
grep '0.0000E+00' ma_plot019.top | cut -c5-100 > ${plot_unwgt[16]}_tmp.dat
grep '0.0000E+00' ma_plot020.top | cut -c5-100 > ${plot_unwgt[17]}_tmp.dat
grep '0.0000E+00' ma_plot021.top | cut -c5-100 > ${plot_unwgt[18]}_tmp.dat
grep '0.0000E+00' ma_plot022.top | cut -c5-100 > ${plot_unwgt[19]}_tmp.dat
grep '0.0000E+00' ma_plot023.top | cut -c5-100 > ${plot_unwgt[20]}_tmp.dat
grep '0.0000E+00' ma_plot024.top | cut -c5-100 > ${plot_unwgt[21]}_tmp.dat
grep '0.0000E+00' ma_plot025.top | cut -c5-100 > ${plot_unwgt[22]}_tmp.dat
grep '0.0000E+00' ma_plot026.top | cut -c5-100 > ${plot_unwgt[23]}_tmp.dat
grep '0.0000E+00' ma_plot027.top | cut -c5-100 > ${plot_unwgt[24]}_tmp.dat
grep '0.0000E+00' ma_plot028.top | cut -c5-100 > ${plot_unwgt[25]}_tmp.dat
grep '0.0000E+00' ma_plot029.top | cut -c5-100 > ${plot_unwgt[26]}_tmp.dat
grep '0.0000E+00' ma_plot030.top | cut -c5-100 > ${plot_unwgt[27]}_tmp.dat
grep '0.0000E+00' ma_plot031.top | cut -c5-100 > ${plot_unwgt[28]}_tmp.dat
grep '0.0000E+00' ma_plot032.top | cut -c5-100 > ${plot_unwgt[29]}_tmp.dat

read wgt < ../../uwgt.dat

i=0
while [ $i -lt ${nplot} ]; do
    touch ${plot_unwgt[$i]}.dat
    while read x y z; do
	yr=${y%%E*}
	ye=${y##*E}
	iye=${ye##*+0}
	iye=${iye##*-0}
	iye=${iye##*+}
	iye=${iye##*-}
	yy=`echo "${yr}*10^${iye}" | bc` 
	nevent=`echo "scale=3; ${yy} / ${wgt}" | bc`
	sqrtnevent=`echo "sqrt(${nevent})" | bc`
	dy=`echo "${sqrtnevent}*${wgt}" | bc`
	echo $x $y ${dy} >> ${plot_unwgt[$i]}.dat
    done < ${plot_unwgt[$i]}_tmp.dat
    i=`expr $i + 1`
done