#!/bin/bash

# script file to make .dat file in MG/ME/Myproc/Events/ (KM: 2011/12/02)

run=$1
nch=$2
ndim=$3
ndist=$4
chstart=$5
chp1=$6

nchp1=`expr $nch + 1`
ndimp1=`expr $ndim + 1`
ndistp1=`expr $ndist + 1`

cd rslt_${run}/data

i=$chstart
while [ $i -lt $chp1 ]; do
    grep '0.0000E+00' x1_$i.top | cut -c3-100 > x1_$i.dat
    grep '0.0000E+00' x2_$i.top | cut -c3-100 > x2_$i.dat
    grep '0.0000E+00' m2summ2_$i.top | cut -c3-100 > m2summ2_$i.dat
    j=1
    while [ $j -lt $ndimp1 ]; do
	grep '0.0000E+00' z${j}_$i.top | cut -c3-100 > z${j}_$i.dat
	j=`expr $j + 1`
    done
    j=1
    while [ $j -lt $ndistp1 ]; do
	grep '0.0000E+00' dist${j}_$i.top | cut -c3-100 > dist${j}_$i.dat
	j=`expr $j + 1`
    done
    i=`expr $i + 1`
done