#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: mkgnu.sh [run1] [run2] ... [runn]"
    echo ""
    exit
fi
if [[ $1 = "" ]]; then
    echo "input run name"
    echo "exiting..."
    exit
fi    

narg=$#
input=()
for i in $* ;do
  input=("${input[@]}" $i)
done

gnu_file=plots.gnu
plot_dir=plots

read ipdf < ipdf.dat
read nch isingle_ch < nch.dat
nchp1=`expr $nch + 1`
read ndim < ndim.dat
nplot=$ndim
read chstart nchp1 < chrange.dat

k=$chstart
while [ $k -lt $nchp1 ]; do
    i=0
    while [ $i -lt $nplot ];do
        j=`expr $i + 1`
        plot[$i]="z${j}_$k"
        xlabel[$i]="z$j"
        i=`expr $i + 1`
    done

    nvplots=`expr $nplot \/ 2`
    nvplots=`expr $nvplots + 1`
    height=`expr 3 \* $nvplots`
    echo "set terminal postscript eps enhanced 'Times-Roman' color 17 size 8in,${height}in" > $gnu_file
    echo "#set logscale x" >> $gnu_file
    echo "#set logscale y" >> $gnu_file 
    echo "#set format x '%L'" >> $gnu_file
    echo "#set format y '10^{%L}'" >> $gnu_file
    echo "#set xtics (0.001,0.01)" >> $gnu_file
    echo "#set ytics (1,10,1E2)" >> $gnu_file
    echo "#set tics scale 2" >> $gnu_file
    echo "set grid" >> $gnu_file
    echo "#set key at 1.0E3,1.0E7 samplen 2" >> $gnu_file
    echo "#set key spacing 1.5" >> $gnu_file
    echo "#set xlabel 'cost' offset -1,0" >> $gnu_file
    echo "#set ylabel 'log_{/=10 10} L (Mpc)' offset 1,0" >> $gnu_file
    echo "#set xrange [-1:1]" >> $gnu_file
    echo "#set yrange [1E-5:2E8]" >> $gnu_file
    echo "" >> $gnu_file
    echo "set output '${plot_dir}/intvar_$k.eps'" >> $gnu_file
    
    echo "set multiplot layout ${nvplots},2" >> $gnu_file
    j=0
    while [ $j -lt $nplot ];do
	echo "set xlabel '{/=25 ${xlabel[$j]}}' offset -1,0" >> $gnu_file
	echo "set ylabel '{/=25 dXsec / ${xlabel[$j]}}' offset 1,0" >> $gnu_file
	echo "plot \\" >> $gnu_file
	nargp1=`expr $narg + 1`
	i=1
	while [ $i -lt $nargp1 ];do
	    im1=`expr $i - 1`
	    echo "'rslt_${input[$im1]}/data/${plot[$j]}.dat' u 1:2 t 'key${i}' w histep lt $i lw 3, \\" >> $gnu_file
	    i=`expr $i + 1`
	done
	echo "'rslt_${input[0]}/data/${plot[$j]}.dat' every::0::0 notitle w histep lt 1 lw 0" >> $gnu_file
	j=`expr $j + 1` 
    done
    
#echo "set nomultiplot" >> $gnu_file
    echo "unset multiplot" >> $gnu_file
    echo "" >> $gnu_file
    echo "reset" >> $gnu_file
    
    gnuplot plots.gnu

    k=`expr $k + 1`
done