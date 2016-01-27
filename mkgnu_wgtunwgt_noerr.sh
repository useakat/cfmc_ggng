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

read ndist < ndist.dat
gnu_file=plots.gnu
plot_dir=plots
nplot=$ndist

i=0
while [ $i -ne $nplot ];do
    j=`expr $i + 1`
    plot_wgt[$i]="dist$j"
    i=`expr $i + 1`
done

source ./plot_info.sh

nvplots=`expr $nplot \/ 2`
nvplots=`expr $nvplots + 1`
height=`expr 3 \* $nvplots`
echo "set terminal postscript eps enhanced 'Times-Roman' color 17" > $gnu_file
echo "#set logscale x" >> $gnu_file
echo "set logscale y" >> $gnu_file 
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
j=0
while [ $j -ne $nplot ];do
    echo "set output '${plot_dir}/plot_${plot_unwgt[$j]}.eps'" >> $gnu_file
    echo "set multiplot" >> $gnu_file
    echo "set xlabel '{/=25 ${xlabel[$j]}}' offset -1,0" >> $gnu_file
    echo "set ylabel '{/=25 dXsec / d${xlabel[$j]}}' offset 1,0" >> $gnu_file
    echo "plot \\" >> $gnu_file
    nargp1=`expr $narg + 1`
    i=1
    while [ $i -ne $nargp1 ];do
	im1=`expr $i - 1`
	echo "'bases_plots/${input[$im1]}/${plot_wgt[$j]}.dat' u (\$1+${xshift[$j]}):2 t '${input[$im1]}_{wgt}' w histep lt $i lw 4, \\" >> $gnu_file
	echo "'Events/${input[$im1]}/${plot_unwgt[$j]}.dat' u 1:(\$2/${binsize[$j]}) t '${input[$im1]}_{unwgt}' w histep lt $i lw 2, \\" >> $gnu_file
	i=`expr $i + 1`
    done
    echo "'bases_plots/${input[0]}/${plot_wgt[$j]}.dat' every::0::0 notitle w histep lt 1 lw 0" >> $gnu_file
#echo "set nomultiplot" >> $gnu_file
    echo "unset multiplot" >> $gnu_file
    j=`expr $j + 1` 
done

echo "" >> $gnu_file
echo "reset" >> $gnu_file

gnuplot plots.gnu