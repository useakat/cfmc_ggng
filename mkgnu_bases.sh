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
nplot=`expr $ndist + 3`

plot=("${plot[@]}" "x1")
plot=("${plot[@]}" "x2")
plot=("${plot[@]}" "m2summ2")
i=3
while [ $i -ne $nplot ];do
    j=`expr $i - 2`
    plot[$i]="dist$j"
    i=`expr $i + 1`
done

# xlabel=("x1" "x2" "M2/SumM2" "pT (g1)" "pT (g2)" "pT (g3)" "pT (g4)" "y (g1)" "y (g2)" "y (g3)" "y (g4)" "phi (g1)" \
# "phi (g2)" "phi (g3)" "phi (g4)" "dR (g1g2)" "dR (g1g3)" "dR (g1g4)" "dR(g2g3)" "dR (g2g4)" "dR (g3g4)" "m (g1g2)" "m (g1g3)" \
# "m (g1g4)" "m (g2g3)" "m (g2g4)" "m (g3g4)" "kT (g1g2)" "kT (g1g3)" "kT (g1g4)" "kT (g2g3)" "kT (g2g4)" "kT (g3g4)")

source ./plot_info_bases.sh

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
echo "set output '${plot_dir}/plot_distw.eps'" >> $gnu_file

echo "set multiplot layout ${nvplots},2" >> $gnu_file
j=0
while [ $j -ne $nplot ];do
    if [ $j == "3" ]; then
	echo "set logscale y" >> $gnu_file
	echo "set format y '10^{%L}'" >> $gnu_file
    fi
    echo "set xlabel '{/=25 ${xlabel[$j]}}' offset -1,0" >> $gnu_file
    echo "set ylabel '{/=25 dXsec / ${xlabel[$j]}}' offset 1,0" >> $gnu_file
    echo "plot \\" >> $gnu_file
    nargp1=`expr $narg + 1`
    i=1
    while [ $i -ne $nargp1 ];do
	im1=`expr $i - 1`
	echo "'rslt_${input[$im1]}/data/${plot[$j]}.dat' u 1:2 t '${input[$im1]} (wgted)' w histep lt $i lw 3, \\" >> $gnu_file
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