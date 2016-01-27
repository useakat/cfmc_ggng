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

rm -rf plots/*

narg=$#
input=()
for i in $* ;do
  input=("${input[@]}" $i)
done

gnu_file=plots.gnu
plot_dir=plots
read ndist < ndist.dat
nplot=$ndist

plot=("pt_g1" "pt_g2" "pt_g3" "y_g1" "y_g2" "y_g3" "phi_g1" "phi_g2" "phi_g3" "dr_g1g2" "dr_g1g3" \
"dr_g2g3" "m_g1g2" "m_g1g3" "m_g2g3" "kt_g1g2" "kt_g1g3" "kt_g2g3")

xlabel=("pT (g1)" "pT (g2)" "pT (g3)" "y (g1)" "y (g2)" "y (g3)" "phi (g1)" "phi (g2)" "phi (g3)" \
"dR (g1g2)" "dR (g1g3)" "dR (g2g3)" "m (g1g2)" "m (g1g3)" "m (g2g3)" "kT (g1g2)" "kT (g1g3)" \
"kT (g2g3)")

echo "set terminal postscript eps enhanced 'Times-Roman' color 17" > $gnu_file
echo "#set logscale x" >> $gnu_file
echo "set logscale y" >> $gnu_file 
echo "#set format x '%L'" >> $gnu_file
echo "set format y '%L'" >> $gnu_file
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
    echo "set output '${plot_dir}/plot_${plot[$j]}.eps'" >> $gnu_file
    echo "set xlabel '${xlabel[$j]}' offset -1,0" >> $gnu_file
    echo "set multiplot" >> $gnu_file
    echo "plot \\" >> $gnu_file
    
    nargp1=`expr $narg + 1`
    i=1
    while [ $i -ne $nargp1 ];do
	im1=`expr $i - 1`
	echo "'Events/${input[$im1]}/${plot[$j]}.dat' u 1:2 t '${input[$im1]}' w histep lt $i lw 3 ,\\" >> $gnu_file
	i=`expr $i + 1`
    done
    echo "'Events/${input[0]}/${plot[$j]}.dat' every::0::0 notitle w histep lt 1 lw 0" >> $gnu_file
    echo "set nomultiplot" >> $gnu_file
    echo "" >> $gnu_file

    j=`expr $j + 1` 
done

echo "reset" >> $gnu_file

gnuplot plots.gnu