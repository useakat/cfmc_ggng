#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: fortran_temp.sh [type]"
    echo ""
    echo " [type] template type: pro or sub"
    echo ""
    exit
fi

next=$1
nini=$2
nfin=$3
ipdf=$4

out_file="Cards/process.dat"
#if [ -e $out_file ];then
#    rm -rf $out_file
#fi
#touch $out_file

if [ $ipdf -eq 0 ];then
    ini1=21
    ini2=21
else
    ini1=2212
    ini2=2212
fi

echo $ini1 $ini2 1 1 > $out_file
echo 21 -1 0 0 >> $out_file
echo 21 -1 0 0 >> $out_file
i=1
while [ $i -le $nfin ];do
    echo 21 1 1 2 >> $out_file
    i=`expr $i + 1`
done
    



    