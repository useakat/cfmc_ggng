#!/bin/bash

run=$1
nch=$2

cd rslt_$run

if [ -e uwgteff.dat ];then
    rm -f uwgteff.dat
    touch uwgteff.dat
else
    touch uwgteff.dat
fi

nchp1=`expr $nch + 1`
total_ntry=0
total_nevent=0
i=1
while [ $i -ne $nchp1 ];do
    eff=`grep "Generation efficiency" log_spring_$i`
    eff=${eff##*=}
    eff=${eff%%P*}
    nevent=`grep "generated events" log_spring_$i`
    nevent=${nevent##*=}
    nevent2=`expr $nevent \* 100`
    ntry=`divide.sh $nevent2 $eff 0`
    miss=`grep "miss-generation" log_spring_$i`
    miss=${miss##*=}
    miss=${miss%%times}
    echo $i $nevent $eff $ntry $miss >> uwgteff.dat

    total_ntry=`expr $total_ntry + $ntry`
    total_nevent=`expr $total_nevent + $nevent`
    i=`expr $i + 1`
done
total_nevent2=`expr $total_nevent \* 100`
uwgteff=`divide.sh $total_nevent2 $total_ntry 3`
echo $total_nevent " events are generated."
#echo "unweighting efficiency: " $uwgteff "%" 
echo "" >> uwgteff.dat
echo "total unweighting efficiency= " $uwgteff >> uwgteff.dat

