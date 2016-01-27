#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: combine_lhe.sh [output] file] [lhe1] [lhe2] (lhe3) ..."
    echo "lhe1 is the file name part of lhe1_unweighted_events.lhe.gz"
    echo "Event numbers should be same for all the combined lhe files"
    exit
fi

ext="unweighted_events.lhe"

output=${1}_${ext}
nfile=`expr $# - 1`
input=()
for i in $* ;do
  input=("${input[@]}" $i)
done
i=1
while [ $i -le $nfile ];do
    gunzip -c ${input[$i]}_${ext}.gz > tmp_$i
    if [ $i -eq 1 ];then
	sed -e '$d' tmp_$i > tmp2_${i}
    elif [ $i -eq $nfile ];then
	sed -e '1,9d' tmp_$i > tmp2_${i}
    else
	sed -e '1,9d' -e '$d' tmp_$i > tmp2_${i}
    fi
    rm -rf tmp_$i
    i=`expr $i + 1`
done


# i=1
# while [ $i -le $nfile ];do
#     touch tmp3_${i}
#     event_flag=0
#     while read line ;do
# 	if [ $event_flag -eq 1 ];then
# 	    colm1=`echo ${line} | cut -d " " -f 1`
# 	    colm2=`echo ${line} | cut -d " " -f 2`
# 	    colm3=`echo ${line} | cut -d " " -f 3`
# 	    colm4=`echo ${line} | cut -d " " -f 4`
# 	    colm5=`echo ${line} | cut -d " " -f 5`
# 	    colm6=`echo ${line} | cut -d " " -f 6`
# 	    num=${colm3%%E*}
# 	    exp=${colm3##${num}}
# 	    colm3_2=`echo "scale=7; ${num}/${nfile}" | bc | sed 's/^\./0./'`
# 	    echo $colm1 $colm2 ${colm3_2}$exp $colm4 $colm5 $colm6 >> tmp3_${i}
# 	    event_flag=0
# 	else
# 	    echo $line >> tmp3_${i}
# 	fi
# 	if [ "$line" == "<event>" ];then
# 	    event_flag=1
# 	fi
#     done < tmp2_${i}    

#     i=`expr $i + 1`
# done


touch ${output}
i=1
while [ $i -le $nfile ];do
    cat tmp2_${i} >> ${output}
    i=`expr $i + 1`
done

gzip ${output}

rm -rf tmp2_* tmp3_*