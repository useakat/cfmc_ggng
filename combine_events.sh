#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "combine_events.sh [run]"
    echo ""
    exit
fi
run=$1
read nch isingle_ch < nch.dat
chp1=`expr $nch + 1`

outfile=${run}_unweighted_events.lhe

cd rslt_$run
#sed -e '/<\/LesHouchesEvents>/d' ${run}_1.lhe > $outfile
sed -e '/<\/LesHouchesEvents>/d' 1.lhe > $outfile
i=2
while [ $i -lt $chp1 ];do
    sed -e '/<LesHouchesEvents version="1.0">/,/ <\/init>/d' \
	-e '/<\/LesHouchesEvents>/d' ${i}.lhe > tmp.lhe
    cat tmp.lhe >> $outfile
    i=`expr $i + 1`
done
rm -rf tmp.lhe
echo "</LesHouchesEvents>" >> $outfile

mv ${outfile} ../Events/.
cd ..
cd Events
gzip -f $outfile
rm -rf $outfile