#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "gg4g_combine_ch.sh [run]"
    echo ""
    exit
fi
run=$1

run_dir=run_$run

nch1=24
nch2=36
nch3=48
nch4=60
nch5=66

nch1p1=`expr $nch1 + 1`
nch2p1=`expr $nch2 + 1`
nch3p1=`expr $nch3 + 1`
nch4p1=`expr $nch4 + 1`
nch5p1=`expr $nch5 + 1`

outfile_t=${run}_t_unweighted_events.lhe
outfile_1s1=${run}_1s1_unweighted_events.lhe
outfile_1s2=${run}_1s2_unweighted_events.lhe
outfile_1s3=${run}_1s3_unweighted_events.lhe
outfile_2s=${run}_2s_unweighted_events.lhe

cd ${run_dir}

sed -e '/<\/LesHouchesEvents>/d' ${run}_1.lhe > $outfile_t
i=2
while [ $i -ne $nch1p1 ];do
    sed -e '/<LesHouchesEvents version="1.0">/,/ <\/init>/d' \
	-e '/<\/LesHouchesEvents>/d' ${run}_$i.lhe > tmp.lhe
    cat tmp.lhe >> $outfile_t
    i=`expr $i + 1`
done
rm -rf tmp.lhe
echo "</LesHouchesEvents>" >> $outfile_t

sed -e '/<\/LesHouchesEvents>/d' ${run}_${nch1p1}.lhe > $outfile_1s1
i=`expr $nch1p1 + 1`
while [ $i -ne $nch2p1 ];do
    sed -e '/<LesHouchesEvents version="1.0">/,/ <\/init>/d' \
	-e '/<\/LesHouchesEvents>/d' ${run}_$i.lhe > tmp.lhe
    cat tmp.lhe >> $outfile_1s1
    i=`expr $i + 1`
done
rm -rf tmp.lhe
echo "</LesHouchesEvents>" >> $outfile_1s1

sed -e '/<\/LesHouchesEvents>/d' ${run}_${nch2p1}.lhe > $outfile_1s2
i=`expr $nch2p1 + 1`
while [ $i -ne $nch3p1 ];do
    sed -e '/<LesHouchesEvents version="1.0">/,/ <\/init>/d' \
	-e '/<\/LesHouchesEvents>/d' ${run}_$i.lhe > tmp.lhe
    cat tmp.lhe >> $outfile_1s2
    i=`expr $i + 1`
done
rm -rf tmp.lhe
echo "</LesHouchesEvents>" >> $outfile_1s2

sed -e '/<\/LesHouchesEvents>/d' ${run}_${nch3p1}.lhe > $outfile_1s3
i=`expr $nch3p1 + 1`
while [ $i -ne $nch4p1 ];do
    sed -e '/<LesHouchesEvents version="1.0">/,/ <\/init>/d' \
	-e '/<\/LesHouchesEvents>/d' ${run}_$i.lhe > tmp.lhe
    cat tmp.lhe >> $outfile_1s3
    i=`expr $i + 1`
done
rm -rf tmp.lhe
echo "</LesHouchesEvents>" >> $outfile_1s3

sed -e '/<\/LesHouchesEvents>/d' ${run}_${nch4p1}.lhe > $outfile_2s
i=`expr $nch4p1 + 1`
while [ $i -ne $nch5p1 ];do
    sed -e '/<LesHouchesEvents version="1.0">/,/ <\/init>/d' \
	-e '/<\/LesHouchesEvents>/d' ${run}_$i.lhe > tmp.lhe
    cat tmp.lhe >> $outfile_2s
    i=`expr $i + 1`
done
rm -rf tmp.lhe
echo "</LesHouchesEvents>" >> $outfile_2s

mv $outfile_t ../Events/.
mv $outfile_1s1 ../Events/.
mv $outfile_1s2 ../Events/.
mv $outfile_1s3 ../Events/.
mv $outfile_2s ../Events/.
cd ..
cd Events
gzip $outfile_t
gzip $outfile_1s1
gzip $outfile_1s2
gzip $outfile_1s3
gzip $outfile_2s
rm -rf $outfile_t $outfile_1s1 $outfile_1s2 $outfile_1s3 $outfile_2s  