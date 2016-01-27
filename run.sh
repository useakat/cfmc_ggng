#!/bin/bash
####
####  Usage
####
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: run.sh [run_name] [number_of_events] [serial/parallel] [run_mode] [gen_mode] [mail]"
    echo ""
    echo " [serial/parallel] 0:serial 1:parallel"
    echo " [run_mode] 0:optim+int 1:optim 2:int"
    echo " [gen_mode] 2:weighted 3:unweighted"
    echo " [mail] whether inform when the job finishes or not, 0:no 1:inform"
    exit
fi

####
#### Input arguments
####
if [[ $1 = "" ]]; then
    echo "input run name"
    read run
    echo "input the number of events to be generated"
    read genevents
    echo "input parallel/serial  0:serial 1:parallel"
    read mode
    echo "input run mode 0:optim+int 1:optim 2:int"
    read type
    echo "inpyt event generation mode 2:weighted 3:unweighted"
    read genmode
    echo "input whether inform when the job finishes or not 0:no 1:inform"
    read mail
else
    run=$1
    genevents=$2
    mode=$3
    type=$4
    genmode=$5
    mail=$6
fi    

####
#### Cluster Setting
####
cluster=kekcc
#cluster=ingrid
#que=e # For jobs < 10 min
que=s # For jobs < 3 hours
#que=l # For jobs < 24 hours
#que=h # For heavy jobs

if [ $cluster == "kekcc" ]; then
    job_system=bsub
elif [ $cluster == "ingrid" ]; then
    job_system=qsub
fi

####
#### Initialization
####
make clean >/dev/null 2>&1
make 
make spring 
#rm -rf *.dat
rm -rf plots/*

start_time=`date '+%s'`
echo `date '+%T'`

echo $run
echo '      integer len_run_name' > len_run_name.inc
echo '      parameter (len_run_name='${#run}')' >> len_run_name.inc 

make get_numch 
./get_numch
read nch isingle_ch < nch.dat
read ipdf < ipdf.dat
read next nini nfin < next.dat
if [ $isingle_ch -eq 0 ];then
    chstart=1
    chp1=`expr $nch + 1`
else
    chstart=$isingle_ch
    chp1=`expr $chstart + 1`
fi

echo $chstart $chp1 > chrange.dat

run_dir=rslt_$run
out_dir=${run_dir}/output
data_dir=${run_dir}/data
bases_data_dir=bases_data

./make_process_dat.sh $next $nini $nfin $ipdf

####
#### Grid Optimization Phase
####
if [ $type == 0 -o $type == 1 ]; then 
    echo ""
    echo "Integrating ..."

### Setting output directory
    if [ -e $run_dir ]; then
	rm -rf ${run_dir}/*
	mkdir ${out_dir}
	mkdir ${data_dir}
    else
	mkdir ${run_dir}
	mkdir ${out_dir}
	mkdir ${data_dir}
    fi

### Job submission
    i=$chstart
    job=int
    jobname="int"$RANDOM
    while [ $i -lt $chp1 ]; do
#	cp -rf maxwgts${i}.dat maxwgts_ref${i}.dat
	ich=$i
	idiv=1
	seed=$RANDOM
##  Serial submission
	if [ $mode == 0 ]; then
	    ./run_bases $run $ich $idiv $seed
## Parallel submission
	elif [ $mode == 1 ]; then
	    echo "#!/bin/bash" > ajob$i
	    echo "#PBS -N $job" >> ajob$i
	    echo "start_time=\`date '+%s'\`" >> ajob$i
	    echo "rm -rf wait.ajob$i" >> ajob$i
	    echo "touch run.ajob$i" >> ajob$i
	    echo "./run_bases $run $ich $idiv $seed" >> ajob$i
	    echo "rm -rf run.ajob$i" >> ajob$i
	    echo "touch done.ajob$i" >> ajob$i   
	    echo "end_time=\`date '+%s'\`" >> ajob$i
	    echo "timediff=\`expr \${end_time} - \${start_time}\`" >> ajob$i
	    echo "HH=\`expr \${timediff} / 3600\`" >> ajob$i
	    echo "SS=\`expr \${timediff} % 3600\`" >> ajob$i
	    echo "MM=\`expr \${SS} / 60\`" >> ajob$i
	    echo "SS=\`expr \${SS} % 60\`" >> ajob$i
	    echo "elapsed_time=\${HH}:\${MM}:\${SS}" >> ajob$i
	    echo "echo \$timediff \$elapsed_time > btime_${ich}.dat" >> ajob$i
	    chmod +x ajob$i
	    touch wait.ajob$i 
	    if [ $job_system == "bsub" ]; then
		bsub -q $que -J $jobname ./ajob$i 1>/dev/null
	    elif [ $job_system == "qsub" ]; then
		qsub -N $jobname ./ajob$i 1>/dev/null
	    fi
	fi
	i=`expr $i + 1`
    done
### Monitoring the parallel jobs
    if [ $mode == 1 ]; then
	./monitor
	rm -rf *ajob*  
    fi

####
#### Combining the submitted job results
####
    make get_nevent >/dev/null 2>&1
    ./get_nevent $genevents $chstart $chp1
    read xsec_total < xsec_total.dat
    cat log_xsec.dat
    echo ""

    cat params.dat > ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    echo "Cross section" >> ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    cat log_xsec.dat >> ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    echo "Number of events to be generated" >> ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    cat log_nevent.dat >> ${run_dir}/summary.txt

    read ndim < ndim.dat
    read ndist < ndist.dat

    cd ${bases_data_dir}
    mv *.top ../${data_dir}/.
    cd ..

    ./makedat_bases.sh $run $nch $ndim $ndist $chstart $chp1
    make combine_dat >/dev/null 2>&1
    ./combine_dat $run $ndim $chstart $chp1
    ./mkgnu_z.sh $run
    if [ $type == 1 ]; then
	./mkgnu_bases.sh $run
    fi

####
####  Combining Time measurements
####
    if [ $mode -eq 1 ];then
	rm -rf btime.dat
	touch btime.dat
	sum_btime=0
	i=$chstart
	while [ $i -lt $chp1 ]; do
	    read etimes etime < btime_${i}.dat
	    echo $i $etimes $etime >> btime.dat
	    sum_btime=`expr $sum_btime + $etimes`
	    i=`expr $i + 1`
	done
	elapsed_time=`./sec2hms.sh $sum_btime`
	echo $sum_btime $elapsed_time >> btime.dat
	echo $sum_btime $elapsed_time > total_btime.dat
	echo "Serial integration time: " $elapsed_time
	echo
	echo "Serial integration time: " $elapsed_time >> ${run_dir}/summary.txt
	rm -rf btime_*
    fi

    mv xsec_* $run_dir
    mv nevt_* $run_dir
fi

####
####  Event Generation Phase
####
if [ $type == 0 -o $type == 2 ]; then 
    if [ $isingle_ch -ne 0 ];then
	echo "ERROR: isignle_ch is not 0"
	echo "exit event generation..."
	exit
    fi
    echo "Generating Events ..."

### Preparing the information for event generation
    cp -f ${run_dir}/xsec_* .
    cp -f ${run_dir}/nevt_* .
    make get_nevent >/dev/null 2>&1
    ./get_nevent $genevents $chstart $chp1

#### Initializing the events variable
    i=$chstart
    while [ $i -ne $chp1 ]; do
	events[$i]=`cat nevt_${i}.dat`
	i=`expr $i + 1`
    done

### Job submission
    i=$chstart
    while [ $i -ne $chp1 ]; do    
	ich=$i
	idiv=1
	seed=$RANDOM
	jobname="evt"$RANDOM
##  Serial submission
	if [ $mode == 0 ]; then
	    ./run_spring $run $ich $idiv $seed $genmode ${events[$i]}
## Parallel submission
	elif [ $mode == 1 ]; then
	    echo "#!/bin/bash" > ajob$i 
	    echo "#PBS -N $job" >> ajob$i
	    echo "start_time=\`date '+%s'\`" >> ajob$i
	    echo "rm -rf wait.ajob$i" >> ajob$i
	    echo "touch run.ajob$i" >> ajob$i  
	    echo "./run_spring $run $ich $idiv $seed $genmode ${events[$i]}" >> ajob$i
	    echo "rm -rf run.ajob$i" >> ajob$i
	    echo "touch done.ajob$i" >> ajob$i   
	    echo "end_time=\`date '+%s'\`" >> ajob$i
	    echo "timediff=\`expr \${end_time} - \${start_time}\`" >> ajob$i
	    echo "HH=\`expr \${timediff} / 3600\`" >> ajob$i
	    echo "SS=\`expr \${timediff} % 3600\`" >> ajob$i
	    echo "MM=\`expr \${SS} / 60\`" >> ajob$i
	    echo "SS=\`expr \${SS} % 60\`" >> ajob$i
	    echo "elapsed_time=\${HH}:\${MM}:\${SS}" >> ajob$i
	    echo "echo \$timediff \$elapsed_time > stime_${ich}.dat" >> ajob$i
	    chmod +x ajob$i
	    touch wait.ajob$i 
	    if [ $job_system == "bsub" ]; then
		bsub -q $que -J $jobname ./ajob$i 1>/dev/null
	    elif [ $job_system == "qsub" ]; then
		qsub -N $jobname ./ajob$i 1>/dev/null
	    fi
	fi
	i=`expr $i + 1`	
    done
### Monitoring the parallel jobs
    if [ $mode == 1 ]; then
	./monitor
	rm -rf *ajob*  
    fi

####
#### Combining the submitted job results
####
    echo "Combining events and creating plots"
    ./combine_plot.sh $run

    echo ""
    ./combine_uwgteff.sh $run $nch
    echo "" >> ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    echo "Unweighting efficiency" >> ${run_dir}/summary.txt
    echo "" >> ${run_dir}/summary.txt
    echo "ch. nevent  eff.  ntry  miss_gen." >> ${run_dir}/summary.txt
    cat ${run_dir}/uwgteff.dat >> ${run_dir}/summary.txt
    cat ${run_dir}/uwgteff.dat

    ./mkgnu_wgtunwgt_each2.sh ${run}

####
####  Combining Time measurements
####
    if [ $mode -eq 1 ];then
	rm -rf stime.dat
	touch stime.dat
	sum_stime=0
	i=$chstart
	while [ $i -lt $chp1 ]; do
	    read etimes etime < stime_${i}.dat
	    echo $i $etimes $etime >> stime.dat
	    sum_stime=`expr $sum_stime + $etimes`
	    i=`expr $i + 1`
	done
	elapsed_time=`./sec2hms.sh $sum_stime`
	echo $sum_stime $elapsed_time >> stime.dat
	echo $sum_stime $elapsed_time > total_stime.dat
	echo "Serial generation time: " $elapsed_time
	echo
	echo "Serial generation time: " $elapsed_time >> ${run_dir}/summary.txt
	echo "" >> ${run_dir}/summary.txt
	rm -rf stime_*
    fi

    rm -f xsec_* nevt_*
fi

####
#### Summarizing the final results
####
read xsec_total < ${run_dir}/xsec_total.dat
echo "Xsec:" $xsec_total "[pb]"
echo

if [ $mode -eq 1 ];then
    if [ $type -eq 0 ];then
	read totalbtime dum < total_btime.dat
	read totalstime dum < total_stime.dat
	serial_time=`expr $totalbtime + $totalstime`
    elif [ $type -eq 1 ];then
	read serial_time dum < total_btime.dat
    elif [ $type -eq 2 ];then
	read serial_time dum < total_stime.dat
    fi
    serial_hmstime=`./sec2hms.sh $serial_time`
    echo "Total serial time: " $serial_hmstime
    echo "" >> ${run_dir}/summary.txt
    echo "Total serial time: " $serial_hmstime >> ${run_dir}/summary.txt
fi
mv btime.dat ${out_dir}
mv stime.dat ${out_dir}
rm -rf total_btime.dat
rm -rf total_stime.dat

end_time=`date '+%s'`
SS=`expr ${end_time} - ${start_time}` 
elapsed_time=`./sec2hms.sh $SS`
echo "Elapsed time:" $elapsed_time

echo "Total time = " $elapsed_time >> ${run_dir}/summary.txt
echo

cp ${run_dir}/summary.txt ${out_dir}/. >/dev/null 2>&1
cp ${run_dir}/log_bases_* ${out_dir}/. >/dev/null 2>&1
cp ${run_dir}/log_spring_* ${out_dir}/. >/dev/null 2>&1

cp -rf plots/* ${out_dir}/.
cp -rf ${out_dir} .

rm -f log_xsec.dat log_nevent.dat

echo `date '+%T'`

####
####  Sending a notification mail for the completion of the calculations
####
if [ $mail -eq 1 ]; then
    if [ $cluster == "kekcc" ]; then
	bsub -q ${que} -J $run -u takaesu@post.kek.jp nulljob.sh >/dev/null 2>&1
    else
	echo "Notification mail cannot be send from this cluster system. Exting..."
	exit
    fi
fi