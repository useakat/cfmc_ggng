#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: mkplotdata.sh [run_name]"
    echo ""
    exit
fi
if [[ $1 = "" ]]; then
    echo "input run name"
    read run
else
    run=$1
fi    

maplot.sh $run
./makedat.sh $run