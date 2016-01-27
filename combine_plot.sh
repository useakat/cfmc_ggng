#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: combine_plot.sh [run]"
    echo ""
    exit
fi

run=$1

./combine_events.sh $run
./maplot.sh $run
./makedat.sh $run