#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: clean_run.sh [run_name]"
    echo ""
    exit
fi

run=$1
rm -rf Events/${run}_plots.html Events/$run rslt_$run Events/${run}_unweighted_events.lhe.gz
