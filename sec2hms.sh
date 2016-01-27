#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: run.sh [sectime]"
    echo ""
    exit
fi
sectime=$1
HH=`expr ${sectime} / 3600` 
SS=`expr ${sectime} % 3600` 
MM=`expr ${SS} / 60` 
SS=`expr ${SS} % 60` 
echo "${HH}:${MM}:${SS}" 