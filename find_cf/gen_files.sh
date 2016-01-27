#!/bin/bash
if [[ "$1" == "-h" ]]; then
    echo ""
    echo "Usage: fortran_temp.sh [type]"
    echo ""
    echo " [type] template type: pro or sub"
    echo ""
    exit
fi

ngluons=$1

make clean
./make_chinfo.sh $ngluons
./make_find_perch.sh $ngluons
make find_perch
./find_perch
    