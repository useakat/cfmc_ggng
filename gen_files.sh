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

cd find_cf

make clean
./make_chinfo.sh $ngluons
./make_find_perch.sh $ngluons
make find_perch
./find_perch
cp opt_info.inc ../.
cd ..

./make_set_histo.sh $ngluons
./make_fill_histo.sh $ngluons
./make_save_histo.sh

make clean
make makefiles
./makefiles $ngluons

nfin=`expr $ngluons - 2`
cp -rf temp/matrix_mg_gg${nfin}g.f matrix_mg.f
cp -rf temp/matrix_mg_per_gg${nfin}g.f matrix_mg_per.f
cp -rf temp/matrix_co1_gg${nfin}g.f matrix_co1.f
cp -rf temp/flow1_gg${nfin}g.f flow1.f
cp -rf temp/matrix_cfqcd_gg${nfin}g.f matrix_cfqcd.f
cp -rf temp/icf_ch_gg${nfin}g.inc icf_ch.inc
cp -rf temp/plot_info_gg${nfin}g.sh plot_info.sh
cp -rf temp/plot_info_bases_gg${nfin}g.sh plot_info_bases.sh
cp -rf temp/makedat_gg${nfin}g.sh makedat.sh
