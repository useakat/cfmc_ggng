set terminal postscript eps enhanced 'Times-Roman' color 17
#set logscale x
set logscale y
#set format x '%L'
set format y '%L'
#set xtics (0.001,0.01)
#set ytics (1,10,1E2)
#set tics scale 2
set grid
#set key at 1.0E3,1.0E7 samplen 2
#set key spacing 1.5
#set xlabel 'cost' offset -1,0
#set ylabel 'log_{/=10 10} L (Mpc)' offset 1,0
#set xrange [-1:1]
#set yrange [1E-5:2E8]

set output 'wgted_plots/plots/pt_em.eps'
set xlabel 'pT_{e-}' offset -1,0
set multiplot
plot \
'Events/me.1m/pt_em.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/pt_em.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/pt_em.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/pt_ep.eps'
set xlabel 'pT_{e+}' offset -1,0
set multiplot
plot \
'Events/me.1m/pt_ep.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/pt_ep.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/pt_ep.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/pt_mum.eps'
set xlabel 'pT_{mu-}' offset -1,0
set multiplot
plot \
'Events/me.1m/pt_mum.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/pt_mum.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/pt_mum.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/pt_mup.eps'
set xlabel 'pT_{mu+}' offset -1,0
set multiplot
plot \
'Events/me.1m/pt_mup.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/pt_mup.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/pt_mup.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/y_em.eps'
set xlabel 'y_{e-}' offset -1,0
set multiplot
plot \
'Events/me.1m/y_em.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/y_em.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/y_em.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/y_ep.eps'
set xlabel 'y_{e+}' offset -1,0
set multiplot
plot \
'Events/me.1m/y_ep.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/y_ep.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/y_ep.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/y_mum.eps'
set xlabel 'y_{mu-}' offset -1,0
set multiplot
plot \
'Events/me.1m/y_mum.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/y_mum.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/y_mum.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/y_mup.eps'
set xlabel 'y_{mu+}' offset -1,0
set multiplot
plot \
'Events/me.1m/y_mup.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/y_mup.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/y_mup.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/dr_ee.eps'
set xlabel 'dR_{ee}' offset -1,0
set multiplot
plot \
'Events/me.1m/dr_ee.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/dr_ee.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/dr_ee.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/dr_e1mu1.eps'
set xlabel 'dR_{e1mu1}' offset -1,0
set multiplot
plot \
'Events/me.1m/dr_e1mu1.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/dr_e1mu1.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/dr_e1mu1.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/dr_e1mu2.eps'
set xlabel 'dR_{e1mu2}' offset -1,0
set multiplot
plot \
'Events/me.1m/dr_e1mu2.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/dr_e1mu2.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/dr_e1mu2.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/dr_e2mu1.eps'
set xlabel 'dR_{e2mu1}' offset -1,0
set multiplot
plot \
'Events/me.1m/dr_e2mu1.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/dr_e2mu1.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/dr_e2mu1.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/dr_e2mu2.eps'
set xlabel 'dR_{e2mu2}' offset -1,0
set multiplot
plot \
'Events/me.1m/dr_e2mu2.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/dr_e2mu2.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/dr_e2mu2.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

set output 'wgted_plots/plots/dr_mumu.eps'
set xlabel 'dR_{mumu}' offset -1,0
set multiplot
plot \
'Events/me.1m/dr_mumu.dat' u 1:2 t 'key1' w histep lt 1 lw 3 ,\
'Events/ycut.1m/dr_mumu.dat' u 1:2 t 'key2' w histep lt 2 lw 3 ,\
'Events/me.1m/dr_mumu.dat' every::0::0 notitle w histep lt 1 lw 0
set nomultiplot

reset
