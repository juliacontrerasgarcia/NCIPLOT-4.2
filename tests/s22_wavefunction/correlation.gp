reset
set encoding iso_8859_1
set title 'Integral vs. IE'
set term pngcairo font "sans,24" enhanced size 1240,1240
set output 'correlation.png'

set key autotitle columnheader
set key horiz
set key inside box opaque vert right top maxrows 5
set mxtics
unset mx2tics
unset my2tics
set mytics
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set style line 6 lt rgb '#006600' lw 1 pt 6 pi 5
set style line 15 lt rgb '#006600' lw 1 pt 2 pi 5
set tics nomirror
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid ls 12 back
set style fill transparent solid 0.3 noborder
set format x "%.2f"
set format y "%.2f"

f1(x) = m1*x + n1
stats 'correlation.dat' u (-1*$2):3 name "S1" 
fit f1(x) 'correlation.dat' u (-1*$2):3 via m1,n1
set label 1 sprintf("r = %4.2f",S1_correlation)  at graph 0.1, graph 0.85

plot 'correlation.dat' u (-1*$2):3 w p title 'n=2.0' pt 1 ps 1.5 lc rgb "red",\
 f1(x) w l ls 1 lc rgb "red" lw 2 title 'Lineal fit'



set output 'correlation.png' 
replot

exit


