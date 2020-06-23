# Gnuplot script for mapping NCI color code over NCI diagrams
set terminal pngcairo enhanced 
set encoding iso_8859_1
set output 'PLACEHOLDER.png' 
set key 
set ylabel 's(a.u.)' font "Helvetica, 30" 
set xlabel 'sign({/Symbol l}_2){/Symbol r}(a.u.)' font "Helvetica, 30"
set pm3d map
# Define a color gradient palette used by pm3d
set palette defined (-0.04 "blue",0.00 "green", 0.04 "red")
set format y "% .2f"
set format x "% .2f"
set format cb "% -.2f"
set border lw 4
set xtic  -0.06,0.01,0.06 nomirror rotate font "Helvetica"
set ytic   0.0,0.25,1.0 nomirror font "Helvetica"
# set the color bar tics
set cbtic  -0.06,0.01,0.06 nomirror font "Helvetica"
set xrange [-0.06:0.06]  
set yrange [0.0:1.0] 
# set the range of values which are colored using the current palette
set cbrange [-0.06:0.06]
plot 'PLACEHOLDER.dat' u 1:2:1 w p lw 6 palette t '' 
