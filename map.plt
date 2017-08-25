reset
set dgrid3d 300,300
set pm3d map

set xrange [0:300]
set yrange [0:300]
set term png size 1000,500

set output 'output/particles_'.i.'.png'
set multiplot layout 1,2

set title "Electrons"
set cbrange [-1:0]
splot 'output/fields/el_'.i.'.txt' u 1:2:(-$3)

set title "Ions"
set cbrange [0:1]
splot 'output/fields/ion_'.i.'.txt' u 1:2:($3)

unset multiplot

print i
