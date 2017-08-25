reset
set dgrid3d 30,30
set pm3d map

set xrange [0:300]
set yrange [0:300]
set term png size 1500,500

set output 'output/map/plasma_'.i.'.png'
set multiplot layout 1,3

set title "Rho"
set cbrange [-7e-19:7e-19]
splot 'output/fields/rho_'.i.'.txt'

set title "Phi"
set cbrange [-4e-16:4e-16]
splot 'output/fields/phi_'.i.'.txt'

set title "E"
plot 'output/fields/E_'.i.'.txt' u 1:2:(3e10*$3):(3e10*$4) with vectors

print i

unset multiplot
