reset

#iterations=99

nVel=100 #number of intervals
nX=300
nY=300
nB=498
multFactor=200000
velFactor=50
#max = 150 #0.00001

maxVel=1e7 #max value
minVel=0 #min value
maxX=300
minX=0
maxY=300
minY=0
minB=0
maxB=498
widthVel=(maxVel-minVel)/nVel #interval width
widthX=(maxX-minX)/nX
widthY=(maxY-minY)/nY
widthB=1

#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term pngcairo size 1500,1000
set yrange [0:]

#to put an empty boundary around the
#data inside an autoscaled graph.
set offset graph 0.05,0.05,0.05,0.0
set style fill solid 0.5 #fillstyle
set tics out nomirror
#set xlabel "v"
set ylabel "entries"

#f(x) = sqrt((mk/(2*pi*T1))**3)*4*pi*((x-B)**2)*exp((-mk*((x-B)**2))/(2*T1))
#g(x) = sqrt((mk/(2*pi*T2))**3)*4*pi*((x-D)**2)*exp((-mk*((x-D)**2))/(2*T2))

#f(x) = A*((x-B)**2)*exp((-mk*((x-B)**2))/(2*T1))
#g(x) = C*((x-D)**2)*exp((-mk*((x-D)**2))/(2*T2))

#m = 9.1094e-31
#k = 1.38064852e-23
#mk = 6.5979e-9
#T1 = 1e7
#T2 = 1e7
#A = sqrt((mk/(2*pi*T1))**3)*4*pi
#C = sqrt((mk/(2*pi*T1))**3)*4*pi
#B = 1e6
#D = 1e6

#fit f(x) "electronsStart.txt" u (hist(sqrt(($3**2)+($4**2)+($5**2)),width)):(1.0) via T1
#fit g(x) "electronsStop.txt" u (hist(sqrt(($3**2)+($4**2)+($5**2)),width)):(1.0) via T2

#count and plot
set output 'output/hist/stats_'.i.'.png'
set multiplot layout 2,4

set xrange [minX:maxX]
set xtics minX,(maxX-minX)/5,maxX
set boxwidth widthX*0.9
set title "x_{el}"
plot 'output/particles/el_'.i.'.txt' u (hist($1*multFactor,widthX)):(1.0) smooth freq w boxes notitle

set xrange [minY:maxY]
set xtics minY,(maxY-minY)/5,maxY
set boxwidth widthY*0.9
set title "y_{el}"
plot 'output/particles/el_'.i.'.txt' u (hist($2*multFactor,widthY)):(1.0) smooth freq w boxes notitle

set xrange [minB:maxB]
set xtics minB,(maxB-minB)/5,maxB
set boxwidth widthB*0.9
set title "current_{el}"
plot 'output/fields/elIntArr_'.i.'.txt' u 1:2 smooth freq w boxes notitle

set xrange [minVel:maxVel]
set xtics minVel,(maxVel-minVel)/5,maxVel
set boxwidth widthVel*0.9
set title "v^2_{el}"
plot 'output/particles/el_'.i.'.txt' u (hist(sqrt(($4**2)+($5**2)+($6**2)),widthVel)):(1.0) smooth freq w boxes notitle

set xrange [minX:maxX]
set xtics minX,(maxX-minX)/5,maxX
set boxwidth widthX*0.9
set title "x_{ion}"
plot 'output/particles/ion_'.i.'.txt' u (hist($1*multFactor,widthX)):(1.0) smooth freq w boxes notitle

set xrange [minY:maxY]
set xtics minY,(maxY-minY)/5,maxY
set boxwidth widthY*0.9
set title "y_{ion}"
plot 'output/particles/ion_'.i.'.txt' u (hist($2*multFactor,widthY)):(1.0) smooth freq w boxes notitle

set xrange [minB:maxB]
set xtics minB,(maxB-minB)/5,maxB
set boxwidth widthB*0.9
set title "current_{ion}"
plot 'output/fields/ionIntArr_'.i.'.txt' u 1:2 smooth freq w boxes notitle

set xrange [minVel:maxVel]
set xtics minVel,(maxVel-minVel)/5,maxVel
set boxwidth widthVel*0.9
set title "v^2_{ion}"
plot 'output/particles/ion_'.i.'.txt' u (hist(velFactor*sqrt(($4**2)+($5**2)+($6**2)),widthVel)):(1.0) smooth freq w boxes notitle

unset multiplot

print i
