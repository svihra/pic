reset

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
set term png size 1200,600
set yrange [0:*]

#to put an empty boundary around the
#data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set style fill solid 0.5 #fillstyle
set tics out nomirror
#set xlabel "v"
set ylabel "entries"

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
set yrange [0.5:100]
set xtics minB,(maxB-minB)/5,maxB
set logscale y
set arrow from 100, graph 0 to 100, graph 1 nohead
set arrow from 199, graph 0 to 199, graph 1 nohead
set arrow from 299, graph 0 to 299, graph 1 nohead
set arrow from 398, graph 0 to 398, graph 1 nohead
set title "current_{el}"
plot 'output/fields/elIntArr_'.i.'.txt' u 1:($2) w boxes notitle
unset arrow
unset logscale y

set yrange [0:*]
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
set yrange [0.5:100]
set logscale y
set xtics minB,(maxB-minB)/5,maxB
set boxwidth widthB*0.9
set title "current_{ion}"
set arrow from 100, graph 0 to 100, graph 1 nohead
set arrow from 199, graph 0 to 199, graph 1 nohead
set arrow from 299, graph 0 to 299, graph 1 nohead
set arrow from 398, graph 0 to 398, graph 1 nohead
plot 'output/fields/ionIntArr_'.i.'.txt' u 1:($2) w boxes notitle
unset arrow
unset logscale y

set yrange [0:*]
set xrange [minVel:maxVel]
set xtics minVel,(maxVel-minVel)/5,maxVel
set boxwidth widthVel*0.9
set title "v^2_{ion}"
plot 'output/particles/ion_'.i.'.txt' u (hist(velFactor*sqrt(($4**2)+($5**2)+($6**2)),widthVel)):(1.0) smooth freq w boxes notitle

unset multiplot

print i
