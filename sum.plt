reset

tit1="electrons"
tit2="ions"

set logscale y
set term png size 500,500

set output './output/sum.png'

plot "./output/sum.txt" u 1:2 t tit1, "./output/sum.txt" u 1:3 t tit2
