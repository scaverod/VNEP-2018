set terminal png size 800,600
set output 't1gap.png'

unset label

#set title "Update"
set style line 1 lt 1 lw 3 pt 3 linecolor rgb "#000000"

set xlabel "Size 15*S + V"
set ylabel "gap(%)"

set yrange [0:2]

plot "tlarge.dat" u ($1*15 + $2):($6/$7)  t 'gap' w linespoints

set output 't2gap.png'
plot "tlarge2.dat" u ($1*15 + $2):($6/$7)  t 'gap' w linespoints
set output 't3gap.png'
plot "tlarge3.dat" u ($1*15 + $2):($7/$8)  t 'gap' w linespoints
set output 't4gap.png'
plot "tlarge4.dat" u ($1*15 + $2):($7/$8)  t 'gap' w linespoints

set xlabel "Size 15*S + V"
set ylabel "time(s)"

set yrange [0:1000]

set output 't2time.png'
plot "tlarge.dat" u ($1*15 + $2):4  t 'cplex' w linespoints, "tlarge.dat" u ($1*15 + $2):5  t 'cgu' w linespoints
set output 't2time2.png'
plot "tlarge2.dat" u ($1*15 + $2):3  t 'cplex' w linespoints, "tlarge2.dat" u ($1*15 + $2):5  t 'cgu' w linespoints
set output 't2time3.png'
plot "tlarge3.dat" u ($1*15 + $2):4  t 'cplex' w linespoints, "tlarge3.dat" u ($1*15 + $2):6  t 'cgu' w linespoints
set output 't2time4.png'
plot "tlarge4.dat" u ($1*15 + $2):4  t 'cplex' w linespoints, "tlarge4.dat" u ($1*15 + $2):6  t 'cgu' w linespoints
