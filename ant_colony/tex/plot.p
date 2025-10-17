set terminal png size 800,600
set output 'cpu.png'

unset label

#set title "Update"
set style line 1 lt 1 lw 3 pt 3 linecolor rgb "#000000"

set ylabel "CPU"
set xlabel "time"

set datafile separator ";"
set arrow 1 from 0,1549 to 60000,1549 nohead front lt 1 lw 2

plot "ressim2.dat" using 1:3 title 'cpu' with lines

#set output 'rejectrate.png'

set terminal pdf
set output "rejectrate.pdf"

set ylabel "reject rate"
#set format x "10^{%L}"
#set logscale x
plot "ressim2.dat" using 1:2 title 'reject rate' with lines 

#, "greedy.dat" using 1:2 title 'greedy reeject rate' with lines

#set output 'band.png'

#plot "sim2sub3band.dat" using 1:4 title 'bandwidth' with lines
