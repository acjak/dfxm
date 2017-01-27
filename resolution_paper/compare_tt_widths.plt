set term postscript eps color enhanced "Helvetica" 16

#set log y
#set format y "10^{%L}"
#set xrange [0.5:6]
#set yrange [0.0000001:1]

#set ylabel "Normalized intensity"
set xlabel "2{/Symbol q} (degrees)"


set output 'tt_widths.eps'
plot '../output/230117-1418/bot500.txt' using 1:2 with linespoints lt rgb "black" ti "500 px below center", \
'../output/230117-1418/fit_bot500.txt' using 1:2 w l lw 2 lt rgb "red" ti "Gaussian fit", \
'../output/230117-1418/top500.txt' using 1:2 with linespoints lt rgb "black" ti "500 px above center", \
'../output/230117-1418/fit_top500.txt' using 1:2 w l lw 2 lt rgb "blue" ti "Gaussian fit",\
'../output/230117-1418/center.txt' using 1:2 with linespoints lt rgb "black" ti "Center", \
'../output/230117-1418/fit_center.txt' using 1:2 w l lw 2 lt rgb "violet" ti "Gaussian fit" \
