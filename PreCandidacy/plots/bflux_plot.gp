load "paddim_v3_rotation.dat"
set terminal postscript enh col "Times-Roman,30"
set output "bflux_plot.eps"
set xlabel "1/Ro"
set ylabel "|wT|/(Fr^2)" rotate by 0
set format y "10^{%T}"
set log xy
set xrange [0.4:5]
set yrange [0.01:1]
set key bottom left spacing 1.3 font ",20"
#set arrow from 5.5,0.00001 to 5.5,0.02 nohead dt 2 lw 2 lc rgb "blue"
#set arrow from 10,0.00001 to 10,0.02 nohead dt 2 lw 2 lc rgb "red"   

plot sample\
 [i=1:5] '+' using (invRo[i]):(30*tflux[i]):(tferr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:6] '+' using (invRo[i]):(30*tflux[i]):(tferr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(100*tflux[i]):(tferr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1",\
 [i=13:13] '+' using (invRo[i]):(100*tflux[i]):(tferr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle
#[0.5:3] .24*x**(1./4.) dt 2 lw 2 lc rgb "coral" title "(1/Ro)^{1/4}"
#[3:5] 1.08*x**(-1) dt 2 lw 2 lc rgb "dark-violet" title "Ro",\
