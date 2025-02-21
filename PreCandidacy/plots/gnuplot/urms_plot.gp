load "paddim_v3_rotation.dat"
set terminal postscript enh col "Times-Roman,30"
set output "urms_plot.eps"
set ylabel "|{/Bold u}_h|_{rms}" rotate by 0
set xlabel "1/Ro"
set log xy
set xrange [0.4:5]
#set yrange [1:9]
set key top left
set key spacing 1.3 font ",20"
#set arrow from 5.5,1 to 5.5,9 nohead dt 2 lw 2 lc rgb "blue"
#set arrow from 10,1 to 10,9 nohead dt 2 lw 2 lc rgb "red"   

plot sample\
 [i=1:5] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:6] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1", \
 [i=13:15] '+' using (invRo[i]):(sqrt(urms[i]**2 + vrms[i]**2)):(sqrt(uerr[i]**2 + verr[i]**2)) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle, \
 [1:4] 2*x dt 2 lw 2 lc rgb "dark-blue" title "2Ro^{-1}"
#[1:3] .5*x**(-1./2.) dt 2 lw 2 lc rgb "coral" title "Ro^{1/2}"
