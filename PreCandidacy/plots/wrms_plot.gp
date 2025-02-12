load "paddim_v3_rotation.dat"
set terminal postscript enh col "Times-Roman,30"
set output "wrms_plot.eps"
set xlabel "1/Ro"
set ylabel "w_{rms}" rotate by 0
set key  bottom left spacing 1.3 font ",20"
set log xy
set xrange [0.4:5]
set yrange [0.1:1]
#set arrow from 5.5,0.03 to 5.5,1 nohead dt 2 lw 2 lc rgb "blue"
#set arrow from 10,0.03 to 10,1 nohead dt 2 lw 2 lc rgb "red"   

plot sample \
 [i=1:5] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "blue" title "Fr = 0.18",\
 [i=6:6] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" notitle,\
 [i=9:12] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 7 ps 2 lw 2 lc rgb "red" title "Fr = 0.1",\
 [i=13:15] '+' using (invRo[i]):(wrms[i]):(werr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" notitle,\
# [i=1:8] '+' using (invRo[i]):(wturb[i]):(wterr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" title "w_{turb}",\
# [i=1:8] '+' using (invRo[i]):(wlam[i]):(wlerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" title "w_{lam}",\
#[3:5] 1*x**(-1) dt 2 lw 2 lc rgb "black" title "Ro"
#[i=8:12] '+' using (invRo[i]):(wturb[i]):(wterr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "blue" title "w_{turb}"
#[i=8:12] '+' using (invRo[i]):(wlam[i]):(wlerr[i]) w yerrorbars pt 6 ps 2 lw 2 lc rgb "red" title "w_{lam}",\

