
set terminal qt size 900,540 font "Times Roman,10"
set xtics
set xtics nomirror
unset grid
#set title  "Energía media por partícula" font "Times-Roman, 11"
set ylabel "Energía por partícula"font "Times Roman, 12";set xlabel "Temperatura" font "Times Roman, 12"
set origin 0,0;set size 1,1; set xrange[0:5];set yrange[-2:-0.2]
set key left top
set title "Energía por partícula frente a Tª" font "Times Roman, 14"
set multiplot# Plot the data
plot './apartados16/energy_per_particle.dat' using 1:2 t 'N=16' ,\
'./apartados/energy_per_particle.dat' using 1:2 t 'N=32',\
'./apartados64/energy_per_particle.dat' using 1:2 t 'N=64'

set origin 0.47,0.12
unset xlabel;unset ylabel; unset key; set yrange[0:-0.008];set size 0.5,0.6
set grid; #set key
unset title
set ylabel 'Error relativo'
plot './apartados16/energy_per_particle.dat' using 1:3 t 'N=16' ,\
'./apartados/energy_per_particle.dat' using 1:3 t 'N=32',\
'./apartados64/energy_per_particle.dat' using 1:3 t 'N=64'