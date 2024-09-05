set terminal qt size 900,540 font "Times Roman,10"
set xtics
set xtics nomirror
set origin 0,0;set size 1,1; 
set yrange[-1.1:1.1]
set title "Magnetizacón por dominios" font "Times Roman,14"
set ylabel "Magnetizacón inferior     Magnetizacón superior"font "Times Roman, 12";set xlabel "Temperatura" font "Times Roman, 12"
set multiplot
# Plot the data
plot './apartados16/magnet_domains.dat' using 1:2 t 'N=16' ,\
'./apartados/magnet_domains.dat' using 1:2 t 'N=32',\
'./apartados64/magnet_domains.dat' using 1:2 t 'N=64'

plot './apartados16/magnet_domains.dat' using 1:3 t 'N=16' ,\
'./apartados/magnet_domains.dat' using 1:3 t 'N=32',\
'./apartados64/magnet_domains.dat' using 1:3 t 'N=64'
unset key
plot 0 #linea horizontal en 0
