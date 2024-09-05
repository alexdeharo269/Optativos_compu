set terminal qt size 400,300 font "Times Roman,9"
set xtics
set xtics nomirror
set origin 0,0;set size 1,1; 
set yrange[-1.1:1.1]
set title "Magnetización por dominios" font "Times Roman,13"
set ylabel "m"font "Times Roman, 11";set xlabel "Temperatura" font "Times Roman, 11"
set key right bottom 
set multiplot

plot './apartados75/magnet_domains.dat' using 1:2 t 'Dominio positivo' lc rgb 'red',\
'./apartados75/magnet_domains.dat' using 1:3 t 'Dominio negativo' lc rgb 'blue',\
0.5