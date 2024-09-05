set terminal qt size 900,540 font "Times Roman,10"
set xtics
set xtics nomirror
set origin 0,0;set size 1,1; 
#set yrange[-0.5:7]
set title "Susceptibilidad magnética frente a Tª" font "Times Roman,14"
set ylabel "Chi"font "Times Roman, 12";set xlabel "Temperatura" font "Times Roman, 12"
set key
# Plot the data
plot './apartados16/chi.dat' using 1:($2/16*8) t 'N=16' ,\
'./apartados/chi.dat' using 1:($2/32*8) t 'N=32',\
'./apartados64/chi.dat' using 1:($2/64*8) t 'N=64'

