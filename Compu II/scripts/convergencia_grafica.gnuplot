set terminal pngcairo size 800,600
set output '../grafica/convergencia.png'
set xlabel 'n (numero de intervalos)'
set ylabel 'Error absoluto |Vnum - Vana|'
set title 'Convergencia de Simpson vs Trapecio'
set logscale x
set logscale y
set grid
set key top right
plot '../resultados/convergencia.dat' using 1:4 with lines title 'Simpson', '../resultados/convergencia.dat' using 1:5 with lines title 'Trapecio'
