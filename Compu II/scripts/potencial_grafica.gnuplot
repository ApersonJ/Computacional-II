set terminal pngcairo size 800,600
set output '../grafica/potencial_grafica.png'
set xlabel 'r'
set ylabel 'z'
set zlabel 'V(r,z)'
set title 'Potencial electrico'
set dgrid3d 50,50
set pm3d
set hidden3d
set view 60, 30, 1, 1
splot '../resultados/resultados_integracion.dat' using 1:2:3 with pm3d
