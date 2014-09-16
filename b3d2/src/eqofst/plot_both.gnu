set terminal postscript eps noenhanced defaultplex
set output 'plot_T_and_P.eps'
set title "Pressure and Temperature vs Energy density"
set xlabel "energy density (Mev/fm^3)"
set ylabel "Mev"
plot "output_mesh.dat" using 1:2 with lines ti "Pressure", "output_mesh.dat" using 1:3 ti "Temperature"
