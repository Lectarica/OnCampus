set terminal pngcairo size 850,850
set output 'fs_n1p0000_T0p0200_t1p0000_tpm0p2000_mum0p4759.png'
set size ratio -1
set xrange [-pi:pi]
set yrange [-pi:pi]
set xlabel 'k_x'
set ylabel 'k_y'
set grid
set title 'FS (OpenMP): n1p0000_T0p0200_t1p0000_tpm0p2000, mum0p4759'
pi=3.141592653589793
set xtics (-pi, -pi/2, 0, pi/2, pi)
set ytics (-pi, -pi/2, 0, pi/2, pi)
plot 'fs_n1p0000_T0p0200_t1p0000_tpm0p2000_mum0p4759.dat' using 1:2 with points pt 7 ps 0.25 title 'FS points'
