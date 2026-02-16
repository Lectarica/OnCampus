set terminal pngcairo size 1000,700
set output 'band_n1p0000_T0p0200_t1p0000_tpm0p2000_mum0p4759.png'
set xlabel 'k-path'
set ylabel 'E(k) = eps(k) - mu'
set grid
set title 'Band (OpenMP): t-t\047 square, n1p0000_T0p0200_t1p0000_tpm0p2000, mum0p4759'
set arrow from 3.14159, graph 0 to 3.14159, graph 1 nohead dt 2
set arrow from 6.28319, graph 0 to 6.28319, graph 1 nohead dt 2
set xtics ('{/Symbol G}' 0, 'X' 3.14159, 'M' 6.28319, '{/Symbol G}' 10.7261)
plot 'band_n1p0000_T0p0200_t1p0000_tpm0p2000_mum0p4759.dat' using 1:2 with lines title 'E(k)'
