set terminal pngcairo size 1000,700
set output 'dos_n1p0000_T0p0200_t1p0000_tpm0p2000_eta0p0200.png'
set xlabel 'E'
set ylabel 'DOS(E)  (per site, spin included)'
set grid
set title 'DOS (OpenMP): n1p0000_T0p0200_t1p0000_tpm0p2000, eta=0.02'
plot 'dos_n1p0000_T0p0200_t1p0000_tpm0p2000_eta0p0200.dat' using 1:2 with lines title 'DOS'
