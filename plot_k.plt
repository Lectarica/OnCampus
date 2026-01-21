set term pngcairo size 900,800
set output 'lorentzian_k_mag.png'

set size ratio -1
set xlabel 'k_y (shifted index)'
set ylabel 'k_x (shifted index)'
set title '|F(k)| (fftshifted)'

set pm3d map
set logscale cb
set cblabel '|F| (log scale)'

# dat: kx ky |F| なので、描画は ky:kx:mag
splot 'lorentzian_k_mag.dat' u 2:1:3 notitle
