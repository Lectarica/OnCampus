set term pngcairo size 900,800
set output 'lorentzian_real.png'

set size ratio -1
set xlabel 'y'
set ylabel 'x'
set title 'Real-space Lorentzian f(x,y)'

set pm3d map
set cblabel 'f'

# dat: x y f なので、描画は y:x:val にする（見た目を x縦, y横に）
splot 'lorentzian_real.dat' u 2:1:3 notitle
