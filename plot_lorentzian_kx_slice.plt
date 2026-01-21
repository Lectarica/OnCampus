# plot_lorentzian_kx_slice.plt
reset

infile = "lorentzian_k_fft_exact.dat"

set term pngcairo size 1000,700
set output "lorentzian_kx_slice_ky0.png"

set xlabel "k_x"
set ylabel "F(k_x, k_y=0)"
set grid
set key left top
set tics nomirror

# --- IMPORTANT ---
# ここでは「ky が 4列目」と仮定して ky==0 の行だけ抽出しています。
# 数値: |F_fft| = 5列目, 厳密: F_exact = 8列目 という仮定です。
#
# もし列が違うなら、$4, $5, $8 の数字だけ変えればOKです。

fft_curve   = sprintf("< awk '($4==0){print $3, $5}' %s", infile)
exact_curve = sprintf("< awk '($4==0){print $3, $8}' %s", infile)

plot fft_curve   with lines title "|F_{FFT}(k_x,0)|", \
     exact_curve with lines title "Exact"
