#include <fftw3.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

static inline int idx(int ix, int iy, int Ny) { return ix * Ny + iy; }

int main(int argc, char** argv) {
    int Nx = 128;
    int Ny = 128;
    double gamma = 6.0;

    if (argc >= 3) {
        Nx = std::atoi(argv[1]);
        Ny = std::atoi(argv[2]);
    }
    if (argc >= 4) {
        gamma = std::atof(argv[3]);
    }

    fftw_complex* in  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

    // 実空間でローレンチアンを作る
    for (int ix = 0; ix < Nx; ++ix) {
        const double x = ix - Nx / 2.0;
        for (int iy = 0; iy < Ny; ++iy) {
            const double y = iy - Ny / 2.0;
            const double val = (gamma * gamma) / (x * x + y * y + gamma * gamma);
            in[idx(ix, iy, Ny)][0] = val;  // 実部
            in[idx(ix, iy, Ny)][1] = 0.0;  // 虚部
        }
    }

    fftw_plan plan = fftw_plan_dft_2d(
        Nx, Ny,
        in, out,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );
    if (!plan) {
        std::cerr << "fftw_plan_dft_2d failed\n";
        fftw_free(in);
        fftw_free(out);
        return 1;
    }

    // FFTの実行
    fftw_execute(plan);

    // ファイル出力
    {
        std::ofstream ofs("lorentzian_real.dat");
        ofs << "# x y f(x,y)\n";
        for (int ix = 0; ix < Nx; ++ix) {
            const int x = ix - Nx / 2;
            for (int iy = 0; iy < Ny; ++iy) {
                const int y = iy - Ny / 2;
                ofs << x << " " << y << " " << in[idx(ix, iy, Ny)][0] << "\n";
            }
            ofs << "\n";
        }
    }

    // k空間のファイル出力
    {
        std::ofstream ofs("lorentzian_k_mag.dat");
        ofs << "# kx ky |F(k)| (fftshifted indices)\n";

        for (int sx = 0; sx < Nx; ++sx) {
            const int kx = sx - Nx / 2;
            const int ix = (sx + Nx / 2) % Nx;  // shifted -> FFTW index

            for (int sy = 0; sy < Ny; ++sy) {
                const int ky = sy - Ny / 2;
                const int iy = (sy + Ny / 2) % Ny;

                const double re = out[idx(ix, iy, Ny)][0];
                const double im = out[idx(ix, iy, Ny)][1];
                const double mag = std::sqrt(re * re + im * im);

                ofs << kx << " " << ky << " " << mag << "\n";
            }
            ofs << "\n";
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    std::cout << "Done.\n"
              << "Outputs: lorentzian_real.dat, lorentzian_k_mag.dat\n";
    return 0;
}
