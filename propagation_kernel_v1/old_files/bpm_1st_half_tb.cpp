// bpm_1st_half_tb.cpp
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>                    // for fprintf
#include "bpm_split_step.h"         // defines DIM and bpm_1st_half, includes <hls_x_complex.h>
#include "test_utils.h"

const int NDX = 46;
const int NDY = 46;
const float DX = 1.5e-05f;
const float DY = 1.5e-05f;
const float DZ = 1e-4f;
const float k   = 7853981.6339f;
const float n0  = 1.0f;
const float n2  = 2.5e-20f;


int bpm_1sth_only_test() {
    std::fprintf(stdout, "SIMULATION START *********************   \n");
    std::fprintf(stdout, " -_-\n");
    // Host buffers
    static std::complex<float>    PHI_m_std[DIM][DIM];
    static std::complex<float>    PHI_auxNL_std[DIM][DIM];
    // HLS buffers
    static hls::x_complex<float>  PHI_m[DIM][DIM];
    static hls::x_complex<float>  PHI_auxNL[DIM][DIM];
    static hls::x_complex<float>  PHI_half[DIM][DIM];

    // 1) Read inputs
    read_complex_matrix(
      "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_m0.dat",
      PHI_m_std);
    read_complex_matrix(
      "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_aux.dat",
      PHI_auxNL_std);

    // 2) Print host‐side inputs
    printMatrix("PHI_m_std",     PHI_m_std);
    printMatrix("PHI_auxNL_std", PHI_auxNL_std);

    // 3) Convert to HLS types
    for (int j = 0; j <= NDY; ++j)
      for (int i = 0; i <= NDX; ++i) {
        PHI_m    [i][j] = hls::x_complex<float>(
                            PHI_m_std[i][j].real(),
                            PHI_m_std[i][j].imag());
        PHI_auxNL[i][j] = hls::x_complex<float>(
                            PHI_auxNL_std[i][j].real(),
                            PHI_auxNL_std[i][j].imag());
      }

    // 4) Call kernel
    bpm_1st_half(PHI_m, PHI_auxNL, k, n0, n2, NDX, NDY, DX, DY, DZ, PHI_half);

    // 5) Print HLS output
    printMatrix("PHI_half", PHI_half);

    // 6) Write simulation output
    const char* out_fn =
      "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_half_csim.dat";
    write_complex_matrix(out_fn, PHI_half);

    // 7) Read Python reference and print & metric
    static std::complex<float> PHI_ref[DIM][DIM];
    {
      std::ifstream in(
        "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_half_ref.dat"
      );
      if (!in) {
        std::fprintf(stderr, "Error: cannot open reference file\n");
        return 1;
      }
      for (int j = 0; j <= NDY; ++j)
        for (int i = 0; i <= NDX; ++i) {
          float re, im;
          in >> re >> im;
          PHI_ref[i][j] = {re, im};
        }
    }
    printMatrix("PHI_ref", PHI_ref);

    // 8) Compute and print error metrics
    double sum_sq = 0.0, max_err = 0.0;
    for (int j = 0; j <= NDY; ++j) {
      for (int i = 0; i <= NDX; ++i) {
        double dre = double(PHI_half[i][j].real()) - PHI_ref[i][j].real();
        double dim = double(PHI_half[i][j].imag()) - PHI_ref[i][j].imag();
        double err = std::hypot(dre, dim);
        sum_sq += err*err;
        if (err > max_err) max_err = err;
      }
    }
    double rms = std::sqrt(sum_sq / double((NDX+1)*(NDY+1)));
    std::fprintf(stdout,
      "\nMax abs error = %.6e, RMS error = %.6e\n",
      max_err, rms
    );
    fflush(stdout);

    return 0;
}
