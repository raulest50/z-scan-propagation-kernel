// bpm_1st_half_tb.cpp
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include "bpm_split_step.h"  // defines DIM and bpm_1st_half, includes <hls_x_complex.h>

const int NDX = 46;
const int NDY = 46;
const float DX = 1.5e-05f;
const float DY = 1.5e-05f;
const float DZ = 1e-4f;
const float k   = 7853981.6339f;
const float n0  = 1.0f;
const float n2  = 2.5e-20f;

// Read a DIM×DIM complex matrix from disk into a std::complex array
void read_complex_matrix(const char* fn, std::complex<float> M[DIM][DIM]) {
    std::ifstream in(fn);
    if (!in) {
        std::cerr << "Error: cannot open " << fn << std::endl;
        std::exit(1);
    }
    for (int j = 0; j <= NDY; ++j) {
        for (int i = 0; i <= NDX; ++i) {
            float re, im;
            in >> re >> im;
            M[i][j] = {re, im};
        }
    }
    in.close();
}

// Write a DIM×DIM hls::x_complex matrix to disk
void write_complex_matrix(const char* fn, const hls::x_complex<float> M[DIM][DIM]) {
    std::ofstream out(fn);
    if (!out) {
        std::cerr << "Error: cannot write " << fn << std::endl;
        std::exit(1);
    }
    for (int j = 0; j <= NDY; ++j) {
        for (int i = 0; i <= NDX; ++i) {
            out << M[i][j].real() << " "
                << M[i][j].imag() << "\n";
        }
    }
    out.close();
}

int main() {
    // 1) Allocate host‐side std::complex buffers and HLS x_complex buffers
    static std::complex<float>    PHI_m_std[DIM][DIM];
    static std::complex<float>    PHI_auxNL_std[DIM][DIM];
    static hls::x_complex<float>  PHI_m[DIM][DIM];
    static hls::x_complex<float>  PHI_auxNL[DIM][DIM];
    static hls::x_complex<float>  PHI_half[DIM][DIM];

    // 2) Read inputs from absolute paths
    read_complex_matrix(
      "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_m0.dat",
      PHI_m_std);
    read_complex_matrix(
      "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_aux.dat",
      PHI_auxNL_std);

    // 3) Debug print first 3×3 of PHI_m_std
    std::cout << "PHI_m_std[0..2][0..2]:\n";
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            auto &c = PHI_m_std[i][j];
            std::cout << "(" << c.real() << "," << c.imag() << ") ";
        }
        std::cout << "\n";
    }
    std::cout << std::flush;

    // 4) Convert to HLS types
    for (int j = 0; j <= NDY; ++j) {
        for (int i = 0; i <= NDX; ++i) {
            PHI_m[i][j]     = hls::x_complex<float>(PHI_m_std[i][j].real(),
                                                     PHI_m_std[i][j].imag());
            PHI_auxNL[i][j] = hls::x_complex<float>(PHI_auxNL_std[i][j].real(),
                                                     PHI_auxNL_std[i][j].imag());
        }
    }

    // 5) Invoke the kernel
    bpm_1st_half(
      PHI_m, PHI_auxNL,
      k, n0, n2,
      NDX, NDY, DX, DY, DZ,
      PHI_half
    );

    // 6) Debug print first 3×3 of PHI_half
    std::cout << "PHI_half[0..2][0..2]:\n";
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            auto &c = PHI_half[i][j];
            std::cout << "(" << c.real() << "," << c.imag() << ") ";
        }
        std::cout << "\n";
    }
    std::cout << std::flush;

    // 7) Write out the result
    write_complex_matrix(
      "C:/Vws/z_scan_acceleration_ovr/propagation_kernel_v1/bpm_1st_h_testfiles/phi_half_csim.dat",
      PHI_half
    );

    return 0;  // report success so CSIM passes
}
