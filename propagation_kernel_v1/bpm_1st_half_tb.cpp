// bpm_1st_half_tb.cpp
#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include "bpm_split_step.h"  // or wherever bpm_1st_half is declared

// Dimensions and parameters (must match your kernel)
const int NDX = 46;
const int NDY = 46;
const float DX = 1.5e-05;
const float DY = 1.5e-05;
const float DZ = 0.0001;
const float k   = 7853981.6339;
const float n0  = 1.0;
extern float n2[DIM][DIM];  // adjust linkage as needed

// Utility: read a complex matrix from a file (real imag pairs)
void read_complex_matrix(const char* filename,
                         std::complex<float> M[DIM][DIM]) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        std::exit(1);
    }
    for (int j = 0; j <= NDY; ++j) {
        for (int i = 0; i <= NDX; ++i) {
            float re, im;
            in >> re >> im;
            M[i][j] = std::complex<float>(re, im);
        }
    }
    in.close();
}

// Utility: write a complex matrix to a file (real imag pairs)
void write_complex_matrix(const char* filename,
                          std::complex<float> M[DIM][DIM]) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: cannot write " << filename << "\n";
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
    // Allocate arrays
    static std::complex<float> PHI_m[DIM][DIM];
    static std::complex<float> PHI_auxNL[DIM][DIM];
    static std::complex<float> PHI_half[DIM][DIM];

    // 1) Read inputs
    read_complex_matrix("phi_m0.dat",  PHI_m);
    read_complex_matrix("phi_aux.dat", PHI_auxNL);

    // 2) Call the HLS kernel
    bpm_1st_half(
        // convert std::complex<> to hls::x_complex<> if needed
        (hls::x_complex<float> (*)[DIM])PHI_m,
        (hls::x_complex<float> (*)[DIM])PHI_auxNL,
        k, n0, n2, NDX, NDY, DX, DY, DZ,
        (hls::x_complex<float> (*)[DIM])PHI_half
    );

    // 3) Write CSIM output
    write_complex_matrix("phi_half_csim.dat", PHI_half);

    // 4) Self-check: compare to golden reference
    int diff_ret = std::system(
        "diff --brief -w phi_half_csim.dat phi_half_ref.dat"
    );
    if (diff_ret == 0) {
        std::cout << "[PASS] Output matches phi_half_ref.dat\n";
        return 0;
    } else {
        std::cerr << "[FAIL] Discrepancy detected against phi_half_ref.dat\n";
        return 1;
    }
}
