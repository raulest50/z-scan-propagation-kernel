// test_utils.cpp
#include "test_utils.h"
#include <fstream>
#include <algorithm>
#include <cmath>

// -----------------------------------------------------------------------------
// printMatrix implementation
// -----------------------------------------------------------------------------
template<typename Cplx>
void printMatrix(const char* name,
                 Cplx M[DIM][DIM],
                 int maxR, int maxC)
{
    int rows = std::min(maxR, DIM);
    int cols = std::min(maxC, DIM);
    std::fprintf(stdout, "\n=== %s [%dx%d] ===\n", name, rows, cols);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            auto &z = M[c][r];
            std::fprintf(stdout, "(%8.5f,%8.5f) ",
                         float(z.real()), float(z.imag()));
        }
        std::fprintf(stdout, "\n");
    }
    std::fflush(stdout);
}

// Explicit instantiations
template void printMatrix<std::complex<float>>(const char*, std::complex<float>[DIM][DIM], int, int);
#ifndef HLS_STUB_H
template void printMatrix<hls::x_complex<float>>(const char*, hls::x_complex<float>[DIM][DIM], int, int);
#endif

// -----------------------------------------------------------------------------
// read_complex_matrix
// -----------------------------------------------------------------------------
void read_complex_matrix(const char* fn,
                         std::complex<float> M[DIM][DIM])
{
    std::ifstream in(fn);
    if (!in) {
        std::fprintf(stderr, "Error: cannot open %s\n", fn);
        std::exit(1);
    }
    for (int j = 0; j < DIM; ++j)
        for (int i = 0; i < DIM; ++i) {
            float re, im;
            in >> re >> im;
            M[i][j] = {re, im};
        }
}

// -----------------------------------------------------------------------------
// write_complex_matrix
// -----------------------------------------------------------------------------
void write_complex_matrix(const char* fn,
                          const hls::x_complex<float> M[DIM][DIM])
{
    std::ofstream out(fn);
    if (!out) {
        std::fprintf(stderr, "Error: cannot write %s\n", fn);
        std::exit(1);
    }
    for (int j = 0; j < DIM; ++j)
        for (int i = 0; i < DIM; ++i)
            out << M[i][j].real() << " " << M[i][j].imag() << "\n";
}
