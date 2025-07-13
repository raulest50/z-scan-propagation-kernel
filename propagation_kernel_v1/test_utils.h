// test_utils.h
#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <complex>
#include <cstdio>
#include "step_propagators.h"   // defines DIM
#include <hls_x_complex.h>     // defines hls::x_complex

// Print the top-left maxR×maxC block of a DIM×DIM complex matrix
// Works for std::complex<float> and hls::x_complex<float>
template<typename Cplx>
void printMatrix(const char* name,
                 Cplx M[DIM][DIM],
                 int maxR = 12, int maxC = 12);

// Read a DIM×DIM ASCII text "real imag" per line into std::complex<float>
void read_complex_matrix(const char* fn,
                         std::complex<float> M[DIM][DIM]);

// Write a DIM×DIM hls::x_complex<float> matrix out as ASCII text
void write_complex_matrix(const char* fn,
                          const hls::x_complex<float> M[DIM][DIM]);

#endif // TEST_UTILS_H