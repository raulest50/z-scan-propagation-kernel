// step_propagators_tb.cpp

#include "step_propagators.h"
#include <complex>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <cassert>

// Helper to read an arbitrary sized matrix from ASCII file.
static int read_complex_matrix_dyn(const char* fn,
                                   std::vector<complex_t>& data)
{
    std::ifstream in(fn);
    if (!in) {
        std::perror("fopen");
        std::exit(1);
    }
    float re, im;
    while (in >> re >> im)
        data.emplace_back(re, im);
    int n = std::sqrt(double(data.size()));
    assert(n*n == (int)data.size());
    return n;
}

// Copy data into a DIMxDIM matrix (zero padded)
static void load_into(const std::vector<complex_t>& src, int n,
                      complex_t dst[DIM][DIM])
{
    for (int j = 0; j < DIM; ++j)
        for (int i = 0; i < DIM; ++i)
            dst[i][j] = complex_t{0.f,0.f};
    for (int j = 0, idx=0; j < n; ++j)
        for (int i = 0; i < n; ++i, ++idx)
            dst[i][j] = src[idx];
}

// Compute RMS error on the top-left NxN block
static float rms_region(const complex_t A[DIM][DIM],
                        const complex_t B[DIM][DIM], int n)
{
    double sum_sq = 0.0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            float dre = A[i][j].real() - B[i][j].real();
            float dim = A[i][j].imag() - B[i][j].imag();
            sum_sq += double(dre*dre + dim*dim);
        }
    return std::sqrt(sum_sq / double(n*n));
}

static void run_adi_test(const char* name,
                         void (*func)(complex_t[DIM][DIM], complex_t[DIM][DIM]),
                         const char* in_file,
                         const char* ref_file)
{
    std::vector<complex_t> buf_in, buf_ref;
    int n = read_complex_matrix_dyn(in_file, buf_in);
    int nref = read_complex_matrix_dyn(ref_file, buf_ref);
    assert(n == nref);

    static complex_t A[DIM][DIM], B[DIM][DIM], R[DIM][DIM];
    load_into(buf_in, n, A);
    load_into(buf_ref, n, R);

    func(A, B);

    float rms = rms_region(B, R, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3f) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

static void run_half_test(const char* name,
                          void (*func)(complex_t[DIM], complex_t[DIM]),
                          const char* in_file,
                          const char* ref_file)
{
    std::vector<complex_t> buf_in, buf_ref;
    int n = read_complex_matrix_dyn(in_file, buf_in);
    int nref = read_complex_matrix_dyn(ref_file, buf_ref);
    assert(n == nref);

    static complex_t A[DIM], B[DIM], R[DIM];
    for (int i = 0; i < DIM; ++i) A[i] = complex_t{0.f,0.f};
    for (int i = 0; i < n; ++i) A[i] = buf_in[i];
    for (int i = 0; i < n; ++i) R[i] = buf_ref[i];
    for (int i = n; i < DIM; ++i) R[i] = complex_t{0.f,0.f};

    func(A, B);

    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        float dre = B[i].real() - R[i].real();
        float dim = B[i].imag() - R[i].imag();
        sum_sq += double(dre*dre + dim*dim);
    }
    float rms = std::sqrt(sum_sq / double(n));
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3f) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

static void run_full_step(const char* name,
                          const char* init_file,
                          const char* ref_file,
                          int steps)
{
    std::vector<complex_t> buf_in, buf_ref;
    int n = read_complex_matrix_dyn(init_file, buf_in);
    int nref = read_complex_matrix_dyn(ref_file, buf_ref);
    assert(n == nref);

    static complex_t A[DIM][DIM], B[DIM][DIM], tmp[DIM][DIM];
    load_into(buf_in, n, A);
    load_into(buf_ref, n, tmp);

    for (int s = 0; s < steps; ++s) {
        propagation_step(A, B);
        for (int j = 0; j < DIM; ++j)
            for (int i = 0; i < DIM; ++i)
                A[i][j] = B[i][j];
    }

    float rms = rms_region(A, tmp, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3f) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

int main() {
    std::printf("Running step propagator tests...\n");

    run_adi_test("ADI X",
                 adi_x,
                 "propagation_kernel_v1/validationData/in.dat",
                 "propagation_kernel_v1/validationData/adi_x_out.dat");

    run_adi_test("ADI Y",
                 adi_y,
                 "propagation_kernel_v1/validationData/in.dat",
                 "propagation_kernel_v1/validationData/adi_y_out.dat");

    run_half_test("Half nonlinear",
                  half_nonlinear,
                  "propagation_kernel_v1/validationData/in.dat",
                  "propagation_kernel_v1/validationData/half_nonlinear_out.dat");

    run_half_test("Half linear absorption",
                  half_linear_absorption,
                  "propagation_kernel_v1/validationData/in.dat",
                  "propagation_kernel_v1/validationData/half_linear_absorption_out.dat");

    run_half_test("Half 2-photon absorption",
                  half_2photon_absorption,
                  "propagation_kernel_v1/validationData/in.dat",
                  "propagation_kernel_v1/validationData/half_2photon_absorption_out.dat");

    run_full_step("Propagation step 1",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/initial_field.dat",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/full_step_within_tissue_step_1.dat",
                  1);

    run_full_step("Propagation step 40",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/initial_field.dat",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/full_step_within_tissue_step_40.dat",
                  40);

    run_full_step("Propagation step 80",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/initial_field.dat",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/full_step_within_tissue_step_80.dat",
                  80);

    run_full_step("Propagation step 120",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/initial_field.dat",
                  "propagation_kernel_v1/validationData/full_step_within_tissue/full_step_within_tissue_step_120.dat",
                  120);

    return 0;
}

