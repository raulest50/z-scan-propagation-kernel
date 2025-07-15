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

// Helper to read a 1-D complex vector from ASCII file.
static int read_complex_vector(const char* fn,
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
    return data.size();
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

// Compute RMS error between two vectors
static float rms_vector(const complex_t* A,
                        const complex_t* B,
                        int n)
{
    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        float dre = A[i].real() - B[i].real();
        float dim = A[i].imag() - B[i].imag();
        sum_sq += double(dre*dre + dim*dim);
    }
    return std::sqrt(sum_sq / double(n));
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

static void half_nonlin_ops_array(complex_t in[DIM], complex_t out[DIM]) {
    hls::stream<complex_t> s_in, s_out;
    for (int i = 0; i < DIM; ++i) {
        s_in.write(in[i]);
        half_nonlin_ops(s_in, s_out);
        out[i] = s_out.read();
    }
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

// Validate compute_b_vector using reference vectors
static void run_bvec_test(const char* name,
                          const char* x0_file,
                          const char* ref_file)
{
    std::vector<complex_t> buf_x0, buf_ref;
    int n = read_complex_vector(x0_file, buf_x0);
    int nref = read_complex_vector(ref_file, buf_ref);
    assert(n == nref);

    static complex_t X0[DIM], B[DIM], R[DIM];
    for (int i = 0; i < DIM; ++i) X0[i] = complex_t{0.f,0.f};
    for (int i = 0; i < n; ++i) X0[i] = buf_x0[i];
    for (int i = 0; i < n; ++i) R[i]  = buf_ref[i];
    for (int i = n; i < DIM; ++i) R[i] = complex_t{0.f,0.f};

    compute_b_vector(complex_t{2.0f,0.0f},
                     complex_t{2.0f,0.0f},
                     complex_t{2.0f,0.0f},
                     complex_t{-1.0f,0.0f},
                     X0, B);

    float rms = rms_vector(B, R, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3f) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

// Validate custom_thomas_solver using reference vectors
static void run_thomas_test(const char* name,
                            const char* b_file,
                            const char* ref_file)
{
    std::vector<complex_t> buf_b, buf_ref;
    int n = read_complex_vector(b_file, buf_b);
    int nref = read_complex_vector(ref_file, buf_ref);
    assert(n == nref);

    static complex_t Bvec[DIM], X[DIM], R[DIM];
    for (int i = 0; i < DIM; ++i) Bvec[i] = complex_t{0.f,0.f};
    for (int i = 0; i < DIM; ++i) X[i]    = complex_t{0.f,0.f};
    for (int i = 0; i < n; ++i) Bvec[i] = buf_b[i];
    for (int i = 0; i < n; ++i) R[i]    = buf_ref[i];
    for (int i = n; i < DIM; ++i) R[i] = complex_t{0.f,0.f};

    custom_thomas_solver(complex_t{2.0f,0.0f},
                         complex_t{2.0f,0.0f},
                         complex_t{2.0f,0.0f},
                         complex_t{-1.0f,0.0f},
                         Bvec, X);

    float rms = rms_vector(X, R, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3f) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

int main() {
    std::printf("Running step propagator tests...\n");

    run_adi_test("ADI X",
                 adi_x,
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\in.dat",
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\adi_x_out.dat");

    run_adi_test("ADI Y",
                 adi_y,
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\in.dat",
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\adi_y_out.dat");

    run_half_test("Half nonlin ops",
                  half_nonlin_ops_array,
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\in.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\half_2photon_absorption_out.dat");


    run_bvec_test("compute_b_vector",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\thomas\\bvec_x0.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\thomas\\bvec_result.dat");

    run_thomas_test("custom_thomas_solver",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\thomas\\thomas_b.dat",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\thomas\\thomas_x_expected.dat");

    run_full_step("Propagation step 1",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\full_step_within_tissue_step_1.dat",
                  1);

    run_full_step("Propagation step 40",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\full_step_within_tissue_step_40.dat",
                  40);

    run_full_step("Propagation step 80",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\full_step_within_tissue_step_80.dat",
                  80);

    run_full_step("Propagation step 120",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validationData\\full_step_within_tissue\\full_step_within_tissue_step_120.dat",
                  120);

    return 0;
}
