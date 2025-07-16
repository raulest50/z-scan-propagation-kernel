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

    // Modificado para leer valores como ap_fixed en lugar de float
    data_t re, im;
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

    // Modificado para leer valores como ap_fixed en lugar de float
    data_t re, im;
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
            dst[i][j] = complex_t{0, 0};
    for (int j = 0, idx=0; j < n; ++j)
        for (int i = 0; i < n; ++i, ++idx)
            dst[i][j] = src[idx];
}

// Compute RMS error on the top-left NxN block
static double rms_region(const complex_t A[DIM][DIM],
                        const complex_t B[DIM][DIM], int n)
{
    double sum_sq = 0.0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            double dre = A[i][j].real() - B[i][j].real();
            double dim = A[i][j].imag() - B[i][j].imag();
            sum_sq += dre*dre + dim*dim;
        }
    return std::sqrt(sum_sq / double(n*n));
}

// Compute RMS error between two vectors
static double rms_vector(const complex_t* A,
                        const complex_t* B,
                        int n)
{
    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double dre = A[i].real() - B[i].real();
        double dim = A[i].imag() - B[i].imag();
        sum_sq += dre*dre + dim*dim;
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

    double rms = rms_region(B, R, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3) std::printf("  [FAILED]\n");
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
    for (int i = 0; i < DIM; ++i) A[i] = complex_t{0, 0};
    for (int i = 0; i < n; ++i) A[i] = buf_in[i];
    for (int i = 0; i < n; ++i) R[i] = buf_ref[i];
    for (int i = n; i < DIM; ++i) R[i] = complex_t{0, 0};

    func(A, B);

    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double dre = B[i].real() - R[i].real();
        double dim = B[i].imag() - R[i].imag();
        sum_sq += dre*dre + dim*dim;
    }
    double rms = std::sqrt(sum_sq / double(n));
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3) std::printf("  [FAILED]\n");
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

    double rms = rms_region(A, tmp, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

// Validate compute_b_vector using reference vectors
static void run_bvec_test(const char* name,
                          const char* dp_file,
                          const char* dp1_file,
                          const char* dp2_file,
                          const char* do_file,
                          const char* x0_file,
                          const char* ref_file)
{
    std::vector<complex_t> buf_dp, buf_dp1, buf_dp2, buf_do, buf_x0, buf_ref;

    // Leer parámetros individuales
    read_complex_vector(dp_file, buf_dp);
    read_complex_vector(dp1_file, buf_dp1);
    read_complex_vector(dp2_file, buf_dp2);
    read_complex_vector(do_file, buf_do);

    int n = read_complex_vector(x0_file, buf_x0);
    int nref = read_complex_vector(ref_file, buf_ref);
    assert(n == nref);

    static complex_t X0[DIM], B[DIM], R[DIM];
    for (int i = 0; i < DIM; ++i) X0[i] = complex_t{0, 0};
    for (int i = 0; i < n; ++i) X0[i] = buf_x0[i];
    for (int i = 0; i < n; ++i) R[i]  = buf_ref[i];
    for (int i = n; i < DIM; ++i) R[i] = complex_t{0, 0};

    compute_b_vector(buf_dp[0],
                     buf_dp1[0],
                     buf_dp2[0],
                     buf_do[0],
                     X0, B);

    double rms = rms_vector(B, R, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

// Validate custom_thomas_solver using reference vectors
static void run_thomas_test(const char* name,
                            const char* dp_file,
                            const char* dp1_file,
                            const char* dp2_file,
                            const char* do_file,
                            const char* b_file,
                            const char* ref_file)
{
    std::vector<complex_t> buf_dp, buf_dp1, buf_dp2, buf_do, buf_b, buf_ref;

    // Leer parámetros individuales
    read_complex_vector(dp_file, buf_dp);
    read_complex_vector(dp1_file, buf_dp1);
    read_complex_vector(dp2_file, buf_dp2);
    read_complex_vector(do_file, buf_do);

    int n = read_complex_vector(b_file, buf_b);
    int nref = read_complex_vector(ref_file, buf_ref);
    assert(n == nref);

    static complex_t Bvec[DIM], X[DIM], R[DIM];
    for (int i = 0; i < DIM; ++i) Bvec[i] = complex_t{0, 0};
    for (int i = 0; i < DIM; ++i) X[i]    = complex_t{0, 0};
    for (int i = 0; i < n; ++i) Bvec[i] = buf_b[i];
    for (int i = 0; i < n; ++i) R[i]    = buf_ref[i];
    for (int i = n; i < DIM; ++i) R[i] = complex_t{0, 0};

    custom_thomas_solver(buf_dp[0],
                         buf_dp1[0],
                         buf_dp2[0],
                         buf_do[0],
                         Bvec, X);

    double rms = rms_vector(X, R, n);
    std::printf("%s RMS error: %e\n", name, rms);
    if (rms >= 1e-3) std::printf("  [FAILED]\n");
    else             std::printf("  [OK]\n");
}

int main() {
    std::printf("Running step propagator tests with ap_fixed data type...\n");

    // Actualizado para usar los archivos de validation_data_main
    run_adi_test("ADI X",
                 adi_x,
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\individual_ops_in.dat",
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\adi_x_out.dat");

    run_adi_test("ADI Y",
                 adi_y,
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\individual_ops_in.dat",
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\adi_y_out.dat");

    // Prueba para operadores no lineales combinados
    run_adi_test("Half Nonlinear Ops Combined",
                 [](complex_t in[DIM][DIM], complex_t out[DIM][DIM]) {
                     // Implementación para aplicar los operadores no lineales combinados
                     // Esta es una simplificación, la implementación real dependerá de cómo
                     // se aplican estos operadores en el código
                     for (int j = 0; j < DIM; j++) {
                         for (int i = 0; i < DIM; i++) {
                             hls::stream<complex_t> s_in, s_out;
                             s_in.write(in[i][j]);
                             half_nonlin_ops(s_in, s_out);
                             out[i][j] = s_out.read();
                         }
                     }
                 },
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\individual_ops_in.dat",
                 "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\half_nonlinear_ops_combined_out.dat");

    // Prueba para compute_b_vector con los nuevos archivos de parámetros
    run_bvec_test("compute_b_vector",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_dp_in.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_dp1_in.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_dp2_in.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_do_in.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_x0_in.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_out.dat");

    // Prueba para custom_thomas_solver con los nuevos archivos de parámetros
    run_thomas_test("custom_thomas_solver",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_dp_in.dat",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_dp1_in.dat",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_dp2_in.dat",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_do_in.dat",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\b_vector_out.dat",
                    "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\thomas_solver_out.dat");

    // Pruebas para propagación completa
    run_full_step("Propagation step 1",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\full_step_within_tissue_step_1.dat",
                  1);

    run_full_step("Propagation step 40",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\full_step_within_tissue_step_40.dat",
                  40);

    run_full_step("Propagation step 80",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\full_step_within_tissue_step_80.dat",
                  80);

    run_full_step("Propagation step 120",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\initial_field.dat",
                  "C:\\Vws\\z_scan_acceleration_ovr\\propagation_kernel_v1\\validation_data_main\\full_step_within_tissue_step_120.dat",
                  120);

    return 0;
}
