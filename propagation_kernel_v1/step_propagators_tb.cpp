#include "step_propagators.h"
#include <complex>
#include <cstdio>
#include <cmath>

static void read_complex_matrix(const char* fn, complex_t M[DIM][DIM]) {
    FILE* f = std::fopen(fn, "r");
    if (!f) {
        std::perror("fopen");
        std::exit(1);
    }
    for (int j = 0; j < DIM; ++j)
        for (int i = 0; i < DIM; ++i) {
            float re, im;
            std::fscanf(f, "%f %f", &re, &im);
            M[i][j] = complex_t{re, im};
        }
    std::fclose(f);
}

static float compare_matrices(const complex_t A[DIM][DIM], const complex_t B[DIM][DIM]) {
    double sum_sq = 0.0;
    for (int j = 0; j < DIM; ++j)
        for (int i = 0; i < DIM; ++i) {
            float dre = A[i][j].real() - B[i][j].real();
            float dim = A[i][j].imag() - B[i][j].imag();
            sum_sq += double(dre*dre + dim*dim);
        }
    return std::sqrt(sum_sq / double(DIM*DIM));
}

int main() {
    static complex_t phi_in[DIM][DIM];
    static complex_t phi_out[DIM][DIM];
    static complex_t phi_ref[DIM][DIM];

    read_complex_matrix("propagation_kernel_v1/validationData/full_step_within_tissue/initial_field.dat", phi_in);
    read_complex_matrix("propagation_kernel_v1/validationData/full_step_within_tissue/full_step_within_tissue_step_1.dat", phi_ref);

    propagation_step(phi_in, phi_out);

    float rms = compare_matrices(phi_out, phi_ref);
    std::printf("RMS error: %e\n", rms);
    return (rms < 1e-3f) ? 0 : 1;
}
