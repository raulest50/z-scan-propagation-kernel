// step_propagators.cpp
#include "step_propagators.h"
#if __has_include(<hls_math.h>)
#  include <hls_math.h>
#else
#  include "hls_stub.h"
#endif

const float eps         = 1e-12f;             // Epsilon para estabilidad
const float k           = 10681416.0f;        // Número de onda
const float dz          = 1e-6f;              // Paso en z
const float dy          = 1.7578e-7f;         // Paso en y
const float dx          = 1.7578e-7f;         // Paso en x
const float n0          = 1.36f;              // Índice de refracción lineal
const float n2          = 3e-20f;             // Índice de refracción no lineal
const float alpha       = 0.3f;               // Coeficiente de absorción lineal
const float beta        = 1e-11f;             // Coeficiente de absorción de dos fotones

// Factor para método ADI en X/Y
complex_t ung = complex_t{0.0f, dz/(4.0f * k * dy * dy)};

// Constantes para efectos no lineal y absorción
const float phase_const   = k * n2 * dz / 2.0f;
const float attenuation   = hls::exp(-alpha * dz / 4.0f);
const float tpa_const     = -beta  * dz / 4.0f;


/**
 * Resuelve un sistema tridiagonal con estructura especial (Thomas).
 */
void custom_thomas_solver(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t do_val,
    complex_t b[DIM],
    complex_t x[DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=dp,dp1,dp2,do_val bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=b  offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=x  offset=slave bundle=GMEM1

    complex_t c_prime[DIM-1];
    complex_t d_prime[DIM];
    #pragma HLS ARRAY_PARTITION variable=c_prime complete
    #pragma HLS ARRAY_PARTITION variable=d_prime complete

    // Forward elimination
    c_prime[0]   = do_val / dp1;
    d_prime[0]   = b[0]   / dp1;
    for (int i = 1; i < DIM-1; i++) {
        #pragma HLS PIPELINE II=1
        complex_t denom = dp - do_val * c_prime[i-1];
        c_prime[i]      = do_val / denom;
        d_prime[i]      = (b[i] - do_val * d_prime[i-1]) / denom;
    }
    d_prime[DIM-1] = (b[DIM-1] - do_val * d_prime[DIM-2]) /
                     (dp2 - do_val * c_prime[DIM-2]);

    // Back substitution
    x[DIM-1] = d_prime[DIM-1];
    for (int i = DIM-2; i >= 0; i--) {
        #pragma HLS PIPELINE II=1
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}

/**
 * Multiplica una matriz tridiagonal especial por un vector.
 */
void compute_b_vector(
    complex_t dp,
    complex_t dp1,
    complex_t dp2,
    complex_t do_val,
    complex_t x0[DIM],
    complex_t b[DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE s_axilite port=dp,dp1,dp2,do_val bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=x0 offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=b   offset=slave bundle=GMEM1

    b[0] = dp1 * x0[0] + do_val * x0[1];
    for (int i = 1; i < DIM-1; i++) {
        #pragma HLS PIPELINE II=1
        b[i] = do_val * x0[i-1] + dp * x0[i] + do_val * x0[i+1];
    }
    b[DIM-1] = do_val * x0[DIM-2] + dp2 * x0[DIM-1];
}

/**
 * ADI en dirección X.
 */
void adi_x(
    complex_t phi[DIM][DIM],
    complex_t phi_inter[DIM][DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=phi       offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=phi_inter offset=slave bundle=GMEM1

    complex_t phi_col[DIM], b_vec[DIM], result[DIM];
    #pragma HLS ARRAY_PARTITION variable=phi_col complete
    #pragma HLS ARRAY_PARTITION variable=b_vec   complete
    #pragma HLS ARRAY_PARTITION variable=result  complete

    for (int col = 0; col < DIM; col++) {
        #pragma HLS DATAFLOW

        // Carga de columna
        for (int row = 0; row < DIM; row++) {
            #pragma HLS PIPELINE II=1
            phi_col[row] = phi[row][col];
        }

        // Ratio frontera
        complex_t ratio0, ratioN;
        if (hls::hypot(phi_col[1].real(), phi_col[1].imag()) < eps) {
            ratio0 = complex_t{1.0f, 0.0f};
        } else {
            ratio0 = phi_col[0] / phi_col[1];
        }
        if (hls::hypot(phi_col[DIM-2].real(), phi_col[DIM-2].imag()) < eps) {
            ratioN = complex_t{1.0f, 0.0f};
        } else {
            ratioN = phi_col[DIM-1] / phi_col[DIM-2];
        }

        // Matriz B
        complex_t dp1_B = (-ung - ung) + complex_t{1.0f,0.0f} + ung * ratio0;
        complex_t dp2_B = (-ung - ung) + complex_t{1.0f,0.0f} + ung * ratioN;
        complex_t dp_B  = (-ung - ung) + complex_t{1.0f,0.0f};
        complex_t do_B  =  ung;
        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_col, b_vec);

        // Matriz A
        complex_t dp1_A = ung + ung + complex_t{1.0f,0.0f} - ung * ratio0;
        complex_t dp2_A = ung + ung + complex_t{1.0f,0.0f} - ung * ratioN;
        complex_t dp_A  = ung + ung + complex_t{1.0f,0.0f};
        complex_t do_A  = -ung;
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b_vec, result);

        // Escritura de resultados
        for (int row = 0; row < DIM; row++) {
            #pragma HLS PIPELINE II=1
            phi_inter[row][col] = result[row];
        }
    }
}

/**
 * ADI en dirección Y.
 */
void adi_y(
    complex_t phi[DIM][DIM],
    complex_t phi_inter[DIM][DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=phi       offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=phi_inter offset=slave bundle=GMEM1

    complex_t phi_row[DIM], b_vec[DIM], result[DIM];
    #pragma HLS ARRAY_PARTITION variable=phi_row complete
    #pragma HLS ARRAY_PARTITION variable=b_vec    complete
    #pragma HLS ARRAY_PARTITION variable=result   complete

    for (int row = 0; row < DIM; row++) {
        #pragma HLS DATAFLOW

        // Carga de fila
        for (int col = 0; col < DIM; col++) {
            #pragma HLS PIPELINE II=1
            phi_row[col] = phi[row][col];
        }

        // Ratio frontera
        complex_t ratio0, ratioN;
        if (hls::hypot(phi_row[1].real(), phi_row[1].imag()) < eps) {
            ratio0 = complex_t{1.0f, 0.0f};
        } else {
            ratio0 = phi_row[0] / phi_row[1];
        }
        if (hls::hypot(phi_row[DIM-2].real(), phi_row[DIM-2].imag()) < eps) {
            ratioN = complex_t{1.0f,0.0f};
        } else {
            ratioN = phi_row[DIM-1] / phi_row[DIM-2];
        }

        // Matriz B
        complex_t dp1_B = -ung - ung + complex_t{1.0f,0.0f} + ung * ratio0;
        complex_t dp2_B = -ung - ung + complex_t{1.0f,0.0f} + ung * ratioN;
        complex_t dp_B  = -ung - ung + complex_t{1.0f,0.0f};
        complex_t do_B  =  ung;
        compute_b_vector(dp_B, dp1_B, dp2_B, do_B, phi_row, b_vec);

        // Matriz A
        complex_t dp1_A = ung + ung + complex_t{1.0f,0.0f} - ung * ratio0;
        complex_t dp2_A = ung + ung + complex_t{1.0f,0.0f} - ung * ratioN;
        complex_t dp_A  = ung + ung + complex_t{1.0f,0.0f};
        complex_t do_A  = -ung;
        custom_thomas_solver(dp_A, dp1_A, dp2_A, do_A, b_vec, result);

        // Escritura de resultados
        for (int col = 0; col < DIM; col++) {
            #pragma HLS PIPELINE II=1
            phi_inter[row][col] = result[col];
        }
    }
}

/**
 * Kerr non-linear effect (half-step).
 */
void half_nonlinear(
    complex_t phi[DIM],
    complex_t phi_out[DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=phi     offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=GMEM1

    for (int i = 0; i < DIM; i++) {
        #pragma HLS PIPELINE II=1
        float abs_sq = phi[i].real()*phi[i].real() + phi[i].imag()*phi[i].imag();
        float angle  = phase_const * abs_sq;
        float sinv, cosv;
        hls::sincos(angle, &sinv, &cosv);
        phi_out[i] = phi[i] * complex_t{cosv, sinv};
    }
}

/**
 * Linear absorption (half-step).
 */
void half_linear_absorption(
    complex_t phi[DIM],
    complex_t phi_out[DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=phi     offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=GMEM1

    for (int i = 0; i < DIM; i++) {
        #pragma HLS PIPELINE II=1
        phi_out[i] = phi[i] * attenuation;
    }
}

/**
 * Two-photon absorption (half-step).
 */
void half_2photon_absorption(
    complex_t phi[DIM],
    complex_t phi_out[DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=phi     offset=slave bundle=GMEM0
    #pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=GMEM1

    for (int i = 0; i < DIM; i++) {
        #pragma HLS PIPELINE II=1
        float abs_sq = phi[i].real()*phi[i].real() + phi[i].imag()*phi[i].imag();
        float att    = hls::exp(tpa_const * abs_sq);
        phi_out[i]   = phi[i] * att;
    }
}

/**
 * Full propagation step: ADI X → half-TPA → half-Kerr → half-linear
 *                     → ADI Y → half-TPA → half-Kerr → half-linear
 */
void propagation_step(
    complex_t phi_in[DIM][DIM],
    complex_t phi_out[DIM][DIM]
) {
    #pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    #pragma HLS INTERFACE m_axi     port=phi_in  offset=slave bundle=G_MEM
    #pragma HLS INTERFACE m_axi     port=phi_out offset=slave bundle=G_MEM
    #pragma HLS bind_storage variable=phi_in  type=ram_1p impl=uram
    #pragma HLS bind_storage variable=phi_out type=ram_1p impl=uram
    #pragma HLS DATAFLOW

    static complex_t tmp1[DIM][DIM], tmp2[DIM][DIM];
    static complex_t tmp3[DIM][DIM], tmp4[DIM][DIM];
    #pragma HLS bind_storage variable=tmp1 type=ram_1p impl=uram
    #pragma HLS bind_storage variable=tmp2 type=ram_1p impl=uram
    #pragma HLS bind_storage variable=tmp3 type=ram_1p impl=uram
    #pragma HLS bind_storage variable=tmp4 type=ram_1p impl=uram
    #pragma HLS ARRAY_PARTITION variable=tmp1 complete dim=2
    #pragma HLS ARRAY_PARTITION variable=tmp2 complete dim=2
    #pragma HLS ARRAY_PARTITION variable=tmp3 complete dim=2
    #pragma HLS ARRAY_PARTITION variable=tmp4 complete dim=2

    // ADI X
    adi_x(phi_in, tmp1);
    // half-TPA X
    half_2photon_absorption(tmp1[0], tmp2[0]); // unroll manual para cada fila
    // half-Kerr X
    half_nonlinear(tmp2[0], tmp3[0]);
    // half-linear X
    half_linear_absorption(tmp3[0], tmp4[0]);

    // ADI Y
    adi_y(tmp4, tmp1);
    // half-TPA Y
    half_2photon_absorption(tmp1[0], tmp2[0]);
    // half-Kerr Y
    half_nonlinear(tmp2[0], tmp3[0]);
    // half-linear Y
    half_linear_absorption(tmp3[0], phi_out[0]);
}
